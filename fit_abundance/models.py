import statsmodels.api as sm
import piecewise_regression
import numpy as np
import pandas as pd
import os
from .calibrators import CALIBRATORS

def fit_models(
    x_array,
    y_array,
    ey_array,
    name,
    criterion,
    calibrator,
    save_model_selection,
    *,
    n_boot=200,
    min_distance_to_edge=0.1,
    min_distance_between_breakpoints=0.15,
    start_values=None
    ):
    
    """
    Fit oxygen and nitrogen abundance gradients using linear and piecewise models.

    Three models are tested:

        1. Simple linear gradient
        2. Broken gradient with one breakpoint
        3. Broken gradient with two breakpoints

    The best model is selected using the Akaike Information Criterion (AIC).
    For small samples the corrected AIC (AICc) is used.
    
    The "average" module combines the parameters from the different fits based on the
    weight of each AIC fit.

    Parameters
    ----------
    x_array : array-like
        Galactocentric radius (normalized by effective radius).
    y_array : array-like
        Oxygen abundance (12 + log(O/H)).
    ey_array : array-like
        Abundance uncertainty.
    name : str
        Name of the galaxy.
    criterion : str
        Identifier of the selection criterion.
    calibrator : int
        Abundance calibrator identifier.
    save_model_selection : bool
        If True, saves the table containing the statistical information of the model fits to a CSV file.
    n_boot : int
        Number of bootstrap iterations to be performed. By default, n_boot = 200.
    min_distance_to_edge : float
        Minimum distance from edge of data to a breakpoint, as a proportion of the data range
    min_distance_between_breakpoints : float
        The minimum required distance between breakpoints, as a proportion of the data range.
    start_values : List(float)
        Initial guesses for breakpoint positions.

    Returns
    -------
    dict
        Dictionary containing the fitted models and the best model index.
    """

    x = np.asarray(x_array)
    y = np.asarray(y_array)
    ey = np.asarray(ey_array)
    
    # --- REMOVE NaNs ---
    mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(ey)
    
    x = x[mask]
    y = y[mask]
    ey = ey[mask]
    
    if len(x) < 10:
        print("Warning: fewer than 10 data points available for fitting.")
        return None
        
    # CASE 1 fit: simple linear regression
    X = sm.add_constant(x)
    model = sm.OLS(y, X)
    results = model.fit()
    a_1 = results.params[1]
    ea_1 = results.bse[1]
    b_1 = results.params[0]
    eb_1 = results.bse[0]
    rss_linear = results.ssr
    
    if start_values is None:
        start_values = [0.5, 1.5]
    
    # CASE 2 fit: 1 breakpoint
    fit2 = piecewise_regression.main.Fit(x, y, n_boot=n_boot, n_breakpoints=1,
                                         min_distance_to_edge=min_distance_to_edge)
    rss_1break = fit2.get_results()["rss"] if fit2.get_results()["converged"] else np.inf
    
    # CASE 3 fit: 2 breakpoints
    fit3 = piecewise_regression.main.Fit(x, y, n_boot=n_boot, n_breakpoints=2,
                                         min_distance_to_edge=min_distance_to_edge, 
                                         min_distance_between_breakpoints=min_distance_between_breakpoints,
                                         start_values=start_values)
    rss_2break = fit3.get_results()["rss"] if fit3.get_results()["converged"] else np.inf

    # Functions for AIC
    def llf_(X, rss):
        nobs = float(X.shape[0])
        if rss == np.inf or rss <= 0:
            return -np.inf
        llf = -0.5 * nobs * (np.log(2*np.pi) + np.log(rss/nobs) + 1)
        return llf
        
    def aic_final(X, rss, k):
        nobs = float(X.shape[0])
        llf = llf_(X, rss)
        
        if llf == -np.inf:
            return np.inf
        
        aic = -2*llf + 2*k
        
        if (nobs / k) < 40:
            aic += (2 * k * (k + 1)) / (nobs - k - 1)
        
        return aic
        
    # AIC calculation 
    AIC1 = aic_final(x, rss_linear, 2)
    AIC2 = aic_final(x, rss_1break, 4)
    AIC3 = aic_final(x, rss_2break, 6)
    
    AICs = np.array([AIC1, AIC2, AIC3])
    k = np.array([2,4,6])
    
    valid_aics = AICs[np.isfinite(AICs)]
    if len(valid_aics) == 0:
        print("Warning: All models failed to yield a valid AIC.")
        return None

    # ---------------------------------------------------------
    # AIC selection
    # ---------------------------------------------------------
    min_aic = np.min(valid_aics)

    # Delta calculated for valid values
    delta = np.where(np.isfinite(AICs), AICs - min_aic, np.inf)

    # Akaike weights
    raw_weights = np.where(np.isfinite(delta), np.exp(-0.5 * delta), 0.0)
    weights = raw_weights / np.sum(raw_weights)
  
    w_linear = weights[0]               
    w_nonlinear = np.sum(weights[1:])  

    # Classification
    if w_nonlinear < 0.30:
        profile = "Linear"
        description = "Stable and regular disk evolution following a classical inside-out growth scenario."
    elif 0.30 <= w_nonlinear <= 0.70:
        profile = "Transition"
        description = "Early signature of gas reallocation, likely driven by weak bar dynamics or DIG contamination."
    else: 
        profile = "Broken"
        description = "Strong physical evidence of non-linear features, such as outer flattening or an inner drop."
        
    # Best model selection
    best = np.argmax(weights)
    ratios = weights[best] / (weights + 1e-15)
    ratios = ratios[np.arange(len(ratios)) != best]

    if np.all(ratios > 2):
        best_model = best + 1
    else:
        candidates = np.where(delta <= 2)[0]
        best_model = candidates[np.argmin(k[candidates])] + 1
        
    # --- MODEL AVERAGING ---
    def predict_combined(x_input):
        # 1. Prediction for linear model (Statsmodels)
        X_input = sm.add_constant(x_input, has_constant='add')
        pred_1 = results.predict(X_input)
        
        # 2. Prediction for 1 breakpoint
        if fit2.get_results()["converged"]:
            est2 = fit2.get_results()["estimates"]
            a1 = est2["alpha1"]["estimate"]
            b1 = est2["const"]["estimate"]
            a2 = est2["alpha2"]["estimate"]
            x1 = est2["breakpoint1"]["estimate"]
            y1 = a1 * x1 + b1
            pred_2 = np.where(x_input <= x1, 
                              a1 * x_input + b1, 
                              a2 * (x_input - x1) + y1)
        else:
            pred_2 = np.zeros_like(x_input)
        
        # 3. Prediction for 2 breakpoints
        if fit3.get_results()["converged"]:
            est3 = fit3.get_results()["estimates"]
            a1 = est3["alpha1"]["estimate"]
            b1 = est3["const"]["estimate"]
            a2 = est3["alpha2"]["estimate"]
            a3 = est3["alpha3"]["estimate"]
            x1 = est3["breakpoint1"]["estimate"]
            x2 = est3["breakpoint2"]["estimate"]
            
            y1 = a1 * x1 + b1
            y2 = a2 * (x2 - x1) + y1
            
            condicoes = [
                x_input <= x1,
                (x_input > x1) & (x_input <= x2),
                x_input > x2
            ]
            valores = [
                a1 * x_input + b1,
                a2 * (x_input - x1) + y1,
                a3 * (x_input - x2) + y2
            ]
            pred_3 = np.select(condicoes, valores)
        else:
            pred_3 = np.zeros_like(x_input)
        
        # Weighted average using Akaike weights
        y_combined = (weights[0] * pred_1) + (weights[1] * pred_2) + (weights[2] * pred_3)
        return y_combined
        
    # Predictions for the observed x vector
    y_fitted_combined = predict_combined(x)


    if save_model_selection:
    
        calib_dict = CALIBRATORS

        if calibrator not in calib_dict:
            raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R or 16=NO_F22.")

        calib = CALIBRATORS[calibrator]

        os.makedirs("model_selection", exist_ok=True)
        
        model_names = ["Linear", "1-break", "2-break"]

        df = pd.DataFrame({
            "model": model_names,
            "AIC": AICs,
            "delta_AIC": delta,
            "weight": weights
        })
        
        df["selected"] = False
        df.loc[best_model-1, "selected"] = True

        filepath = os.path.join("model_selection", f"{name}_AIC_{calib}_{criterion}.csv")
        df.to_csv(filepath, index=False)

    return {
        'x': x,
        'y': y,
        'ey': ey,
        'best_case': best_model,
        'fit1_results': results,                       
        'fit1': (a_1, ea_1, b_1, eb_1, rss_linear),
        'fit2': fit2,
        'fit3': fit3,
        'AICs': AICs,
        'delta_AIC': delta,
        'weights': weights,
        'y_fitted_combined': y_fitted_combined,         
        'predict_combined_func': predict_combined,
        'w_nonlinear': w_nonlinear,
        'profile': profile,
        'description': description   
    }
