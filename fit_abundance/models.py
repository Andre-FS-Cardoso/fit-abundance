import statsmodels.api as sm
import piecewise_regression
import numpy as np
import pandas as pd
import os

def fit_models(x_array,
               y_array,
               ey_array,
               name,
               criterion,
               calibrator,
               save_model_selection,
               n_boot=200
               ):
    """
    Fit oxygen abundance gradients using linear and piecewise models.

    Three models are tested:

    1. Simple linear gradient
    2. Broken gradient with one breakpoint
    3. Broken gradient with two breakpoints

    The best model is selected using the Akaike Information Criterion (AIC).
    For small samples the corrected AIC (AICc) is used.

    Parameters
    ----------
    x_array : array-like
        Galactocentric radius (normalized by effective radius).
    y_array : array-like
        Oxygen abundance (12 + log(O/H)).
    ey_array : array-like
        Abundance uncertainty.

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
    a2 = results.params[1]
    ea2 = results.bse[1]
    b0 = results.params[0]
    eb0 = results.bse[0]
    rss_linear = results.ssr
    
    # CASE 2 fit: 1 breakpoint
    fit2 = piecewise_regression.main.Fit(x, y, n_boot=n_boot, n_breakpoints=1, min_distance_to_edge=0.05)
    rss_1break = fit2.get_results()["rss"] if fit2.get_results()["converged"] else np.inf
    
    # CASE 3 fit: 2 breakpoints
    fit3 = piecewise_regression.main.Fit(x, y, n_boot=n_boot, n_breakpoints=2, min_distance_to_edge=0.05, 
                                         min_distance_between_breakpoints=0.20, start_values=[0.5, 1.5])
    rss_2break = fit3.get_results()["rss"] if fit3.get_results()["converged"] else np.inf

    # Functions for AIC
    def llf_(X, rss):
        nobs = float(X.shape[0])
        
        llf = -0.5 * nobs * (np.log(2*np.pi) + np.log(rss/nobs) + 1)
        return llf
        
    def aic_final(X, rss, k):
        nobs = float(X.shape[0])
        
        llf = llf_(X, rss)
        
        aic = -2*llf + 2*k
        
        if (nobs / k) < 40:
            aic += (2 * k * (k + 1)) / (nobs - k - 1)
        
        return aic
        
    AIC1 = aic_final(x, rss_linear, 2)
    AIC2 = aic_final(x, rss_1break, 4)
    AIC3 = aic_final(x, rss_2break, 6)

    # Selection of the best model
    AICs = np.array([AIC1, AIC2, AIC3])
    k = np.array([2,4,6])

    delta = AICs - np.min(AICs)

    weights = np.exp(-0.5*delta)
    weights /= np.sum(weights)

    best = np.argmax(weights)

    ratios = weights[best] / weights

    # ignorar o próprio modelo
    ratios = ratios[np.arange(len(ratios)) != best]

    if np.all(ratios > 2):
        best_model = best + 1
    else:
        # escolher mais simples entre os plausíveis
        candidates = np.where(delta <= 2)[0]
        best_model = candidates[np.argmin(k[candidates])] + 1
        
    if save_model_selection:
    
        calib_dict = {
            1: "PP04_O3N2",
            2: "PP04_N2",
            3: "M13_O3N2",
            4: "M13_N2",
            5: "D16"
        }

        if calibrator not in calib_dict:
            raise ValueError("Invalid calibrator. Use 1=O3N2_PP04, 2=N2_PP04, 3=O3N2_M13, 4=N2_M13, 5=D16.")

        calib = calib_dict[calibrator]

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

        filepath = os.path.join("model_selection", f"{name}_AIC_{criterion}_{calib}.csv")

        df.to_csv(filepath, index=False)

    return {
        'x': x,
        'y': y,
        'ey': ey,
        'best_case': best_model,
        'fit1': (a2, ea2, b0, eb0, rss_linear),
        'fit2': fit2,
        'fit3': fit3,
        'AICs': AICs,
        'delta_AIC': delta,
        'weights': weights
    }

