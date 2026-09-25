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
    n_break=2
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
    n_break : int
        Number of breaks to be fitted. By default, n_break = 2.

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
    
    # CASE 2, 3 or 4 with 1, 2 or 3 breakpoints
    if n_break not in [2, 3]:
        raise ValueError("Default = 2 breaks. You can choose only 3 breaks.")
    else: 
        fits = {}
        rss_values = []

        max_breaks = n_break

        for n in range(1, max_breaks + 1):

            fit = piecewise_regression.main.Fit(
                x,
                y,
                n_boot=n_boot,
                n_breakpoints=n,
                min_distance_to_edge=0.05,
                min_distance_between_breakpoints=0.20
            )

            fits[n] = fit

            if fit.get_results()["converged"]:
                rss_values.append(fit.get_results()["rss"])
            else:
                rss_values.append(np.inf)

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
    AIC2 = aic_final(x, rss_values[0], 4)
    AIC3 = aic_final(x, rss_values[1], 6)
    
    # Selection of the best model
    if n_break == 2:
        AICs = np.array([AIC1, AIC2, AIC3])
        k = np.array([2,4,6])
    else:
        AIC4 = aic_final(x, rss_values[2], 8)
        AICs = np.array([AIC1, AIC2, AIC3, AIC4])
        k = np.array([2,4,6,8])

    # ---------------------------------------------------------
    # AIC selection
    # ---------------------------------------------------------
    delta = AICs - np.min(AICs)

    weights = np.exp(-0.5 * delta)
    weights /= np.sum(weights)

    # modelos plausíveis
    candidates = np.where(delta <= 4)[0]

    if len(candidates) == 1:
        best_model = candidates[0] + 1
    else:
        # escolhe o mais simples
        best_model = candidates[np.argmin(k[candidates])] + 1
        
    if save_model_selection:
    
        calib_dict = CALIBRATORS

        if calibrator not in calib_dict:
            raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R, 16=NO_F22.")

        calib = CALIBRATORS[calibrator]

        os.makedirs("model_selection", exist_ok=True)
        
        if n_break == 2:
            model_names = ["Linear", "1-break", "2-break"]
        else:
            model_names = ["Linear", "1-break", "2-break", "3-break"]

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

    if n_break == 2:
        return {
            'x': x,
            'y': y,
            'ey': ey,
            'best_case': best_model,
            'fit1': (a2, ea2, b0, eb0, rss_linear),
            'fit2': fits[1],
            'fit3': fits[2],
            'AICs': AICs,
            'delta_AIC': delta,
            'weights': weights
        }
    else:
        return {
            'x': x,
            'y': y,
            'ey': ey,
            'best_case': best_model,
            'fit1': (a2, ea2, b0, eb0, rss_linear),
            'fit2': fits[1],
            'fit3': fits[2],
            'fit4': fits[3],
            'AICs': AICs,
            'delta_AIC': delta,
            'weights': weights
        }
