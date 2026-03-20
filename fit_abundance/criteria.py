import numpy as np
import pandas as pd
import os

def points(name,
           criterion,
           r,
           OH,
           OH_err,
           EWHa,
           Ha6562_cor,
           OIII5006_cor,
           NII6583_cor,
           calibrator,
           save_oh_criteria
           ):
    """
    Select HII regions based on commonly used diagnostic criteria.

    The function filters regions using emission-line ratios and
    equivalent width thresholds following several prescriptions
    from the literature.

    Implemented criteria
    --------------------
    ST06 : Stasińska et al. (2006)
    KA03 : Kauffmann et al. (2003)
    KE01 : Kewley et al. (2001)
    KE6A : Kewley et al. (2001) with EW(Hα) ≥ 6 Å
    CF11 : Cid Fernandes et al. (2011)
    EP20 : Espinosa-Ponce et al. (2020)

    Parameters
    ----------
    name : str
        Name of the galaxy.
    criterion : str
        Identifier of the selection criterion.
    r : array-like
        Galactocentric radius normalized by the effective radius.
    OH : array-like
        Oxygen abundance (12 + log(O/H)).
    OH_err : array-like
        Abundance uncertainty.
    EWHa : array-like
        Equivalent width of Hα (Å).
    Ha6562_cor, OIII5006_cor, NII6583_cor : array-like
        Extinction-corrected emission lines.
    calibrator : int
        Abundance calibrator identifier.
    save_table : bool
        If True, save the filtered dataset to a CSV file.

    Returns
    -------
    x : numpy.ndarray
        Radius of selected regions.
    y : numpy.ndarray
        Oxygen abundance of selected regions.
    yerr : numpy.ndarray
        Abundance uncertainty.
    """
    
    # Convert inputs to numpy arrays
    r = np.asarray(r)
    OH = np.asarray(OH)
    OH_err = np.asarray(OH_err)
    EWHa = np.asarray(EWHa)
    Ha6562_cor = np.asarray(Ha6562_cor)
    OIII5006_cor = np.asarray(OIII5006_cor)
    NII6583_cor = np.asarray(NII6583_cor)

    # --- Mask of valid values
    mask_valid = (
        np.isfinite(OH) &
        np.isfinite(EWHa) &
        np.isfinite(OH_err)
    )

    # --- Apply initial mask to all arrays
    OH = OH[mask_valid]
    OH_err = OH_err[mask_valid]
    EWHa = EWHa[mask_valid]
    Ha6562_cor = Ha6562_cor[mask_valid]
    NII6583_cor = NII6583_cor[mask_valid]
    OIII5006_cor = OIII5006_cor[mask_valid]
    r = r[mask_valid]

    # --- Calculate useful ratios
    log_NII_Ha = NII6583_cor - Ha6562_cor
    log_OIII_Hb = OIII5006_cor
    EW = EWHa

    # --- Specific filters by criterion
    if criterion is None or criterion.lower() == 'none':
        mask = np.ones_like(OH, dtype=bool)
        
    elif criterion == 'ST06':
        mask = (log_NII_Ha <= -0.40) & (log_OIII_Hb <= ((-30.787 + 1.1358 * log_NII_Ha + 0.27297 * log_NII_Ha**2)
                                                 * np.tanh(5.7409 * log_NII_Ha) - 31.093))
    elif criterion == 'KA03':
        mask = (log_OIII_Hb <= (0.61 / (log_NII_Ha - 0.05) + 1.3))
               
    elif criterion == 'KE01':
        mask = (log_OIII_Hb <= (0.61 / (log_NII_Ha - 0.47) + 1.19))
               
    elif criterion == 'KE6A':
        mask = (log_OIII_Hb <= (0.61 / (log_NII_Ha - 0.47) + 1.19)) & \
               (EW >= 6.)
               
    elif criterion == 'CF11':
        mask = (EW >= 3.) & (log_NII_Ha <= -0.4)
        
    elif criterion == 'EP20':
        mask = (log_OIII_Hb <= (0.13 / (log_NII_Ha - 0.003) + 0.57))
    else:
        raise ValueError(f"Criterion '{criterion}' not recognized.")

    # --- Apply final mask
    x = np.array(r[mask])
    y = np.array(OH[mask])
    yerr = np.array(OH_err[mask])
    
    # --- Save output table in CSV (optional)
    results = pd.DataFrame({
        'r': x,
        'OH': y,
        'eOH': yerr
    })
        
    if save_oh_criteria:

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

        os.makedirs("oh_criterions", exist_ok=True)

        filepath = os.path.join("oh_criterions", f"{name}_OH_{criterion}_{calib}.csv")

        results.to_csv(filepath, index=False)

    return x, y, yerr

