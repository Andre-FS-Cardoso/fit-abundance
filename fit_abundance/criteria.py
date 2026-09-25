import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from .calibrators import CALIBRATORS

def points(
    name,
    criterion,
    r,
    OH,
    OH_err,
    EWHa,
    Ha6562_cor,
    OIII5006_cor,
    NII6583_cor,
    calibrator,
    save_oh_criteria,
    save_BPT,
    show_graph_bpt
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
    save_oh_criteria : bool
        If True, saves the abundances for a given criterion, along with their respective errors, to a CSV file.
    save_BPT : bool
        If True, save the BPT diagram plot.
    show_graph_bpt:
        If True, displays the BPT diagram plot.

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
    
    # --------------------------- #
    # --- Data to BPT diagram --- #
    # --------------------------- #
    mask_bpt = (
        np.isfinite(EWHa) &
        np.isfinite(Ha6562_cor) &
        np.isfinite(OIII5006_cor) &
        np.isfinite(NII6583_cor)
    )
    x_bpt = NII6583_cor[mask_bpt] - Ha6562_cor[mask_bpt]
    y_bpt = OIII5006_cor[mask_bpt]
    ew_bpt = EWHa[mask_bpt]
    # --------------------------- #

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

        calib_dict = CALIBRATORS

        if calibrator not in calib_dict:
            raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R, 16=NO_F22.")

        calib = CALIBRATORS[calibrator]

        os.makedirs("abundances", exist_ok=True)

        filepath = os.path.join("abundances", f"{name}_OH_{calib}_{criterion}.csv")

        results.to_csv(filepath, index=False)
        
        if save_BPT:
        
            os.makedirs("graphs", exist_ok=True)
            filepath = os.path.join("graphs", f"BPT_{name}.png")
            if not os.path.exists(filepath):
        
                fig, ax = plt.subplots(figsize=(6.2, 5.2))

                X_KE01 = np.linspace(-1.7, 0.4, 100)
                Y_KE01 = 0.61/(X_KE01-0.47) + 1.19
                ax.plot(X_KE01, Y_KE01, c='black', lw=1., ls=':', label = 'KE01')

                X_KA03 = np.linspace(-1.7, 0.0, 100)
                Y_KA03 = 0.61/(X_KA03-0.05) + 1.3
                ax.plot(X_KA03, Y_KA03, c='black', lw=1., ls='--', label = 'KA03')

                X_ST06 = np.linspace(-1.7, -0.25, 100)
                Y_ST06 = (-30.787 + 1.1358*X_ST06 + 0.27297*X_ST06**2)*np.tanh(5.7409*X_ST06) - 31.093
                ax.plot(X_ST06, Y_ST06, c='black', lw=1., ls='-.', label = 'ST06')

                X_EP20 = np.linspace(-1.7, -0.05, 100)
                Y_EP20 = 0.13/(X_EP20 - 0.003) + 0.57
                ax.plot(X_EP20, Y_EP20, c='black', lw=1., ls='-', label = 'EP20')

                plot1 = ax.scatter(x_bpt, y_bpt, s=25.0, marker="s", c=ew_bpt, cmap='Spectral', alpha=0.7)
                cbar = plt.colorbar(plot1, pad=0.01)
                cbar.set_label(label=r'EW(H$\alpha$) [$\AA$]', fontsize=12)
                cbar.ax.tick_params(labelsize=12)
                plot1.set_clim(3.0, 13.0)

                ax.set_xticks(np.arange(-1.5, 0.6, step=0.5), ['-1.5', '-1.0', '-0.5', '0.0', '0.5'], fontsize=12)
                ax.set_xlabel(r'log([NII]$\lambda$6583/H$\alpha$)', fontsize=12)

                ax.set_yticks(np.arange(-1.5, 1.5, step=0.5), ['-1.5', '-1.0', '-0.5', '0.0', '0.5', '1.0'], fontsize=12)
                ax.set_ylabel(r'log([OIII]$\lambda$5007/H$\beta$)', fontsize=12)

                ax.text(-1.6, 1.25, f'{name}', fontsize=16)

                ax.set_xlim([-1.7, 0.7])
                ax.set_ylim([-1.8, 1.5])

                ax.minorticks_on()
                ax.tick_params(which='major', direction='in', length=4.0, width=0.7, colors='black',
                               grid_color='gray', grid_alpha=0.9)
                ax.tick_params(which='minor', direction='in', length=2.0, width=0.5, colors='black',
                               grid_color='gray', grid_alpha=0.9)

                ax.legend(loc=1, ncols=2, fontsize=12)
                
                plt.savefig(filepath, transparent=False, facecolor='w', edgecolor='w')
                
        # Show the plot or not
        if show_graph_bpt:
            plt.show(block=True)
        else:
            plt.close(fig)
                    

    return x, y, yerr
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

