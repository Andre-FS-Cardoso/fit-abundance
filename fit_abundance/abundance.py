import numpy as np
import pandas as pd
import os

def abundance(name,
              r,
              HIIREGID,
              EWHa,
              Hb4861,
              eHb4861,
              Ha6562,
              eHa6562,
              OIII5006,
              eOIII5006,
              NII6583,
              eNII6583,
              SII6716,
              eSII6716,
              SII6730,
              eSII6730,
              calibrator,
              save_abundance
              ):
    """
    Compute oxygen abundances for HII regions using strong-line calibrators.

    The function performs extinction correction of emission-line fluxes,
    computes abundance-sensitive indices, and derives oxygen abundances
    using several widely used calibrators.

    Implemented calibrators
    -----------------------
    1 : PP04 (O3N2)
    2 : PP04 (N2)
    3 : M13 (O3N2)
    4 : M13 (N2)
    5 : D16

    Parameters
    ----------
    name : str
        Name of the galaxy.
    r : array-like
        Galactocentric distance normalized by the effective radius.
    HIIREGID : array-like
        Identifier of each HII region.
    EWHa : array-like
        Equivalent width of Hα (Å).

    Hb4861, Ha6562, OIII5006, NII6583, SII6716, SII6730 : array-like
        Emission-line fluxes.

    eHb4861, eHa6562, eOIII5006, eNII6583, eSII6716, eSII6730 : array-like
        Flux uncertainties.

    calibrator : int
        Identifier of the abundance calibrator.

    Returns
    -------
    oh : numpy.ndarray
        Oxygen abundance (12 + log(O/H)).
    eoh : numpy.ndarray
        Uncertainty in oxygen abundance.
    Ha6562_cor, OIII5006_cor, NII6583_cor : numpy.ndarray
        Extinction-corrected emission lines.
    """
    
    # Convert inputs to numpy arrays
    r = np.asarray(r)
    HIIREGID = np.asarray(HIIREGID)
    EWHa = np.asarray(EWHa)
    
    Hb4861 = np.asarray(Hb4861)
    Ha6562 = np.asarray(Ha6562)
    OIII5006 = np.asarray(OIII5006)
    NII6583 = np.asarray(NII6583)
    SII6716 = np.asarray(SII6716)
    SII6730 = np.asarray(SII6730)
    
    eHb4861 = np.asarray(eHb4861)
    eHa6562 = np.asarray(eHa6562)
    eOIII5006 = np.asarray(eOIII5006)
    eNII6583 = np.asarray(eNII6583)
    eSII6716 = np.asarray(eSII6716)
    eSII6730 = np.asarray(eSII6730)    

    # --- Extinction factor (Cavichia et al. 2010)
    def extinction(x):
        return (0.00001 + 0.22707/x + 1.95243/x**2 - 2.67596/x**3 +
                2.6507/x**4 - 1.26812/x**5 + 0.27549/x**6 - 0.02212/x**7)
    
    av_Hb4861 = extinction(4861.32e-4)
    av_Ha6562 = extinction(6562.68e-4)
    av_OIII5006 = extinction(5006.84e-4)
    av_NII6583 = extinction(6583.41e-4)
    av_SII6716 = extinction(6716.39e-4)
    av_SII6730 = extinction(6730.74e-4)

    # --- Color excess (Cavichia 2008)
    excess = (np.log10(2.86) - np.log10(Ha6562/Hb4861)) / (0.4 * (av_Ha6562 - av_Hb4861))

    # --- Extinction correction and error propagation
    def flux_cor(flux_line, e_flux_line, flux_Hb, e_flux_Hb, av_line, av_Hb, excess):
    
        flux_corr = np.log10(flux_line / flux_Hb) + 0.4 * excess * (av_line - av_Hb)
        
        e_flux_corr = (1/np.log(10)) * np.sqrt(
            (e_flux_line/flux_line)**2 +
            (e_flux_Hb/flux_Hb)**2
        )
        
        return flux_corr, e_flux_corr

    Ha6562_cor, eHa6562_cor = flux_cor(Ha6562, eHa6562, Hb4861, eHb4861, av_Ha6562, av_Hb4861, excess)
    OIII5006_cor, eOIII5006_cor = flux_cor(OIII5006, eOIII5006, Hb4861, eHb4861, av_OIII5006, av_Hb4861, excess)
    NII6583_cor, eNII6583_cor = flux_cor(NII6583, eNII6583, Hb4861, eHb4861, av_NII6583, av_Hb4861, excess)
    SII6716_cor, eSII6716_cor = flux_cor(SII6716, eSII6716, Hb4861, eHb4861, av_SII6716, av_Hb4861, excess)
    SII6730_cor, eSII6730_cor = flux_cor(SII6730, eSII6730, Hb4861, eHb4861, av_SII6730, av_Hb4861, excess)

    # --- Index
    
    ## O3N2 (Alloin et al. 1979)
    O3N2_index = OIII5006_cor + Ha6562_cor - NII6583_cor
    eO3N2_index = np.sqrt(eOIII5006_cor**2 + eHa6562_cor**2 + eNII6583_cor**2)
    ## N2 (T. Storchi-Bergmann et al. 1994)
    N2_index = NII6583_cor - Ha6562_cor
    eN2_index = np.sqrt(eNII6583_cor**2 + eHa6562_cor**2)

    # --- Calibrators
    
    ## PP04 calibrator with O3N2 (Pettini & Pagel 2004)
    O3N2_PP04 = np.where((-1.0 <= O3N2_index) & (O3N2_index <= 1.9), 8.73 - 0.32*O3N2_index, np.nan)
    eO3N2_PP04 = np.where(np.isfinite(O3N2_PP04), 0.32*eO3N2_index, np.nan)
    
    ## PP04 calibrator with N2 (Pettini & Pagel 2004)
    N2_PP04 = np.where((-2.5 <= N2_index) & (N2_index <= -0.3), 8.90 + 0.57*N2_index, np.nan)
    eN2_PP04 = np.where(np.isfinite(N2_PP04), 0.57*eN2_index, np.nan)
    
    ## M13 calibrator with O3N2 (Marino et al. 2013)
    O3N2_M13 = np.where((-1.1 <= O3N2_index) & (O3N2_index <= 1.7), 8.533 - 0.214*O3N2_index, np.nan)
    eO3N2_M13 = np.where(np.isfinite(O3N2_M13), 0.214*eO3N2_index, np.nan)
    
    ## M13 calibrator with N2 (Marino et al. 2013)
    N2_M13 = np.where((-1.6 <= N2_index) & (N2_index <= -0.2), 8.743 + 0.462*N2_index, np.nan)
    eN2_M13 = np.where(np.isfinite(N2_M13), 0.462*eN2_index, np.nan)
    
    ## D16 calibrator (Dopita et al. 2016)
    # To calculate the flux SII6716,30 = SII6716 + SII6730, it is necessary to revert to the linear flux values 10**(log(F)).
    F1 = 10**SII6716_cor
    F2 = 10**SII6730_cor
    sigma1 = np.log(10) * F1 * eSII6716_cor
    sigma2 = np.log(10) * F2 * eSII6730_cor
    # sum and error propagation
    F_total = F1 + F2
    sigma_total = np.sqrt(sigma1**2 + sigma2**2)
    # Calculate the final flux and error of the sum of the SII6716,30 lines
    SII_total_cor = np.log10(F_total)
    eSII_total_cor = sigma_total / (F_total * np.log(10))
    # Here the actual calculation of the D16 calibrator and error propagation begins
    D16 = 8.77 + (NII6583_cor - SII_total_cor) + 0.264 * (NII6583_cor - Ha6562_cor)
    eD16 = np.sqrt(eNII6583_cor**2 + eSII_total_cor**2 + (0.264**2) * (eNII6583_cor**2 + eHa6562_cor**2))
    
    # --- Quality masks
    with np.errstate(invalid='ignore', divide='ignore'):
        rHb4861 = (eHb4861 < 0.997 * Hb4861) & (Hb4861 > 0)
        rHa6562 = (eHa6562 < 0.997 * Ha6562) & (Ha6562 > 0)
        rOIII5006 = (eOIII5006 < 0.997 * OIII5006) & (OIII5006 > 0)
        rNII6583 = (eNII6583 < 0.997 * NII6583) & (NII6583 > 0)
        rSII6716 = (eSII6716 < 0.997 * SII6716) & (SII6716 > 0)
        rSII6730 = (eSII6730 < 0.997 * SII6730) & (SII6730 > 0)

    # --- Apply quality masks according to the index
    def apply_mask(oh, eoh, mask):
        oh_filtered = np.where(mask, oh, np.nan)
        eoh_filtered = np.where(mask, np.abs(eoh), np.nan)
        return oh_filtered, eoh_filtered

    if calibrator in [1, 3]:  # O3N2 index
        mask_ok = rHb4861 & rHa6562 & rOIII5006 & rNII6583
        if calibrator == 1:
            oh, eoh = apply_mask(O3N2_PP04, eO3N2_PP04, mask_ok)
        else:
            oh, eoh = apply_mask(O3N2_M13, eO3N2_M13, mask_ok)
    elif calibrator in [2, 4]:  # N2 index
        mask_ok = rHa6562 & rNII6583
        if calibrator == 2:
            oh, eoh = apply_mask(N2_PP04, eN2_PP04, mask_ok)
        else:
            oh, eoh = apply_mask(N2_M13, eN2_M13, mask_ok)
    elif calibrator == 5:  # D16
        mask_ok = rHa6562 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(D16, eD16, mask_ok)
    else:
        raise ValueError("Invalid calibrator. Use 1=O3N2_PP04, 2=N2_PP04, 3=O3N2_M13, 4=N2_M13, 5=D16.")

    # --- Salve CSV
    results = pd.DataFrame({
        'HIIREGID': HIIREGID,
        'r': r,
        'EWHa6562': EWHa,
        'Ha6562_cor': Ha6562_cor,
        'eHa6562_cor': eHa6562_cor,
        'OIII5006_cor': OIII5006_cor,
        'eOIII5006_cor': eOIII5006_cor,
        'NII6583_cor': NII6583_cor,
        'eNII6583_cor': eNII6583_cor,
        'SII6716_cor': SII6716_cor,
        'eSII6716_cor': eSII6716_cor,
        'SII6730_cor': SII6730_cor,
        'eSII6730_cor': eSII6730_cor,
        'OH_PP04_O3N2': O3N2_PP04,
        'eOH_PP04_O3N2': eO3N2_PP04,
        'OH_PP04_N2': N2_PP04,
        'eOH_PP04_N2': eN2_PP04,
        'OH_M13_O3N2': O3N2_M13,
        'eOH_M13_O3N2': eO3N2_M13,
        'OH_M13_N2': N2_M13,
        'eOH_M13_N2': eN2_M13,
        'OH_D16': D16,
        'eOH_D16': eD16
    })
    
    
    if save_abundance:
    
        os.makedirs("abundances", exist_ok=True)
        
        # --- Save the corrected flux in a table
        filepath = os.path.join("abundances", f"{name}_abundances.csv")

        results.to_csv(filepath, index=False)

    return oh, eoh, Ha6562_cor, OIII5006_cor, NII6583_cor

