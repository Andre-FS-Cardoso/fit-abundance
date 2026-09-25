import numpy as np
import pandas as pd
import os

def abundance(
    name,
    r,
    HIIREGID,
    EWHa,
    Hb4861,
    eHb4861,
    Ha6562,
    eHa6562,
    OII3727,
    eOII3727,
    OIII4958,
    eOIII4958,
    OIII5006,
    eOIII5006,
    NII6548,
    eNII6548,
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
    Compute oxygen and nitrogen abundances for HII regions using strong-line calibrators.

    The function performs extinction correction of emission-line fluxes,
    computes abundance-sensitive indices, and derives oxygen and nitrogen
    abundances using several widely used calibrators.

    Implemented calibrators
    -----------------------
      Oxygen abundances:
    
         1 : PP04 - O3N2        Pettini & Pagel 2004
         2 : PP04 - N2          Pettini & Pagel 2004
         3 : PP04 - N2_poly     Pettini & Pagel 2004
         4 : M13  - O3N2        Marino et al. 2013
         5 : M13  - N2          Marino et al. 2013
         6 : D16                Dopita et al. 2016
         7 : T04  - R23         Tremonti et al. 2004
         8 : KD02 - N2O2        Kewley and Dopita 2002
         9 : P10  - ONS         Pilyugin et al. 2010
        10 : P10  - ON          Pilyugin et al. 2010
        11 : PM11 - NS          Pilyugin and Mattsson 2011
        12 : PG16 - R           Pilyugin and Grebel 2016
        13 : PG16 - S           Pilyugin and Grebel 2016
    
      Nitrogen abundance:
    
        14 : NH_PG16 - R        Pilyugin and Grebel 2016
        15 : NO_PG16 - R        Pilyugin and Grebel 2016
        16 : NO_F22  - N2O2     Florido et al. 2022

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
    Hb4861, Ha6562, OII3727, OIII4958, OIII5006, NII6548, NII6583, SII6716, SII6730 : array-like
        Emission-line fluxes.
    eHb4861, eHa6562m eOII3727, eOIII4958, eOIII5006, eNII6548, eNII6583, eSII6716, eSII6730 : array-like
        Flux uncertainties.
    calibrator : int
        Identifier of the abundance calibrator.
    save_abundance : bool
        If True, saves all corrected fluxes and abundances, along with their respective errors, to a CSV file.

    Returns
    -------
    oh : numpy.ndarray
        Oxygen and nitrogen abundance.
    eoh : numpy.ndarray
        Uncertainty in oxygen and nitrogen abundance.
    Ha6562_cor, OIII5006_cor, NII6583_cor : numpy.ndarray
        Extinction-corrected emission lines.
    """
    
    # Convert inputs to numpy arrays
    r = np.asarray(r)
    HIIREGID = np.asarray(HIIREGID)
    EWHa = np.asarray(EWHa)
    
    Hb4861 = np.asarray(Hb4861)
    Ha6562 = np.asarray(Ha6562)
    OII3727 = np.asarray(OII3727)
    OIII4958 = np.asarray(OIII4958)
    OIII5006 = np.asarray(OIII5006)
    NII6548 = np.asarray(NII6548)
    NII6583 = np.asarray(NII6583)
    SII6716 = np.asarray(SII6716)
    SII6730 = np.asarray(SII6730)
    
    eHb4861 = np.asarray(eHb4861)
    eHa6562 = np.asarray(eHa6562)
    eOII3727 = np.asarray(eOII3727)
    eOIII4958 = np.asarray(eOIII4958)
    eOIII5006 = np.asarray(eOIII5006)
    eNII6548 = np.asarray(eNII6548)
    eNII6583 = np.asarray(eNII6583)
    eSII6716 = np.asarray(eSII6716)
    eSII6730 = np.asarray(eSII6730)

    def propagate_error(func, values, errors, rel_step=1e-6):
        """
        Generic error propagation using numerical derivatives (finite differences).

        func: function f(x1, x2, ..., xn) that returns an array
        values: list of arrays [x1, x2, ..., xn]
        errors: list of arrays [ex1, ex2, ..., exn]
        rel_step: relative step for numerical differentiation

        Returns: error (array)
        """
        values = [np.asarray(v) for v in values]
        errors = [np.asarray(e) for e in errors]

        f0 = func(*values)
        sigma2 = np.zeros_like(f0, dtype=float)

        valid = np.isfinite(f0)

        for i, (xi, si) in enumerate(zip(values, errors)):

            step = rel_step * (np.abs(xi) + 1e-12)

            x_plus  = values.copy()
            x_minus = values.copy()

            x_plus[i]  = xi + step
            x_minus[i] = xi - step

            f_plus  = func(*x_plus)
            f_minus = func(*x_minus)

            dfdxi = (f_plus - f_minus) / (2.0 * step)

            term = (dfdxi * si)**2
            term[~valid] = 0.0

            sigma2 += term

        sigma = np.sqrt(sigma2)
        sigma[~valid] = np.nan

        return sigma

    # --- Extinction factor (Cavichia et al. 2010)
    def extinction(x):
        return (0.00001 + 0.22707/x + 1.95243/x**2 - 2.67596/x**3 +
                2.6507/x**4 - 1.26812/x**5 + 0.27549/x**6 - 0.02212/x**7)
    
    av_Hb4861 = extinction(4861.32e-4)
    av_Ha6562 = extinction(6562.68e-4)
    av_OII3727 = extinction(3727.40e-4)
    av_OIII4958 = extinction(4958.91e-4)
    av_OIII5006 = extinction(5006.84e-4)
    av_NII6548 = extinction(6548.08e-4)
    av_NII6583 = extinction(6583.41e-4)
    av_SII6716 = extinction(6716.39e-4)
    av_SII6730 = extinction(6730.74e-4)

    # --- Color excess (Cavichia 2008)
    def f_excess(Ha, Hb):
        return (np.log10(2.86) - np.log10(Ha/Hb)) / (0.4 * (av_Ha6562 - av_Hb4861))
    excess = f_excess(Ha6562, Hb4861)
    e_excess = propagate_error(f_excess, [Ha6562, Hb4861], [eHa6562, eHb4861])

    # --- Extinction correction and error propagation
    def f_flux(flux_line, flux_Hb, excess, av_line):
        return np.log10(flux_line / flux_Hb) + 0.4 * excess * (av_line - av_Hb4861)
        
    Ha6562_cor = f_flux(Ha6562, Hb4861, excess, av_Ha6562)
    eHa6562_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_Ha6562), [Ha6562, Hb4861, excess], [eHa6562, eHb4861, e_excess])
    
    OII3727_cor = f_flux(OII3727, Hb4861, excess, av_OII3727)
    eOII3727_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_OII3727), [OII3727, Hb4861, excess], [eOII3727, eHb4861, e_excess])
    
    OIII4958_cor = f_flux(OIII4958, Hb4861, excess, av_OIII4958)
    eOIII4958_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_OIII4958), [OIII4958, Hb4861, excess], [eOIII4958, eHb4861, e_excess])
    
    OIII5006_cor = f_flux(OIII5006, Hb4861, excess, av_OIII5006)
    eOIII5006_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_OIII5006), [OIII5006, Hb4861, excess], [eOIII5006, eHb4861, e_excess])
    
    NII6548_cor = f_flux(NII6548, Hb4861, excess, av_NII6548)
    eNII6548_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_NII6548), [NII6548, Hb4861, excess], [eNII6548, eHb4861, e_excess])
    
    NII6583_cor = f_flux(NII6583, Hb4861, excess, av_NII6583)
    eNII6583_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_NII6583), [NII6583, Hb4861, excess], [eNII6583, eHb4861, e_excess])
    
    SII6716_cor = f_flux(SII6716, Hb4861, excess, av_SII6716)
    eSII6716_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_SII6716), [SII6716, Hb4861, excess], [eSII6716, eHb4861, e_excess])
    
    SII6730_cor = f_flux(SII6730, Hb4861, excess, av_SII6730)
    eSII6730_cor= propagate_error(lambda fl, fhb, ex: f_flux(fl, fhb, ex, av_SII6730), [SII6730, Hb4861, excess], [eSII6730, eHb4861, e_excess])

    ## ============================================================================================================================================================================ ##
    # --- Index
    ## ============================================================================================================================================================================ ##
    
    ## O3N2 (Alloin et al. 1979)
    O3N2_index = OIII5006_cor + Ha6562_cor - NII6583_cor
    eO3N2_index = np.sqrt(eOIII5006_cor**2 + eHa6562_cor**2 + eNII6583_cor**2)
    
    ## N2 (T. Storchi-Bergmann et al. 1994)
    N2_index = NII6583_cor - Ha6562_cor
    eN2_index = np.sqrt(eHa6562_cor**2 + eNII6583_cor**2)
    
    ## N2O2 (Kewley & Dopita 2002)
    N2O2 = NII6583_cor - OII3727_cor
    eN2O2 = np.sqrt(eNII6583_cor**2 + eOII3727_cor**2)
    
    # R2 index (Pilyugin et al. 2010)
    R2_total = 10**OII3727_cor
    sigma1_R2 = np.log(10) * 10**OII3727_cor * eOII3727_cor
    sigma_R2_total = np.sqrt(sigma1_R2**2)
    R2 = np.log10(R2_total)
    eR2 = sigma_R2_total / (R2_total * np.log(10))
    
    # N2 index (Pilyugin et al. 2010)
    N2_total = 10**NII6548_cor + 10**NII6583_cor
    sigma1_N2 = np.log(10) * 10**NII6548_cor * eNII6548_cor
    sigma2_N2 = np.log(10) * 10**NII6583_cor * eNII6583_cor
    sigma_N2_total = np.sqrt(sigma1_N2**2 + sigma2_N2**2)
    N2 = np.log10(N2_total)
    eN2 = sigma_N2_total / (N2_total * np.log(10))
    
    # S2 index (Pilyugin et al. 2010)
    S2_total = 10**SII6716_cor + 10**SII6730_cor
    sigma1_S2 = np.log(10) * 10**SII6716_cor * eSII6716_cor
    sigma2_S2 = np.log(10) * 10**SII6730_cor * eSII6730_cor
    sigma_S2_total = np.sqrt(sigma1_S2**2 + sigma2_S2**2)
    S2 = np.log10(S2_total)
    eS2 = sigma_S2_total / (S2_total * np.log(10))
    
    # R3 index (Pilyugin et al. 2010)
    R3_total = 10**OIII4958_cor + 10**OIII5006_cor
    sigma1_R3 = np.log(10) * 10**OIII4958_cor * eOIII4958_cor
    sigma2_R3 = np.log(10) * 10**OIII5006_cor * eOIII5006_cor
    sigma_R3_total = np.sqrt(sigma1_R3**2 + sigma2_R3**2)
    R3 = np.log10(R3_total)
    eR3 = sigma_R3_total / (R3_total * np.log(10))
    
    # P index (Pilyugin et al. 2010)
    def f_P(R3, R2):
        return R3 / (R3 + R2)
    P = f_P(R3_total, R2_total)
    eP = propagate_error(f_P, [R3_total, R2_total], [sigma_R3_total, sigma_R2_total])
    
    # R23 index (Pilyugin and Mattsson 2011)
    R23 = np.log10(R3_total + R2_total)
    eR23 = np.sqrt(sigma_R3_total**2 + sigma_R2_total**2) / (R23 * np.log(10))
    
    ## ============================================================================================================================================================================ ##
    # --- Calibrators
    ## ============================================================================================================================================================================ ##
    
    ###################################################
    ################ Oxygen Abundance #################
    ###################################################
    
    ## --------------------------------------------- ##
    ## PP04 calibrator with O3N2 (Pettini & Pagel 2004)
    ## --------------------------------------------- ##
    def f_PP04_O3N2(x):
        return 8.73 - 0.32*x
    mask = (-1.0 <= O3N2_index) & (O3N2_index <= 1.9)
    PP04_O3N2 = np.where(mask, f_PP04_O3N2(O3N2_index), np.nan)
    ePP04_O3N2 = np.where(mask, propagate_error(f_PP04_O3N2, [O3N2_index], [eO3N2_index]), np.nan)
    
    ## ------------------------------------------- ##
    ## PP04 calibrator with N2 (Pettini & Pagel 2004)
    ## ------------------------------------------- ##
    def f_PP04_N2(x):
        return 8.90 + 0.57*x
    mask = (-2.5 <= N2_index) & (N2_index <= -0.3)
    PP04_N2 = np.where(mask, f_PP04_N2(N2_index), np.nan)
    ePP04_N2 = np.where(mask, propagate_error(f_PP04_N2, [N2_index], [eN2_index]), np.nan)
    
    ## ------------------------------------------- ##
    ## PP04 calibrator with N2 (Pettini & Pagel 2004)
    ## ------------------------------------------- ##
    def f_PP04_N2_poly(x):
        return 9.37 + 2.03*x + 1.26*x**2 + 0.32*x**3
    mask = (-2.5 <= N2_index) & (N2_index <= -0.3)
    PP04_N2_poly = np.where(mask, f_PP04_N2_poly(N2_index), np.nan)
    ePP04_N2_poly = np.where(mask, propagate_error(f_PP04_N2_poly, [N2_index], [eN2_index]), np.nan)
    
    ## ------------------------------------------ ##
    ## M13 calibrator with O3N2 (Marino et al. 2013)
    ## ---------------------------------------- - ##
    def f_M13_O3N2(x):
        return 8.533 - 0.214*x
    mask = (-1.1 <= O3N2_index) & (O3N2_index <= 1.7)
    M13_O3N2 = np.where(mask, f_M13_O3N2(O3N2_index), np.nan)
    eM13_O3N2 = np.where(mask, propagate_error(f_M13_O3N2, [O3N2_index], [eO3N2_index]), np.nan)
    
    ## ---------------------------------------- ##
    ## M13 calibrator with N2 (Marino et al. 2013)
    ## ---------------------------------------- ##
    def f_M13_N2(x):
        return 8.743 + 0.462*x
    mask = (-1.6 <= N2_index) & (N2_index <= -0.2)
    M13_N2 = np.where(mask, f_M13_N2(N2_index), np.nan)
    eM13_N2 = np.where(mask, propagate_error(f_M13_N2, [N2_index], [eN2_index]), np.nan)
    
    ## -------------------------------- ##
    ## D16 calibrator (Dopita et al. 2016)
    ## -------------------------------- ##
    # To calculate the flux SII6716,30 = SII6716 + SII6730, it is necessary to revert to the linear flux values 10**(log(F)).
    # sum and error propagation
    Flux_SII6716_cor = 10**SII6716_cor
    Flux_SII6730_cor = 10**SII6730_cor
    F_total = Flux_SII6716_cor + Flux_SII6730_cor
    sigma1 = np.log(10) * Flux_SII6716_cor * eSII6716_cor
    sigma2 = np.log(10) * Flux_SII6730_cor * eSII6730_cor
    sigma_total = np.sqrt(sigma1**2 + sigma2**2)
    # Calculate the final flux and error of the sum of the SII6716,30 lines
    SII_total_cor = np.log10(F_total)
    eSII_total_cor = sigma_total / (F_total * np.log(10))
    # Here the actual calculation of the D16 calibrator and error propagation begins
    def f_D16(x, y, z):
        return 8.77 + (x - y) + 0.264 * (x - z)
    D16 = f_D16(NII6583_cor, SII_total_cor, Ha6562_cor)
    eD16 = propagate_error(f_D16, [NII6583_cor, SII_total_cor, Ha6562_cor], [eNII6583_cor, eSII_total_cor, eHa6562_cor])
    
    ## ------------------------------------------- ##
    ## T04 calibrator with R23 index (Tremonti et al. 2004)
    ## ------------------------------------------- ##
    # To calculate R23 index, it is necessary to revert to the linear flux values 10**(log(F)).
    Flux_O = 10**OII3727_cor + 10**OIII4958_cor + 10**OIII5006_cor
    sigma1 = np.log(10) * 10**OII3727_cor * eOII3727_cor
    sigma2 = np.log(10) * 10**OIII4958_cor * eOIII4958_cor
    sigma3 = np.log(10) * 10**OIII5006_cor * eOIII5006_cor
    eFlux_O_total = np.sqrt(sigma1**2 + sigma2**2 + sigma3**2)
    # Calculate the final flux and error of R23 index
    R23 = np.log10(Flux_O)
    eR23 = eFlux_O_total / (Flux_O * np.log(10))
    # Here the actual calculation of the T04 calibrator and error propagation begins
    def f_T04(x):
        return 9.185 - 0.313*x - 0.264*x**2 - 0.321*x**3
    T04 = f_T04(R23)
    eT04 = propagate_error(f_T04, [R23], [eR23])
    
    ## ---------------------------------------------- ##
    ##  KD02 calibrator with N2O2 index (R) (Kewley and Dopita 2002)
    ## ---------------------------------------------- ##
    def f_KD02(x):
        return np.log10(1.54020 + 1.26602*x + 0.167977*x**2) + 8.93
    mask = (N2O2 > -0.97) # O/H > 8.6
    KD02 = np.where(mask, f_KD02(N2O2), np.nan)
    eKD02 = np.where(mask, propagate_error(f_KD02, [N2O2], [eN2O2]), np.nan)
    
    ## ---------------------------------------------------------------------------- ##
    ##  P10_ONS and P10_ON calibrators with ONS and ON index (Pilyugin et al. 2010) ##
    ## ---------------------------------------------------------------------------- ##
    
    # P10_ONS calibrator
    def f_P10_ONS(P, R3, R2, N2, S2, eP, eR3, eR2, eN2, eS2):
    
        OH = np.full_like(N2, np.nan)
        eOH = np.full_like(N2, np.nan)
    
        mask1 = (N2 > -0.1)
        mask2 = (N2 < -0.1) & ((N2 - S2) > -0.25)
        mask3 = (N2 < -0.1) & ((N2 - S2) < -0.25)
    
        OH[mask1] = 8.277 + 0.657*P[mask1] - 0.399*R3[mask1] - 0.061*(N2[mask1] - R2[mask1]) + 0.005*(S2[mask1] - R2[mask1])
        eOH[mask1] = np.sqrt((0.657 * eP[mask1])**2 + (-0.399 * eR3[mask1])**2 + (-0.061 * eN2[mask1])**2 + (0.056 * eR2[mask1])**2 + (0.005 * eS2[mask1])**2)
        
        OH[mask2] = 8.816 - 0.733*P[mask2] + 0.454*R3[mask2] + 0.710*(N2[mask2] - R2[mask2]) - 0.337*(S2[mask2] - R2[mask2])
        eOH[mask2] = np.sqrt((-0.733 * eP[mask2])**2 + (0.454 * eR3[mask2])**2 + (0.710 * eN2[mask2])**2 + (-0.373 * eR2[mask2])**2 + (-0.337 * eS2[mask2])**2)
        
        OH[mask3] = 8.774 - 1.855*P[mask3] + 1.517*R3[mask3] + 0.304*(N2[mask3] - R2[mask3]) + 0.328*(S2[mask3] - R2[mask3])
        eOH[mask3] = np.sqrt((-1.855 * eP[mask3])**2 + (1.517 * eR3[mask3])**2 + (0.304 * eN2[mask3])**2 + (-0.632 * eR2[mask3])**2 + (0.328 * eS2[mask3])**2)
        
        return OH, eOH
    P10_ONS, eP10_ONS = f_P10_ONS(P, R3, R2, N2, S2, eP, eR3, eR2, eN2, eS2)

    # P10_ON calibrator
    def f_P10_ON(R3, R2, N2, S2, eR3, eR2, eN2):
    
        OH = np.full_like(N2, np.nan)
        eOH = np.full_like(N2, np.nan)
    
        mask1 = (N2 > -0.1)
        mask2 = (N2 < -0.1) & ((N2 - S2) > -0.25)
        mask3 = (N2 < -0.1) & ((N2 - S2) < -0.25)
    
        OH[mask1] = 8.606 - 0.105*R3[mask1] - 0.410*R2[mask1] - 0.150*(N2[mask1] - R2[mask1])
        eOH[mask1] = np.sqrt((-0.105 * eR3[mask1])**2 + (-0.260 * eR2[mask1])**2 + (-0.150 * eN2[mask1])**2)
        
        OH[mask2] = 8.642 + 0.077*R3[mask2] + 0.411*R2[mask2] + 0.601*(N2[mask2] - R2[mask2])
        eOH[mask2] = np.sqrt((0.077 * eR3[mask2])**2 + (-0.190 * eR2[mask2])**2 + (0.601 * eN2[mask2])**2)
        
        OH[mask3] = 8.013 + 0.905*R3[mask3] + 0.602*R2[mask3] + 0.751*(N2[mask3] - R2[mask3])
        eOH[mask3] = np.sqrt((0.905 * eR3[mask3])**2 + (-0.149 * eR2[mask3])**2 + (0.751 * eN2[mask3])**2)
        
        return OH, eOH
    P10_ON, eP10_ON = f_P10_ON(R3, R2, N2, S2, eR3, eR2, eN2)
    
    ## ----------------------------------------------------------------- ##
    ##  PM11 calibrator with NS index (Pilyugin and Mattsson 2011) ##
    ## ----------------------------------------------------------------- ##
    def f_PM11(R3, N2, S2, eR3, eN2, eS2):
        
        OH = np.full_like(N2, np.nan)
        eOH = np.full_like(N2, np.nan)
    
        mask1 = (N2 > -0.1)
        mask2 = (N2 < -0.1) & ((N2 - S2) > -0.25)
        mask3 = (N2 < -0.1) & ((N2 - S2) < -0.25)
        
        OH[mask1] = 8.454 - 0.216*R3[mask1] - 0.362*S2[mask1] - 0.101*(N2[mask1] - S2[mask1])
        eOH[mask1] = np.sqrt((-0.216 * eR3[mask1])**2 + (-0.101 * eN2[mask1])**2 + (-0.261 * eS2[mask1])**2)
        
        OH[mask2] = 8.456 + 0.082*R3[mask2] + 0.391*N2[mask2] + 0.290*(N2[mask2] - S2[mask2])
        eOH[mask2] = np.sqrt((0.082 * eR3[mask2])**2 + (0.681 * eN2[mask2])**2 + (-0.290 * eS2[mask2])**2)
        
        OH[mask3] = 7.881 + 0.929*R3[mask3] + 0.650*N2[mask3] + 0.025*(N2[mask3] - S2[mask3])
        eOH[mask3] = np.sqrt((0.929 * eR3[mask3])**2 + (0.675 * eN2[mask3])**2 + (-0.025 * eS2[mask3])**2)
        
        return OH, eOH
    PM11, ePM11 = f_PM11(R3, N2, S2, eR3, eN2, eS2)
    
    ## --------------------------------------------------------------- ##
    ##  PG16 calibrators with R and S index (Pilyugin and Grebel 2016) ##
    ## --------------------------------------------------------------- ##
    
    # PG16_R calibrator
    def f_PG16_R(R3, R2, N2, eR3, eR2, eN2):
    
        OH = np.full_like(N2, np.nan)
        eOH = np.full_like(N2, np.nan)
    
        mask1 = (N2 >= -0.6)
        mask2 = (N2 < -0.6)
        
        def f_PG16_R_upper(R3, R2, N2):
            return 8.589 + 0.022*(R3 - R2) + 0.399*N2 + (-0.137 + 0.164*(R3 - R2) + 0.589*N2) * R2
        OH[mask1] = f_PG16_R_upper(R3[mask1], R2[mask1], N2[mask1])
        eOH[mask1] = propagate_error(f_PG16_R_upper, [R3[mask1], R2[mask1], N2[mask1]], [eR3[mask1], eR2[mask1], eN2[mask1]])
        
        def f_PG16_R_lower(R3, R2, N2):
            return 7.932 + 0.944*(R3 - R2) + 0.695*N2 + (0.970 - 0.291*(R3 - R2) - 0.019*N2) * R2
        OH[mask2] = f_PG16_R_lower(R3[mask2], R2[mask2], N2[mask2])
        eOH[mask2] = propagate_error(f_PG16_R_lower, [R3[mask2], R2[mask2], N2[mask2]], [eR3[mask2], eR2[mask2], eN2[mask2]])
        
        return OH, eOH
    PG16_R, ePG16_R = f_PG16_R(R3, R2, N2, eR3, eR2, eN2)
    
    # PG16_S calibrator
    def f_PG16_S(R3, S2, N2, eR3, eS2, eN2):
    
        OH = np.full_like(N2, np.nan)
        eOH = np.full_like(N2, np.nan)
    
        mask1 = (N2 >= -0.6)
        mask2 = (N2 < -0.6)
        
        def f_PG16_S_upper(R3, S2, N2):
            return 8.424 + 0.030*(R3 - S2) + 0.751*N2 + (-0.349 + 0.182*(R3 - S2) + 0.508*N2) * S2
        OH[mask1] = f_PG16_S_upper(R3[mask1], S2[mask1], N2[mask1])
        eOH[mask1] = propagate_error(f_PG16_S_upper, [R3[mask1], S2[mask1], N2[mask1]], [eR3[mask1], eS2[mask1], eN2[mask1]])
        
        def f_PG16_S_lower(R3, S2, N2):
            return 8.072 + 0.789*(R3 - S2) + 0.726*N2 + (1.069 - 0.170*(R3 - S2) - 0.022*N2) * S2
        OH[mask2] = f_PG16_S_lower(R3[mask2], S2[mask2], N2[mask2])
        eOH[mask2] = propagate_error(f_PG16_S_lower, [R3[mask2], S2[mask2], N2[mask2]], [eR3[mask2], eS2[mask2], eN2[mask2]])
        
        return OH, eOH
    PG16_S, ePG16_S = f_PG16_S(R3, S2, N2, eR3, eS2, eN2)

    ## ----------------------------------------------------------------- ##
    
    ##  Curti et al. 2020 - usa ajuste do chi^2

    ##  Kobulnicky and Kewley 2004 - usa metalicidade
    
    ##  Ho et al. 2019 - usa machine learning
    
    ##  Thomas et al. 2018 - usa forma bayesiana
    
    ##  Pérez-Montero 2014 - usa temperatura eletrônica
    
    ## ----------------------------------------------------------------- ##

    #####################################################
    ################ Nitrogen Abundance #################
    #####################################################
    
    ## --------------------------------------------------------- ##
    ##  PG16 calibrators of NH and NO (Pilyugin and Grebel 2016) ##
    ## --------------------------------------------------------- ##
    
    # NH_PG16_R calibrator
    def f_NH_PG16_R(R3, R2, N2, eR3, eR2, eN2):
    
        NH = np.full_like(N2, np.nan)
        eNH = np.full_like(N2, np.nan)
    
        mask1 = (N2 >= -0.6)
        mask2 = (N2 < -0.6)
        
        def f_NH_PG16_R_upper(R3, R2, N2):
            return 7.939 + 0.135*(R3 - R2) + 1.217*N2 + (-0.765 + 0.166*(R3 - R2) + 0.449*N2) * R2
        NH[mask1] = f_NH_PG16_R_upper(R3[mask1], R2[mask1], N2[mask1])
        eNH[mask1] = propagate_error(f_NH_PG16_R_upper, [R3[mask1], R2[mask1], N2[mask1]], [eR3[mask1], eR2[mask1], eN2[mask1]])
        
        def f_NH_PG16_R_lower(R3, R2, N2):
            return 7.476 + 0.879*(R3 - R2) + 1.451*N2 + (-0.011 - 0.327*(R3 - R2) - 0.064*N2) * R2
        NH[mask2] = f_NH_PG16_R_lower(R3[mask2], R2[mask2], N2[mask2])
        eNH[mask2] = propagate_error(f_NH_PG16_R_lower, [R3[mask2], R2[mask2], N2[mask2]], [eR3[mask2], eR2[mask2], eN2[mask2]])
        
        return NH, eNH
    NH_PG16_R, eNH_PG16_R = f_NH_PG16_R(R3, R2, N2, eR3, eR2, eN2)
    
    # NO_PG16_R calibrator
    def f_NO_PG16_R(R2, N2, eR2, eN2):
    
        NO = np.full_like(N2, np.nan)
        eNO = np.full_like(N2, np.nan)
        
        def f_NO_PG16_R(R2, N2):
            return -0.657 - 0.201*N2 + (0.742 - 0.075*N2) * (N2 - R2)
        NO = f_NO_PG16_R(R2, N2)
        eNO = propagate_error(f_NO_PG16_R, [R2, N2], [eR2, eN2])
        
        return NO, eNO
    NO_PG16_R, eNO_PG16_R = f_NO_PG16_R(R2, N2, eR2, eN2)  
    
    ## ----------------------------------------------------------------- ##
    ##  NO_F22 calibrator with N2O2 index (Florido et al. 2022) ##
    ## ----------------------------------------------------------------- ##
    def f_NO_F22(a1, a2, a3, x):
        return a1 * x**2 + a2 * x - a3
    mask = (N2O2 > -1.74) & (N2O2 < 0.62)
    NO_F22 = np.where(mask, f_NO_F22(-0.102, 0.528, 0.634, N2O2), np.nan)
    eNO_F22 = np.where(mask, propagate_error(f_NO_F22, [-0.102, 0.528, 0.634, N2O2], [0.018, 0.019, 0.006, eN2O2]), np.nan)
    

    
    
    # --- Quality masks
    with np.errstate(invalid='ignore', divide='ignore'):
        rHb4861 = (eHb4861 < 0.997 * Hb4861) & (Hb4861 > 0)
        rHa6562 = (eHa6562 < 0.997 * Ha6562) & (Ha6562 > 0)
        rOII3727 = (eOII3727 < 0.997 * OII3727) & (OII3727 > 0)
        rOIII4958 = (eOIII4958 < 0.997 * OIII4958) & (OIII4958 > 0)
        rOIII5006 = (eOIII5006 < 0.997 * OIII5006) & (OIII5006 > 0)
        rNII6548 = (eNII6548 < 0.997 * NII6548) & (NII6548 > 0)
        rNII6583 = (eNII6583 < 0.997 * NII6583) & (NII6583 > 0)
        rSII6716 = (eSII6716 < 0.997 * SII6716) & (SII6716 > 0)
        rSII6730 = (eSII6730 < 0.997 * SII6730) & (SII6730 > 0)

    # --- Apply quality masks according to the index
    def apply_mask(oh, eoh, mask):
        oh_filtered = np.where(mask, oh, np.nan)
        eoh_filtered = np.where(mask, np.abs(eoh), np.nan)
        return oh_filtered, eoh_filtered

    if calibrator in [1, 4]:  # O3N2 index
        mask_ok = rHb4861 & rHa6562 & rOIII5006 & rNII6583
        if calibrator == 1:
            oh, eoh = apply_mask(PP04_O3N2, ePP04_O3N2, mask_ok)
        else:
            oh, eoh = apply_mask(M13_O3N2, eM13_O3N2, mask_ok)
            
    elif calibrator in [2, 3, 5]:  # N2 index
        mask_ok = rHa6562 & rNII6583
        if calibrator == 2:
            oh, eoh = apply_mask(PP04_N2, ePP04_N2, mask_ok)
        elif calibrator == 3:
            oh, eoh = apply_mask(PP04_N2_poly, ePP04_N2_poly, mask_ok)
        else:
            oh, eoh = apply_mask(M13_N2, eM13_N2, mask_ok)
            
    elif calibrator == 6:  # D16
        mask_ok = rHa6562 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(D16, eD16, mask_ok)
        
    elif calibrator == 7:  # T04
        mask_ok = rOII3727 & rOIII4958 & rOIII5006
        oh, eoh = apply_mask(T04, eT04, mask_ok)
        
    elif calibrator == 8:  # KD02
        mask_ok = rNII6583 & rOII3727
        oh, eoh = apply_mask(KD02, eKD02, mask_ok)
        
    elif calibrator == 9:  # P10_ONS
        mask_ok = rOII3727 & rOIII4958 & rOIII5006 & rNII6548 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(P10_ONS, eP10_ONS, mask_ok)
        
    elif calibrator == 10:  # P10_ON
        mask_ok = rOII3727 & rOIII4958 & rOIII5006 & rNII6548 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(P10_ON, eP10_ON, mask_ok)
        
    elif calibrator == 11:  # PM11
        mask_ok = rOIII4958 & rOIII5006 & rNII6548 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(PM11, ePM11, mask_ok)
        
    elif calibrator == 12:  # PG16_R
        mask_ok = rOII3727 & rOIII4958 & rOIII5006 & rNII6548 & rNII6583
        oh, eoh = apply_mask(PG16_R, ePG16_R, mask_ok)
        
    elif calibrator == 13:  # PG16_S
        mask_ok = rOIII4958 & rOIII5006 & rNII6548 & rNII6583 & rSII6716 & rSII6730
        oh, eoh = apply_mask(PG16_S, ePG16_S, mask_ok)
        
    elif calibrator == 14:  # NH_PG16_R
        mask_ok = rOII3727 & rOIII4958 & rOIII5006 & rNII6548 & rNII6583
        oh, eoh = apply_mask(NH_PG16_R, eNH_PG16_R, mask_ok)
        
    elif calibrator == 15:  # NO_PG16_R
        mask_ok = rOII3727 & rNII6548 & rNII6583
        oh, eoh = apply_mask(NO_PG16_R, eNO_PG16_R, mask_ok)
        
    elif calibrator == 16:  # NO_F22
        mask_ok = rOII3727 & rNII6583
        oh, eoh = apply_mask(NO_F22, eNO_F22, mask_ok)
        
    else:
        raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R, 16=NO_F22.")

    # --- Salve CSV
    results = pd.DataFrame({
        'HIIREGID': HIIREGID,
        'r': r,
        'EWHa6562': EWHa,
        'Ha6562_cor': Ha6562_cor,
        'eHa6562_cor': eHa6562_cor,
        'OII3727_cor': OII3727_cor,
        'eOII3727_cor': eOII3727_cor,
        'OIII4958_cor': OIII4958_cor,
        'eOIII4958_cor': eOIII4958_cor,
        'OIII5006_cor': OIII5006_cor,
        'eOIII5006_cor': eOIII5006_cor,
        'NII6548_cor': NII6548_cor,
        'eNII6548_cor': eNII6548_cor,
        'NII6583_cor': NII6583_cor,
        'eNII6583_cor': eNII6583_cor,
        'SII6716_cor': SII6716_cor,
        'eSII6716_cor': eSII6716_cor,
        'SII6730_cor': SII6730_cor,
        'eSII6730_cor': eSII6730_cor,
        'OH_PP04_O3N2': PP04_O3N2,
        'eOH_PP04_O3N2': ePP04_O3N2,
        'OH_PP04_N2': PP04_N2,
        'eOH_PP04_N2': ePP04_N2,
        'OH_PP04_N2_poly': PP04_N2_poly,
        'eOH_PP04_N2_poly': ePP04_N2_poly,
        'OH_M13_O3N2': M13_O3N2,
        'eOH_M13_O3N2': eM13_O3N2,
        'OH_M13_N2': M13_N2,
        'eOH_M13_N2': eM13_N2,
        'OH_D16': D16,
        'eOH_D16': eD16,
        'OH_T04': T04,
        'eOH_T04': eT04,
        'OH_KD02': KD02,
        'eOH_KD02': eKD02,
        'OH_P10_ONS': P10_ONS,
        'eOH_P10_ONS': eP10_ONS,
        'OH_P10_ON': P10_ON,
        'eOH_P10_ON': eP10_ON,
        'OH_PM11': PM11,
        'eOH_PM11': ePM11,
        'OH_PG16_R': PG16_R,
        'eOH_PG16_R': ePG16_R,
        'OH_PG16_S': PG16_S,
        'eOH_PG16_S': ePG16_S,
        'NH_PG16_R': NH_PG16_R,
        'eNH_PG16_R': eNH_PG16_R,
        'NO_PG16_R': NO_PG16_R,
        'eNO_PG16_R': eNO_PG16_R,
        'NO_F22': NO_F22,
        'eNO_F22': eNO_F22
    })
    
    
    if save_abundance:
    
        os.makedirs("calibrators", exist_ok=True)
        
        # --- Save the corrected flux in a table
        filepath = os.path.join("calibrators", f"{name}_calibrators.csv")

        results.to_csv(filepath, index=False)

    return oh, eoh, Ha6562_cor, OIII5006_cor, NII6583_cor
