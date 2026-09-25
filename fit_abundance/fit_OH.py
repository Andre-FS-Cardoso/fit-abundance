from .distance import distances
from .abundance import abundance
from .criteria import points
from .models import fit_models
from .plot import plot_model

def fit_final(
    name,
    HIIREGID,
    ra,
    ra0,
    dec,
    dec0,
    pa,
    ba,
    d,
    re,
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
    criterion,
    *,
    n_boot=200,
    n_break=2,
    save_abundance = False,
    save_oh_criteria=False,
    save_BPT=False,
    show_graph_bpt=False,
    save_model_selection=False,
    save_graph=False,
    show_graph=False,
    ):
    
    """
    Run the full abundance-gradient analysis pipeline.

    The pipeline performs the following steps:

    1. Compute deprojected galactocentric distances
    2. Apply extinction correction to emission lines
       and compute oxygen abundances
    3. Select HII regions based on spectral criteria
    4. Fit abundance gradients
    5. Plot the best-fit model

    Parameters
    ----------
    name : str
        Name of the galaxy.
    HIIREGID : array-like
        Identifier of each HII region.
    ra : array-like
        Right ascension of the HII regions (degrees).
    ra0 : float
        Right ascension of the galaxy center (degrees).
    dec : array-like
        Declination of the HII regions (degrees).
    dec0 : float
        Declination of the galaxy center (degrees).
    pa : float
        Position angle of the galaxy (degrees).
    ba : float
        Minor-to-major axis ratio.
    d : float
        Distance to the galaxy (Mpc).
    re : float
        Effective radius of the galaxy (kpc).
    EWHa : array-like
        Equivalent width of Hα (Å).
    Hb4861, Ha6562, OII3727, OIII4958, OIII5006, NII6548, NII6583, SII6716, SII6730 : array-like
        Emission-line fluxes.
    eHb4861, eHa6562m eOII3727, eOIII4958, eOIII5006, eNII6548, eNII6583, eSII6716, eSII6730 : array-like
        Flux uncertainties.
    calibrator : int
        Abundance calibrator identifier.
    criterion : str
        Selection criterion for HII regions.
    n_boot : int
        Number of bootstrap iterations to be performed. By default, n_boot = 200.
    n_break : int
        Number of breaks to be fitted. By default, n_break = 2.
    save_abundance : bool
        If True, saves all corrected fluxes and abundances, along with their respective errors, to a CSV file.
    save_oh_criteria : bool
        If True, saves the abundances for a given criterion, along with their respective errors, to a CSV file.
    save_BPT : bool
        If True, save the BPT diagram plot.
    show_graph_bpt:
        If True, displays the BPT diagram plot.
    save_model_selection : bool
        If True, saves the table containing the statistical information of the model fits to a CSV file.
    save_graph : bool
        If True, saves the fit plots.
    show_graph : bool
        If True, displays the fit plot.

    Returns
    -------
    dict or None
        Dictionary containing the fitted gradient parameters,
        or None if the fit cannot be performed.
    """
    
    # ---------------------------------------------------
    # Step 1: distances
    # ---------------------------------------------------
    x = distances(ra, ra0, dec, dec0, pa, ba, d, re)
    
    # ---------------------------------------------------
    # Step 2: abundances
    # ---------------------------------------------------
    y, ey, Ha6562_cor, OIII5006_cor, NII6583_cor = abundance(
        name, x, HIIREGID, EWHa,
        Hb4861, eHb4861,
        Ha6562, eHa6562,
        OII3727, eOII3727,
        OIII4958, eOIII4958,
        OIII5006, eOIII5006,
        NII6548, eNII6548,
        NII6583, eNII6583,
        SII6716, eSII6716,
        SII6730, eSII6730,
        calibrator,
        save_abundance
        )
    
    # ---------------------------------------------------
    # Step 3: apply selection criteria
    # ---------------------------------------------------
    r, oh, eoh = points(
        name, criterion,
        x, y, ey,
        EWHa,
        Ha6562_cor,
        OIII5006_cor,
        NII6583_cor,
        calibrator,
        save_oh_criteria,
        save_BPT,
        show_graph_bpt
        )
    
    # ---------------------------------------------------
    # Step 4: fit gradient models
    # ---------------------------------------------------
    results_dict = fit_models(r, oh, eoh,
        name, criterion, calibrator,
        save_model_selection,
        n_boot=n_boot,
        n_break=n_break
        )
    
    if results_dict is None:
        print(f"Insufficient data for fitting the galaxy {name}.")
        return None
    
    # ---------------------------------------------------
    # Step 5: plot best model
    # ---------------------------------------------------
    output = plot_model(
        results_dict,
        name,
        criterion,
        calibrator,
        save_graph,
        show_graph
        )

    print(f"{name} completed.")

    return output
