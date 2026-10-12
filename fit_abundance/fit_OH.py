from .distance import distances
from .abundance import abundance
from .criteria import points
from .models import fit_models
from .plot import plot_model


def fit_final(
    name=None,
    *,
    HIIREGID=None,
    ra=None,
    ra0=None,
    dec=None,
    dec0=None,
    pa=None,
    ba=None,
    d=None,
    re=None,
    EWHa=None,
    Hb4861=None,
    eHb4861=None,
    Ha6562=None,
    eHa6562=None,
    OII3727=None,
    eOII3727=None,
    OIII4958=None,
    eOIII4958=None,
    OIII5006=None,
    eOIII5006=None,
    NII6548=None,
    eNII6548=None,
    NII6583=None,
    eNII6583=None,
    SII6716=None,
    eSII6716=None,
    SII6730=None,
    eSII6730=None,
    calibrator=None,
    criterion=None,
    analysis='full',
    n_boot=200,
    min_distance_to_edge=0.1,
    min_distance_between_breakpoints=0.15,
    start_values=None,
    plot_mode='best',
    save_abundance=False,
    save_oh_criteria=False,
    save_BPT=False,
    show_graph_bpt=False,
    save_model_selection=False,
    save_graph=False,
    show_graph=False,
):

    """
    Run the full abundance-gradient analysis pipeline.

        - The pipeline performs the following steps:

            1. Compute deprojected galactocentric distances
            2. Apply extinction correction to emission lines
               and compute oxygen and/or nitrogen abundances
            3. Select HII regions based on spectral criteria
            4. Fit abundance gradients
            5. Plot the best-fit model
    
        - The user has the option to choose only some of the following modules:

            'distance'
                Calculate only deprojected galactocentric distances.
            'abundance'
                Calculate only oxygen abundances.
            'distance+abundance'
                Calculate distances and oxygen and/or nitrogen abundances.
            'criteria'
                Calculate distances, abundances, and apply the selection criterion.
            'full'
                Run the complete pipeline, including gradient fitting and plotting.

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
    analysis : str
        Type of module the user wants to compute.
    n_boot : int
        Number of bootstrap iterations to be performed. By default, n_boot = 200.
    min_distance_to_edge : float
        Minimum distance from edge of data to a breakpoint, as a proportion of the data range
    min_distance_between_breakpoints : float
        The minimum required distance between breakpoints, as a proportion of the data range.
    start_values : List(float) or numpy(float)
        Initial guesses for breakpoint positions.
    plot_mode : str
        Type of statistical analysis the user wants to use for fitting the gradient.
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
    
    if calibrator in [14]:
        abund, eabund = 'NH', 'eNH'
    elif calibrator in [15, 16]:
        abund, eabund = 'NO', 'eNO'
    else:
        abund, eabund = 'OH', 'eOH'

    # ---------------------------------------------------
    # General validation
    # ---------------------------------------------------

    if not name:
        raise ValueError("Insert the name of the galaxy.")
        
    if analysis not in ['distance', 'abundance', 'distance+abundance', 'criteria', 'full']:
        raise ValueError("Invalid analysis. Choose 'distance', 'abundance', 'distance+abundance', 'criteria', or 'full'.")

    if plot_mode not in ['best', 'averaged']:
        raise ValueError("Invalid plot_mode. Choose either 'best' or 'averaged'.")

    # ===================================================
    # ANALYSIS = DISTANCE
    # ===================================================

    if analysis == 'distance':

        required_distance = {
            'ra': ra,
            'ra0': ra0,
            'dec': dec,
            'dec0': dec0,
            'pa': pa,
            'ba': ba,
            'd': d,
            're': re,
        }

        missing = [key for key, value in required_distance.items() if value is None]
        if missing:
            raise ValueError(f"Missing arguments for distance calculation: {missing}")

        r = distances(
            ra,
            ra0,
            dec,
            dec0,
            pa,
            ba,
            d,
            re,
        )

        return {
            'analysis': 'distance',
            'r': r,
        }

    # ===================================================
    # ANALYSIS = ABUNDANCE
    # ===================================================

    if analysis == 'abundance':

        required_abundance = {
            'Hb4861': Hb4861,
            'eHb4861': eHb4861,
            'Ha6562': Ha6562,
            'eHa6562': eHa6562,
            'OII3727': OII3727,
            'eOII3727': eOII3727,
            'OIII4958': OIII4958,
            'eOIII4958': eOIII4958,
            'OIII5006': OIII5006,
            'eOIII5006': eOIII5006,
            'NII6548': NII6548,
            'eNII6548': eNII6548,
            'NII6583': NII6583,
            'eNII6583': eNII6583,
            'SII6716': SII6716,
            'eSII6716': eSII6716,
            'SII6730': SII6730,
            'eSII6730': eSII6730,
        }

        missing = [key for key, value in required_abundance.items() if value is None]
        if missing:
            raise ValueError(f"Missing arguments for abundance calculation: {missing}")

        y, ey, Ha6562_cor, OIII5006_cor, NII6583_cor = abundance(
            name,
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
            save_abundance,
            calibrator=None,
            r=None,
            HIIREGID=HIIREGID,
            EWHa=EWHa,
        )

        return {
            'analysis': 'abundance - Make sure to choose save_abundance=True to save all the distances and abundances',
            #'galaxy': name,
            #str(abund): y,
            #str(eabund): ey,
            #'Ha6562_cor': Ha6562_cor,
            #'OIII5006_cor': OIII5006_cor,
            #'NII6583_cor': NII6583_cor,
            #'Make sure to choose save_abundance=True to save all the distances and abundances'
        }
        
    # ===================================================
    # ANALYSIS = DISTANCE + ABUNDANCE
    # ===================================================
        
    if analysis == 'distance+abundance':

        required_inputs = {
            'ra': ra,
            'ra0': ra0,
            'dec': dec,
            'dec0': dec0,
            'pa': pa,
            'ba': ba,
            'd': d,
            're': re,
            'Hb4861': Hb4861,
            'eHb4861': eHb4861,
            'Ha6562': Ha6562,
            'eHa6562': eHa6562,
            'OII3727': OII3727,
            'eOII3727': eOII3727,
            'OIII4958': OIII4958,
            'eOIII4958': eOIII4958,
            'OIII5006': OIII5006,
            'eOIII5006': eOIII5006,
            'NII6548': NII6548,
            'eNII6548': eNII6548,
            'NII6583': NII6583,
            'eNII6583': eNII6583,
            'SII6716': SII6716,
            'eSII6716': eSII6716,
            'SII6730': SII6730,
            'eSII6730': eSII6730,
        }

        missing = [key for key, value in required_inputs.items() if value is None]

        if missing:
            raise ValueError(f"Missing arguments for distance+abundance calculation: {missing}")

        r = distances(
            ra,
            ra0,
            dec,
            dec0,
            pa,
            ba,
            d,
            re,
        )

        oh, eoh, Ha6562_cor, OIII5006_cor, NII6583_cor = abundance(
            name,
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
            save_abundance,
            calibrator=None,
            HIIREGID=HIIREGID,
            EWHa=EWHa,
            r=r,
        )

        return {
            'analysis': 'distance+abundance - Make sure to choose save_abundance=True to save all the distances and abundances',
            #'galaxy': name,
            #str(abund): y,
            #str(eabund): ey,
            #'Ha6562_cor': Ha6562_cor,
            #'OIII5006_cor': OIII5006_cor,
            #'NII6583_cor': NII6583_cor,
        }

    # ===================================================
    # ANALYSIS = CRITERIA
    # ===================================================

    required_distance = {
        'ra': ra,
        'ra0': ra0,
        'dec': dec,
        'dec0': dec0,
        'pa': pa,
        'ba': ba,
        'd': d,
        're': re,
        'Hb4861': Hb4861,
        'eHb4861': eHb4861,
        'Ha6562': Ha6562,
        'eHa6562': eHa6562,
        'OII3727': OII3727,
        'eOII3727': eOII3727,
        'OIII4958': OIII4958,
        'eOIII4958': eOIII4958,
        'OIII5006': OIII5006,
        'eOIII5006': eOIII5006,
        'NII6548': NII6548,
        'eNII6548': eNII6548,
        'NII6583': NII6583,
        'eNII6583': eNII6583,
        'SII6716': SII6716,
        'eSII6716': eSII6716,
        'SII6730': SII6730,
        'eSII6730': eSII6730,
        'calibrator': calibrator,
        'EWHa': EWHa,
        'criterion': criterion,
    }

    missing = [key for key, value in required_distance.items() if value is None]
    if missing:
        raise ValueError(f"Missing arguments for calculation: {missing}")

    x = distances(
        ra,
        ra0,
        dec,
        dec0,
        pa,
        ba,
        d,
        re,
    )

    y, ey, Ha6562_cor, OIII5006_cor, NII6583_cor = abundance(
        name,
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
        save_abundance,
        calibrator=calibrator,
        HIIREGID=HIIREGID,
        EWHa=EWHa,
        r=x,
    )

    r, oh, eoh = points(
        name,
        criterion,
        x,
        y,
        ey,
        EWHa,
        Ha6562_cor,
        OIII5006_cor,
        NII6583_cor,
        calibrator,
        save_oh_criteria,
        save_BPT,
        show_graph_bpt,
    )

    if analysis == 'criteria':

        return {
            'analysis': 'criteria - Make sure to choose save_oh_criteria=True to save all the abundance of the specific calibrator.',
            #'r': r,
            #str(abund): oh,
            #str(eabund): eoh,
        }

    # ===================================================
    # FULL: FIT MODELS
    # ===================================================

    results_dict = fit_models(
        r,
        oh,
        eoh,
        name,
        criterion,
        calibrator,
        save_model_selection,
        n_boot=n_boot,
        min_distance_to_edge=min_distance_to_edge,
        min_distance_between_breakpoints=min_distance_between_breakpoints,
        start_values=start_values,
    )

    if results_dict is None:
        print(f"Insufficient data for fitting the galaxy {name}.")
        return None

    # ---------------------------------------------------
    # Plot
    # ---------------------------------------------------

    output = plot_model(
        results_dict,
        name,
        criterion,
        calibrator,
        save_graph,
        show_graph,
        plot_mode=plot_mode,
    )

    print(f"{name} completed.")

    return output
