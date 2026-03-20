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
    OIII5006,
    eOIII5006,
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
    save_abundance = False,
    save_oh_criteria=False,
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
    ra, dec : array-like
        Coordinates of the HII regions (degrees).
    ra0, dec0 : float
        Coordinates of the galaxy center (degrees).
    pa : float
        Position angle of the galaxy (degrees).
    ba : float
        Minor-to-major axis ratio.
    d : float
        Distance to the galaxy (Mpc).
    re : float
        Effective radius of the galaxy (kpc).
    calibrator : int
        Abundance calibrator identifier.
    criterion : str
        Selection criterion for HII regions.
    save_table : bool
        Save filtered HII region table.
    save_graph : bool
        Save gradient plot.
    show_graph : bool
        Display the plot.

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
        OIII5006, eOIII5006,
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
        save_oh_criteria
        )
    
    # ---------------------------------------------------
    # Step 4: fit gradient models
    # ---------------------------------------------------
    results_dict = fit_models(r, oh, eoh,
        name, criterion, calibrator,
        save_model_selection,
        n_boot=n_boot)
    
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
