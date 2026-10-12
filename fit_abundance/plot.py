import matplotlib.pyplot as plt
import numpy as np
import os
import re

def calculate_abundance_error(calibrator_id, statistical_error):
    """
    Calculates the total error in chemical abundance determination using quadrature sum.
    
    Parameters:
    -----------
    calibrator_id : int
        The identification number of the calibrator (1 to 16).
    statistical_error : float or np.ndarray
        The observational error propagated from emission line fluxes.
        
    Returns:
    --------
    total_error : float or np.ndarray
        The final combined uncertainty (in dex).
    """
    
    # Quadratic sum: sqrt(sigma_stat^2 + sigma_intrin^2)
    # Assuming sigma_intrin = 0.02 dex to account for instrumental errors
    total_error = np.where(
        statistical_error > 0.02, 
        statistical_error, 
        np.sqrt(np.square(statistical_error) + np.square(0.02))
    )
    
    return total_error

def plot_model(
    results_dict,
    name,
    criterion,
    calibrator,
    save_graph,
    show_graph,
    *,
    plot_mode='best'
    ):

    """
    Generate a plot of the oxygen abundance gradient and the best-fit model.

    Parameters
    ----------
    results_dict : dict
        Dictionary containing the fitting results from `fit_models`.
    name : str
        Name of the galaxy.
    criterion : str
        Selection criterion applied to the HII regions.
    calibrator : int
        Abundance calibrator identifier.
    save_graph : bool
        If True, saves the fit plot.
    show_graph : bool
        If True, displays the fit plot.
    plot_mode : str
        Type of statistical analysis the user wants to use for fitting the gradient.

    Returns
    -------
    dict
        Dictionary containing the fitted gradient parameters.
    """
    
    # Selection of the calibrator
    if calibrator == 1:
        calib = 'PP04_O3N2'
    elif calibrator == 2:
        calib = 'PP04_N2'
    elif calibrator == 3:
        calib = 'PP04_N2_poly'
    elif calibrator == 4:
        calib = 'M13_O3N2'
    elif calibrator == 5:
        calib = 'M13_N2'
    elif calibrator == 6:
        calib = 'D16'
    elif calibrator == 7:
        calib = 'T04'
    elif calibrator == 8:
        calib = 'KD02'
    elif calibrator == 9:
        calib = 'P10_ONS'
    elif calibrator == 10:
        calib = 'P10_ON'
    elif calibrator == 11:
        calib = 'PM11'
    elif calibrator == 12:
        calib = 'PG16_R'
    elif calibrator == 13:
        calib = 'PG16_S'
    elif calibrator == 14:
        calib = 'NH_PG16_R'
    elif calibrator == 15:
        calib = 'NO_PG16_R'
    elif calibrator == 16:
        calib = 'NO_F22'
    else:
        raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R or 16=NO_F22.")
    
    x = results_dict['x']
    y = results_dict['y']
    ey_flux = results_dict['ey']
    #Add in quadrature the error from flux propagation and the intrinsic dispersion of the calibration
    ey = calculate_abundance_error(calibrator_id=calibrator, statistical_error=ey_flux)
    best_case = results_dict['best_case']
    #plot_mode = results_dict.get('plot_mode', 'best')
    y_fitted_combined = results_dict['y_fitted_combined']         
    predict_combined  = results_dict['predict_combined_func']
    w_nonlinear       = results_dict['w_nonlinear']
    profile           = results_dict['profile']
    description       = results_dict['description']
    AICs       = results_dict['AICs']
    weights       = results_dict['weights']
    #deltr = results_dict['deltr']

    fig, ax = plt.subplots(figsize=(8,6))
    plotline, caps, bars = ax.errorbar(x, y, ey, fmt='o',
                color='black', 
                markerfacecolor='mediumslateblue',
                markeredgecolor='black',
                elinewidth=0.9,
                markersize=8,
                capsize=3,
                zorder=1,
                alpha=0.7)
            
    for bar in bars:
        bar.set_alpha(0.3)
    for cap in caps:
        cap.set_alpha(0.3)
    
    ymin = 0.97 * np.nanmin(y)
    ymax = 1.03 * np.nanmax(y)
    
    # --- BACKUP variables ---
    #Case 1 (linear):
    a_1 = ea_1 = b_1 = eb_1 = np.nan
    #Case 2 (1-break):
    b_2 = eb_2 = a1_2 = ea1_2 = h1_2 = eh1_2 = a2_2 = ea2_2 = np.nan
    #Case 3 (2-break):
    b_3 = eb_3 = a1_3 = ea1_3 = h1_3 = eh1_3 = a2_3 = ea2_3 = h2_3 = eh2_3 = a3_3 = ea3_3 = np.nan

    # Fix galaxy name to plot
    name_to_plot = str(name)

    match = re.match(r"^(NGC|IC|UGC)\s*(\d+.*)$", name_to_plot, re.IGNORECASE)

    if match:
        # Reconstruct the name with exactly one space between the prefix and the numbers
        prefix = match.group(1).upper() # Ensures standard uppercase (NGC/IC/UGC)
        number = match.group(2)
        name_to_plot = f"{prefix} {number}"
        
    # ---------------------------------------------------
    #Get variables
    # ---------------------------------------------------
    
    #Case 1:
    a_1 = results_dict['fit1'][0]; ea_1 = results_dict['fit1'][1]
    b_1 = results_dict['fit1'][2]; eb_1 = results_dict['fit1'][3]
    
    #Case 2:
    if results_dict['fit2'].get_results()["converged"]:
        est2 = results_dict['fit2'].get_results()["estimates"]
        b_2 = est2["const"]["estimate"]; eb_2 = est2["const"]["se"]
        a1_2 = est2["alpha1"]["estimate"]; ea1_2 = est2["alpha1"]["se"]
        h1_2 = est2["breakpoint1"]["estimate"]; eh1_2 = est2["breakpoint1"]["se"]
        a2_2 = est2["alpha2"]["estimate"]; ea2_2 = est2["alpha2"]["se"]            
    
    #Case 3:    
    if results_dict['fit3'].get_results()["converged"]:
        est3 = results_dict['fit3'].get_results()["estimates"]
        b_3 = est3["const"]["estimate"]; eb_3 = est3["const"]["se"]
        a1_3 = est3["alpha1"]["estimate"]; ea1_3 = est3["alpha1"]["se"]
        h1_3 = est3["breakpoint1"]["estimate"]; eh1_3 = est3["breakpoint1"]["se"]
        a2_3 = est3["alpha2"]["estimate"]; ea2_3 = est3["alpha2"]["se"]
        h2_3 = est3["breakpoint2"]["estimate"]; eh2_3 = est3["breakpoint2"]["se"]
        a3_3 = est3["alpha3"]["estimate"]; ea3_3 = est3["alpha3"]["se"]
    
    # ---------------------------------------------------
    # A: Model Averaging
    # ---------------------------------------------------
    if plot_mode == 'averaged':
        # Generate the points to plot the averaged prediction
        p_smooth = np.linspace(np.min(x), np.max(x), 200)
        y_smooth_combinado = results_dict['predict_combined_func'](p_smooth)
        
        # Plot of model Akaike weighted average
        ax.plot(p_smooth, y_smooth_combinado, color='darkorange', linewidth=2.5, linestyle=':', zorder=3, label='Akaike Averaged')
        ax.text(0.8,0.92,str(profile),fontsize=18,transform=ax.transAxes)
        
        ax.axvline(h1_2, color='gray', linestyle='--', alpha=0.7, zorder=2)
        ax.axvline(h1_3, color='gray', linestyle=':', alpha=0.7, zorder=2)
        ax.axvline(h2_3, color='gray', linestyle=':', alpha=0.7, zorder=2)
            
        print(f"[Model Averaged Mode] Best single structural candidate was Case {best_case}")
        
    # ---------------------------------------------------
    # B: Best model
    # ---------------------------------------------------
    else:
        if best_case == 1:
            p = np.linspace(np.min(x), np.max(x), 100)
            ax.plot(p, a_1*p + b_1, color='red', linewidth=2, zorder=3)
        
        elif best_case == 2:
            fit2 = results_dict['fit2']
            fit2.plot_fit(color='red', linewidth=2, zorder=3)
            fit2.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)
            
            ax.axvline(h1_2, color='gray', linestyle='--', alpha=0.7, zorder=2)
            
            interval0 = fit2.get_results()["estimates"]["breakpoint1"]["confidence_interval"][0]
            interval1 = fit2.get_results()["estimates"]["breakpoint1"]["confidence_interval"][1]
            rect1 = plt.Rectangle((interval0, ymin), (interval1-interval0), ymax-ymin, facecolor="azure", zorder=0.5)
            ax.add_patch(rect1)
            
        else:
            fit3 = results_dict['fit3']
            fit3.plot_fit(color='red', linewidth=2, zorder=3)
            fit3.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)
            
            ax.axvline(h1_3, color='gray', linestyle=':', alpha=0.7, zorder=2)
            ax.axvline(h2_3, color='gray', linestyle=':', alpha=0.7, zorder=2)
            
            interval0 = fit3.get_results()["estimates"]["breakpoint1"]["confidence_interval"][0]
            interval1 = fit3.get_results()["estimates"]["breakpoint1"]["confidence_interval"][1]
            rect1 = plt.Rectangle((interval0, ymin), (interval1-interval0), ymax-ymin, facecolor="azure", zorder=0.5)
            ax.add_patch(rect1)
            
            interval00 = fit3.get_results()["estimates"]["breakpoint2"]["confidence_interval"][0]
            interval11 = fit3.get_results()["estimates"]["breakpoint2"]["confidence_interval"][1]
            rect2 = plt.Rectangle((interval00, ymin), (interval11-interval00), ymax-ymin, facecolor="azure", zorder=0.5)
            ax.add_patch(rect2)
            
    ax.text(0.05,0.92,str(name_to_plot)+' - '+str(calib)+' - '+str(criterion),fontsize=18,transform=ax.transAxes)
    ax.set_xlabel("$r$/r$_e$", size=18)
    if calibrator in [14]:
        ax.set_ylabel("12+log(N/H)", size=18)
    elif calibrator in [15, 16]:
        ax.set_ylabel("log(N/O)", size=18)
    else:
        ax.set_ylabel("12+log(O/H)", size=18)
    ax.set_xlim(0, 1.05*np.max(x))
    ax.set_ylim([ymin, ymax])
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=4.0, width=0.7, colors='black', grid_color='gray', grid_alpha=0.9)
    ax.tick_params(which='minor', direction='in', length=2.0, width=0.5, colors='black', grid_color='gray', grid_alpha=0.9)
    ax.tick_params(axis='both', labelsize=16)
    
    ## Choose to save plot        
    if save_graph:
    
        save_folder = 'graphs'
        os.makedirs(save_folder, exist_ok=True)
        
        # Save the weighted average plot with different name
        suffix = "_averaged" if plot_mode == 'averaged' else ""
        filepath = os.path.join(save_folder, f"{name}_{calib}_{criterion}{suffix}.png")
        plt.savefig(filepath, transparent=False, facecolor='w', edgecolor='w', dpi=300)
        
    # Show the plot or not
    if show_graph:
        #ax.set_title(f'{name} - {calib}')
        plt.show(block=True)
    else:
        plt.close(fig)

    return {
        'galaxy': name,
        'b_1': b_1,
        'eb_1': eb_1,
        'a_1': a_1,
        'ea_1': ea_1,
        'b_2': b_2,
        'eb_2': eb_2,
        'a1_2': a1_2,
        'ea_2': ea1_2,
        'h1_2': h1_2,
        'eh1_2': eh1_2,
        'a2_2': a2_2,
        'ea2_2': ea2_2,
        'b_3': b_3,
        'eb_3': eb_3,
        'a1_3': a1_3,
        'ea_3': ea1_3,
        'h1_3': h1_3,
        'eh1_3': eh1_3,
        'a2_3': a2_3,
        'ea2_3': ea2_3,
        'h2_3': h2_3,
        'eh2_3': eh2_3,
        'a3_3': a3_3,
        'ea3_3': ea3_3,
        'best_case': best_case,
        'AIC_1': np.nan_to_num(AICs[0], posinf=np.nan, neginf=np.nan),
        'AIC_2': np.nan_to_num(AICs[1], posinf=np.nan, neginf=np.nan),
        'AIC_3': np.nan_to_num(AICs[2], posinf=np.nan, neginf=np.nan),
        'weights_1':weights[0],
        'weights_2':weights[1],
        'weights_3':weights[2],
        'w_nonlinear':  w_nonlinear ,
        'profile': profile
    }
