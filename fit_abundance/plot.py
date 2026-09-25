import matplotlib.pyplot as plt
import numpy as np
import os

def plot_model(
    results_dict,
    name,
    criterion,
    calibrator,
    save_graph,
    show_graph
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
        If True, saves the fit plots.
    show_graph : bool
        If True, displays the fit plot.

    Returns
    -------
    dict
        Dictionary containing the fitted gradient parameters.
    """
    
    x = results_dict['x']
    y = results_dict['y']
    ey = results_dict['ey']
    best_case = results_dict['best_case']

    fig, ax = plt.subplots()
    ax.errorbar(x, y, ey, fmt='o',
                color='black', 
                markerfacecolor='mediumslateblue',
                markeredgecolor='black',
                elinewidth=0.5,
                markersize=5,
                capsize=3,
                zorder=1)
    
    ymin = 0.97 * np.nanmin(y)
    ymax = 1.03 * np.nanmax(y)
    
    if best_case == 1:
        a2 = results_dict['fit1'][0]
        ea2 = results_dict['fit1'][1]
        b0 = results_dict['fit1'][2]
        eb0 = results_dict['fit1'][3]
        p = np.linspace(np.min(x), np.max(x), 100)
        ax.plot(p, a2*p + b0, color='red', linewidth=2, zorder=3)
        
        a1 = ea1 = h1 = eh1 = 0.0
        h2 = eh2 = a3 = ea3 = 0.0
    
    elif best_case == 2:
        fit2 = results_dict['fit2']
        fit2.plot_fit(color='red', linewidth=2, zorder=3)
        fit2.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)
        
        b0 = fit2.get_results()["estimates"]["const"]["estimate"]
        eb0 = fit2.get_results()["estimates"]["const"]["se"]
    
        slopes, h_break = [], []
        slope_errors, h_break_errors = [], []
        for i in range(1, best_case + 1):
            alpha = fit2.get_results()["estimates"][f"alpha{i}"]
            slopes.append(alpha["estimate"])
            slope_errors.append(alpha["se"])
        
        for i in range(1, best_case):
            breakpoint_ = fit2.get_results()["estimates"][f"breakpoint{i}"]
            h_break.append(breakpoint_["estimate"])
            h_break_errors.append(breakpoint_["se"])
            
        a1, a2, ea1, ea2 = slopes[0], slopes[1], slope_errors[0], slope_errors[1]
        h1, eh1 = h_break[0], h_break_errors[0]
        h2 = eh2 = a3 = ea3 = 0.0
        
    elif best_case == 3:
        fit3 = results_dict['fit3']
        fit3.plot_fit(color='red', linewidth=2, zorder=3)
        fit3.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)    
        
        b0 = fit3.get_results()["estimates"]["const"]["estimate"]
        eb0 = fit3.get_results()["estimates"]["const"]["se"]
        
        slopes, h_break = [], []
        slope_errors, h_break_errors = [], []
        for i in range(1, best_case + 1):
            alpha = fit3.get_results()["estimates"][f"alpha{i}"]
            slopes.append(alpha["estimate"])
            slope_errors.append(alpha["se"])
        
        for i in range(1, best_case):
            breakpoint_ = fit3.get_results()["estimates"][f"breakpoint{i}"]
            h_break.append(breakpoint_["estimate"])
            h_break_errors.append(breakpoint_["se"])
            
        a1, a2, a3, ea1, ea2, ea3 = slopes[0], slopes[1], slopes[2], slope_errors[0], slope_errors[1], slope_errors[2]
        h1, h2, eh1, eh2 = h_break[0], h_break[1], h_break_errors[0], h_break_errors[1]
        
    else:
        fit4 = results_dict['fit4']
        fit4.plot_fit(color='red', linewidth=2, zorder=3)
        fit4.plot_breakpoints(linewidth=0.8, color='gray', linestyle='-.', zorder=2)    
        
        b0 = fit4.get_results()["estimates"]["const"]["estimate"]
        eb0 = fit4.get_results()["estimates"]["const"]["se"]
        
        slopes, h_break = [], []
        slope_errors, h_break_errors = [], []
        for i in range(1, best_case + 1):
            alpha = fit4.get_results()["estimates"][f"alpha{i}"]
            slopes.append(alpha["estimate"])
            slope_errors.append(alpha["se"])
        
        for i in range(1, best_case):
            breakpoint_ = fit4.get_results()["estimates"][f"breakpoint{i}"]
            h_break.append(breakpoint_["estimate"])
            h_break_errors.append(breakpoint_["se"])
            
        a1, a2, a3, a4, ea1, ea2, ea3, ea4 = slopes[0], slopes[1], slopes[2], slopes[3], slope_errors[0], slope_errors[1], slope_errors[2], slope_errors[3]
        h1, h2, h3, eh1, eh2, eh3 = h_break[0], h_break[1], h_break[2], h_break_errors[0], h_break_errors[1], h_break_errors[2]
        
    ## Add Rectangle color of confidence_interval ##
    confidence_interval = []
    for i in range(1, best_case):
        
        if best_case == 2:
            bp = fit2.get_results()["estimates"][f"breakpoint{i}"]
        elif best_case == 3:
            bp = fit3.get_results()["estimates"][f"breakpoint{i}"]
        elif best_case == 4:
            bp = fit4.get_results()["estimates"][f"breakpoint{i}"]

        x0 = bp["confidence_interval"][0]
        x1 = bp["confidence_interval"][1]

        rect = plt.Rectangle(
            (x0, ymin),
            x1 - x0,
            ymax - ymin,
            facecolor="azure",
            zorder=0.5
        )

        ax.add_patch(rect)
        confidence_interval.append([x0, x1])
        
    ax.set_xlabel("$r$/r$_e$", size=14)
    if calibrator == 14:
        ax.set_ylabel("12+log(N/H)", size=14)
    elif calibrator in [15, 16]:
        ax.set_ylabel("log(N/O)", size=14)
    else:
        ax.set_ylabel("12+log(O/H)", size=14)
    ax.set_xlim(0, 1.05*np.max(x))
    if calibrator not in [15, 16]:
        ax.set_ylim([ymin, ymax])
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=4.0, width=0.7, colors='black', grid_color='gray', grid_alpha=0.9)
    ax.tick_params(which='minor', direction='in', length=2.0, width=0.5, colors='black', grid_color='gray', grid_alpha=0.9)
    
    if save_graph:
    
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
            raise ValueError("Invalid calibrator. Use 1=PP04_O3N2, 2=PP04_N2, 3=PP04_N2_poly, 4=M13_O3N2, 5=M13_N2, 6=D16, 7=T04, 8=KD02, 9=P10_ONS, 10=P10_ON, 11=PM11, 12=PG16_R, 13=PG16_S, 14=NH_PG16_R, 15=NO_PG16_R, 16=NO_F22.")
            
        os.makedirs("graphs", exist_ok=True)
        filepath = os.path.join("graphs", f"{name}_{calib}_{criterion}.png")
        
        ax.set_title(f'{name} - {calib}')
        plt.savefig(filepath, transparent=False, facecolor='w', edgecolor='w')
        
    # Show the plot or not
    if show_graph:
        plt.show(block=True)
    else:
        plt.close(fig)
       
    if best_case in [1, 2, 3]:
        return {
            'galaxy': name,
            'n_regions': len(x),
            'b0': b0,
            'eb0': eb0,
            'a1': a1,
            'ea1': ea1,
            'h1': h1,
            'eh1': eh1,
            'a2': a2,
            'ea2': ea2,
            'h2': h2,
            'eh2': eh2,
            'a3': a3,
            'ea3': ea3,
            'confidence_interval': confidence_interval
        }
    else:
        return {
            'galaxy': name,
            'n_regions': len(x),
            'b0': b0,
            'eb0': eb0,
            'a1': a1,
            'ea1': ea1,
            'h1': h1,
            'eh1': eh1,
            'a2': a2,
            'ea2': ea2,
            'h2': h2,
            'eh2': eh2,
            'a3': a3,
            'ea3': ea3,
            'h3': h3,
            'eh3': eh3,
            'a4': a4,
            'ea4': ea4,
            'confidence_interval': confidence_interval
        }
