import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import colorsys
import os
import matplotlib.patches as mpatches
import FHN_generate_data_exp1 as generate
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize, to_rgba
from matplotlib import cm

def generate_set1_colors(num_colors):
    cmap_set1 = plt.get_cmap('Set1')
    colors = [cmap_set1(i/num_colors) for i in range(num_colors)]
    return colors

def generate_cubehelix_colors(num_colors):
    cubehelix_palette = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
    colors = [to_rgba(c) for c in cubehelix_palette(np.linspace(0, 1, num_colors))]
    colors = ['#%02x%02x%02x' % (int(c[0] * 255), int(c[1] * 255), int(c[2] * 255)) for c in colors]
    return colors

def curve(set_beta,name_fig):
    """Boundary curve that separates the region with and without onset response
    Parameters
    ----------
    set_beta - Required : list with different beta values
    name_fig - Required : name under which he figure is saved
    Returns
    -------
    Save a figure in folder Figures
    """
    plt.figure(figsize = (12,10))
    colors = generate_cubehelix_colors(len(set_beta))
    i = 0
    for beta in set_beta:
        try:
            data= pd.read_csv(f"Data/exp1_FHN_beta_{beta}_log10.csv")
        except:
            print("Missing simulation file" + f"Data/exp1_FHN_beta_{beta}_log10.csv")
            return
        data['rho'] = 10**data['rho_ln']
        data['lambda'] = 10**data['lambda_ln']
        data_new = data[data['num_peaks']==1]
        lambda_ = data_new['lambda'].unique()
        rho = []

        for lam in lambda_:
            minimo = data_new[data_new['lambda']==lam]['rho'].min()
            rho.append(minimo)
        plt.plot(lambda_,rho,color=colors[i], linewidth=4)
        i = i + 1

    norm_current = Normalize(vmin=np.min(set_beta), vmax=np.max(set_beta))
    sm_current = cm.ScalarMappable(norm=norm_current, cmap=sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True))
    sm_current.set_array([])
    cbar = plt.colorbar(sm_current, ticks=np.linspace(np.min(set_beta), np.max(set_beta), 11))
    cbar.set_label(r'  $\beta$', fontsize=25,rotation=0)
    cbar.ax.tick_params(labelsize=20)
    cbar.set_ticks(np.linspace(min(set_beta), max(set_beta), 11))
    
    plt.xscale('log')
    plt.ylim(0,0.7)
    plt.xlim(0.01,1.5)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel(r'$\rho$',fontsize=25)
    plt.xlabel(r'$\lambda$',fontsize=25)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False,size=10)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, size=10)
    plt.savefig('Figures/' + name_fig, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    return

def region_lambda_rho_FHN_exp1(beta, list_rho_plot, list_lambda_plot, markers, name_fig):
    """Region of action potential generation for a given combination of rho and lambda
    Parameters
    ----------
    beta                - Required : value of beta
    list_rho_plot       - Required : list of rho to indicate particular cases
    list_lambda_plot    - Required : list of lambda to indicate particular cases
    markers             - Required : list of markers for particular cases
    name_fig            - Required : name under which he figure is saved
    Returns
    ------
    Save a figure in folder Figures
    """
    try:
        data= pd.read_csv(f"Data/exp1_FHN_beta_{beta}_log10.csv")
    except:
        print("Missing simulation file" + f"Data/exp1_FHN_beta_{beta}_log10.csv")
        return

    data['rho_ln'] = data['rho_ln'].round(4)
    data['lambda_ln'] = data['lambda_ln'].round(4)
    data=data.drop(columns=['frequency']).drop_duplicates()
    data['rho_ln'] = data['rho_ln'].replace(-0.0, 0)
    data['lambda_ln'] = data['lambda_ln'].replace(-0.0, 0)
    data['rho'] = 10**data['rho_ln']
    data['lambda'] = 10**data['lambda_ln']
    data_one =data[data['num_peaks']==1]
    lambda_ = data_one['lambda'].unique()
    rho_min_one = [data_one[data_one['lambda'] == lam]['rho'].min() for lam in lambda_]

    plt.figure(figsize = (10,10))
    plt.fill_between(lambda_,rho_min_one,data_one['rho'].max(), hatch='//', step='pre', edgecolor='black', facecolor='none', linewidth=0)
    
    not_legend = mpatches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', hatch=None)
    one_legend = mpatches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', hatch='//')
    colors = generate_set1_colors(len(list_rho_plot)*len(list_lambda_plot))
    j = 0
    for rho in list_rho_plot:
        for lambda_ in list_lambda_plot:
            plt.plot(lambda_, rho, color=colors[j], marker=markers[j],markersize=20)
            j = j + 1

    plt.ylim(0,0.7)
    plt.xlim(10**(-2),1.5)
    plt.xscale('log')
    plt.xlabel(r'$\lambda$',fontsize=30)
    plt.ylabel(r'$\rho$',fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, size=10)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, size=10)
    plt.legend(handles=[not_legend, one_legend], 
               labels=['Without Onset Activation', 'With Onset Activation'],
                       fontsize=25, bbox_to_anchor=(1, 0.05), 
                       loc='lower right', framealpha=1)
    if not os.path.exists("Figures"):
        os.makedirs("Figures") 
    plt.savefig('Figures/' + name_fig, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    return

def plot_simulations_FHN_exp1(beta, list_rho, list_lambda,markers,name_fig):
    """Plot simulation exp1 for list_rho and list_lambda
    Parameters
    ----------
    beta        - Required : value of beta
    list_rho    - Required : list of rho to indicate particular cases
    list_lambda - Required : list of lambda to indicate particular cases
    markers     - Required : list of markers for particular cases
    name_fig    - Required : name under which he figure is saved
    Returns
    -------
    Save a figure in folder Figures
    """

    colors = generate_set1_colors(len(list_rho)*len(list_lambda))
    legends = []
    j = 0
    for rho in list_rho:
        for lambda_ in list_lambda:
            param= generate.parameters_simulation_FHN_exp1(beta,rho,lambda_)
            df_times, _ = generate.simulation_FHN_exp1(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters"]
                                       , param["initial_conditions"]
                                       , param["beta"])

            x = df_times['environment.t']
            y = df_times['membrane.V_A1']
            plt.plot(x, y, linestyle='-', color=colors[j], linewidth=2)
            line = Line2D([0], [0], color=colors[j], linewidth=0, marker=markers[j], markersize=10)
            text = fr'$\rho={round(rho,2)}$' + '\n' + fr'$\lambda = {lambda_}$'
            legends.append((line, text))
            j = j + 1

    line = Line2D([0], [0], color='black', linestyle='--')
    text = 'Threshold'
    legends.append((line, text))
    plt.axhline(y=1, color='black', linestyle='--', label='Threshold')
    plt.ylim([-3,3])
    plt.xlabel('t',fontsize=15)
    plt.ylabel('V(t)',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend([legend[0] for legend in legends],
               [legend[1] for legend in legends],
               bbox_to_anchor=(1, 1), loc='upper left',framealpha=1,fontsize=15,frameon=False)
    if not os.path.exists("Figures"):
        os.makedirs("Figures") 
    plt.savefig('Figures/' + name_fig, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    return