# This code has been downloaded from: 
# https://github.com/estebanpaduro/qs-simulations
# and it is used for the generation of the simulation data and figures
# in the article
# The impact of high frequency-based stability on the onset of action 
# potentials in neuron models - E. Cerpa, N. Corrales, M. Courdurier, 
# L. E. Medina, E. Paduro

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import fig3_generate_data as generate
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, to_rgba
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


# generate palette
def generate_set1_colors(num_colors):
    cmap_set1 = plt.get_cmap('Set1')
    colors = [cmap_set1(i/num_colors) for i in range(num_colors)]
    return colors

def generate_cubehelix_colors(num_colors):
    cubehelix_palette = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
    colors = [to_rgba(c) for c in cubehelix_palette(np.linspace(0, 1, num_colors))]
    colors = ['#%02x%02x%02x' % (int(c[0] * 255), int(c[1] * 255), int(c[2] * 255)) for c in colors]
    return colors

def curve(set_current, name_fig):
    """Boundary curve that separates the region with and without onset response
    Parameters
    ----------
    set_current - Required : list of current, para ver resultados
    name_fig    - Required : name under which he figure is saved
    Returns
    -------
    Save a figure in folder Figures 
    """
    plt.figure(figsize = (12,10))
    colors = generate_cubehelix_colors(len(set_current))
    i = 0
    for current in set_current[::-1]:
        try: 
            data= pd.read_csv(f"Data/exp2_FHN_I_{current}_log10.csv")
        except:
            print("Missing simulation file " + f"Data/exp2_FHN_I_{current}_log10.csv")
            return
        data['rho'] = 10**data['rho_ln']
        data['delta'] = 10**data['delta_ln']
        data_off = data[data['num_peaks'] == 0]
        delta_off = sorted(list(data_off['delta'].unique()),reverse=True)
        rho_off = list(data_off['rho'].unique())
        delta = []
        rho = []

        for d in delta_off:
            minimo = data_off[data_off['delta'] == d]['rho'].min()
            rho.append(minimo)
            delta.append(d)
        df = pd.DataFrame({'delta':delta, 'rho':rho})
        if min(rho_off) <= data['rho'].min():
            df2 = df[df['rho'] == min(rho_off)].drop_duplicates(subset='rho')
            df1 = df[df['rho'] != min(rho_off)]
            df = pd.concat([df1, df2], ignore_index=True).reset_index(drop=True)
        
        plt.plot(df['delta'],df['rho'], color=colors[i], linewidth=3)
        i = i+1
    
    norm_current = Normalize(vmin=np.min(set_current), vmax=np.max(set_current))
    sm_current = cm.ScalarMappable(norm=norm_current, cmap=sns.cubehelix_palette(start=.5, rot=-.5, reverse= False, as_cmap=True))
    sm_current.set_array([])
    cbar = plt.colorbar(sm_current, ticks=np.linspace(np.min(set_current), np.max(set_current), 11))
    cbar.set_label(r'$I_0$', fontsize=25, rotation=0)
    cbar.ax.tick_params(labelsize=20)
    cbar.set_ticks(np.linspace(min(set_current), max(set_current), 11))
    
    plt.xlabel(r'$\delta$', fontsize=25)
    plt.ylabel(r'$\rho$', fontsize=25)
    plt.xscale('log')
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, size=10)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, size=10)
    


    plt.savefig('Figures/' + name_fig, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    return

def region_delta_rho_FHN_exp2(current, list_rho_plot, list_delta_plot, markers, name_fig):
    """Region of action potential generation for a given combination of rho and delta
    Parameters
    ----------
    current         - Required : value of current
    list_rho_plot   - Required : list of rho to indicate particular cases
    list_delta_plot - Required : list of delta to indicate particular cases
    markers         - Required : list of markers for particular cases
    name_fig        - Required : name under which he figure is saved
    Returns
    ------
    Save a figure in folder Figures
    """
    try: 
        data= pd.read_csv(f"Data/exp2_FHN_I_{current}_log10.csv")
    except:
        print("Missing simulation file" + f"Data/exp2_FHN_I_{current}_log10.csv")
        return
    data['rho_ln'] = data['rho_ln'].round(4)
    data['delta_ln'] = data['delta_ln'].round(4)
    data = data.drop(columns=['frequency','current']).drop_duplicates()
    data['rho_ln'] = data['rho_ln'].replace(-0.0, 0)
    data['delta_ln'] = data['delta_ln'].replace(-0.0, 0)
    data['rho'] = 10**data['rho_ln']
    data['delta'] = 10**data['delta_ln']
    data_one = data[data['num_peaks']==1]
    data_repetitive = data[data['num_peaks']>1]
    delta_one = data_one['delta'].unique()
    rho_max_one = [data_one[data_one['delta'] == lam]['rho'].max() for lam in delta_one]

    plt.figure(figsize = (10,10))
    plt.fill_between(delta_one,data_one['rho'].min(), rho_max_one, hatch='//', step='pre', edgecolor='black', facecolor='none', linewidth=0)
    plt.fill_between(data_repetitive['delta'], data_repetitive['rho'].max(), data_repetitive['rho'], hatch='o', step='pre', edgecolor='black', facecolor='none', linewidth=0)
    
    not_legend = mpatches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', hatch=None)
    one_legend = mpatches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', hatch='//')
    repetitive_legend = mpatches.Circle((0, 0), 0.1, facecolor='white', edgecolor='black', hatch='o')
    colors = generate_set1_colors(len(list_rho_plot)*len(list_delta_plot))
    j = 0
    for rho in list_rho_plot:
        for delta in list_delta_plot:
            plt.plot(delta, rho, color=colors[j], marker=markers[j], markersize=20)
            j = j + 1

    if len(list(data_repetitive['delta'].unique()))>0:
        plt.legend(handles=[not_legend, one_legend, repetitive_legend], 
               labels=['Without Initial Activation', 'With Initial Activation', 'Repetitive Action Potential'],
               fontsize = 25, bbox_to_anchor=(1, 0.05), loc ='lower right', framealpha=1)
    else:
        plt.legend(handles=[not_legend, one_legend], 
               labels=['Without Onset Activation', 'With Onset Activation'],
               fontsize=25, bbox_to_anchor=(1, 0.05), loc ='lower right', framealpha=1)
    
    plt.xscale('log')
    plt.xlabel(r'$\delta$', fontsize=25)
    plt.ylabel(r'$\rho$', fontsize=25)
    plt.xlim(0.005,)
    plt.ylim(bottom=0,top=1.4)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, labeltop=False, size=10)
    plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True, labelright=False, size=10)
    plt.savefig('Figures/' + name_fig, bbox_inches = 'tight', pad_inches = 0.5)
    plt.close()
    return

def plot_simulations_FHN_exp2(current,list_rho, list_delta,markers,name_fig):
    """Plot simulation exp2 for list_rho and list_delta
    Parameters
    ----------
    current     - Required : value of current
    list_rho    - Required : list of rho to indicate particular cases
    list_delta  - Required : list of delta to indicate particular cases
    markers     - Required : list of markers for particular cases
    name_fig    - Required : name under which he figure is saved
    Returns
    -------
    Save a figure in folder Figures
    """

    colors = generate_set1_colors(len(list_rho)*len(list_delta))
    legends = []
    j = 0
    for rho in list_rho:
        for delta in list_delta:
            param = generate.parameters_simulation_FHN_exp2(current,delta,rho)
            df_times, _ = generate.simulation_FHN_exp2(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["parameters_dc"]
                                       , param["initial_conditions"])

            x = df_times['environment.t']
            y = df_times['membrane.V_A1']
            line = Line2D([0], [0], color=colors[j], linewidth=0, marker=markers[j], markersize=10)
            text = fr'$\rho={round(rho,2)}$' + '\n' + fr'$\delta = {round(delta,2)}$'
            legends.append((line, text))
            plt.plot(x, y, linestyle='-', color=colors[j], linewidth=2)
            j = j+1

    line = Line2D([0], [0], color='black', linestyle='--')
    text = 'Threshold'
    legends.append((line, text))
    plt.axhline(y=1, color='black', linestyle='--', label='Threshold')
    plt.ylim([-3,3])
    plt.xlabel('t', fontsize=15)
    plt.ylabel('V(t)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend([legend[0] for legend in legends],[legend[1] for legend in legends],
               bbox_to_anchor=(1, 1), loc='upper left',framealpha=1, fontsize=15, frameon=False)  
    plt.savefig('Figures/' + name_fig, bbox_inches='tight', pad_inches=0.5)
    plt.close()
    return

def main():
    if not os.path.exists("Figures"):
        os.makedirs("Figures") 
    currents = [0.5,0.475,0.45,0.4,0.375,0.35,0.325,0.3,0.275,0.25,0.225,0.2,0.175,0.15,0.125,0.1]
    markers = ['*','o','P','s']

    # This block generate figure 3.a and 3.b using I_0 = 0.2
    list_rho_plot = [0.5,1.25]
    list_delta_plot = [0.02,0.3]
    region_delta_rho_FHN_exp2(0.2, list_rho_plot,list_delta_plot,markers,'Fig3a.png')
    plot_simulations_FHN_exp2(0.2, list_rho_plot, list_delta_plot,markers,'Fig3b.png')

    # This block generate figure 3.c and 3.d using I_0 = 0.4
    list_rho_plot = [0.5,1.25]
    list_delta_plot = [0.3,0.02]
    region_delta_rho_FHN_exp2(0.4, list_rho_plot,list_delta_plot,markers,'Fig3c.png')
    plot_simulations_FHN_exp2(0.4, list_rho_plot, list_delta_plot,markers,'Fig3d.png')

    # This block generate figure 3.e
    curve(currents,name_fig='Fig3e.png')

if __name__ == "__main__":
    main()