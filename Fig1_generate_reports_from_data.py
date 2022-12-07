import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import os
#Visualization for HH

def heatmap_delta_rho_HH(current, frequency):
    """
    Parameters
    ----------
    
    Returns 
    -------

    """

    try: 
        data= pd.read_csv(f"DataHH/f{frequency}_I{current}.csv")
    except:
        print("Missing simulation file" + f"DataHH/f{frequency}_I{current}.csv")
        return
    
    data = pd.read_csv(f'DataHH/f{frequency}_I{current}.csv')
    data.columns = ["index0","index","delta","rho","num_peaks","frequency","current"]
    data['rho'] = data['rho'].round(2)
    data['delta'] = data['delta'].round(4)
    data=data.pivot("rho", "delta", "num_peaks").astype(float).sort_values('rho',ascending=False)
    plt.figure(figsize = (10,6))
    sns.heatmap(data,cmap="GnBu",annot=False,vmin=0, vmax=5)
    plt.title(rf'$I_0$={current} [$\mu$A/$cm^2$] $\omega$={frequency} [kHz]')
    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$\rho$')
    if not os.path.exists("Figure"):
        os.makedirs("Figure")
    if not os.path.exists("Figure/HH"):
        os.makedirs("Figure/HH")    
    plt.savefig(f'Figure/HH/I{current}_F{frequency}.png')
    #plt.show()
    plt.close()

    return 

def generate_polynomial_HH(current, frequency, degree_fit_polynomial, range_delta):
    """
    Parameters
    ----------
    
    Returns 
    -------

    """
    try: 
        data= pd.read_csv(f"DataHH/f{frequency}_I{current}.csv")
    except:
        print("Missing simulation file" + f"DataHH/f{frequency}_I{current}.csv")
        return
    
    #this is a patch to use old simulations, remove later
    data.columns = ["index0","index","delta","rho","num_peaks","frequency","current"]
    data_new = data[data['num_peaks']==0]

    rho = data_new['rho'].unique()
    delta = []

    for r in rho:
        maximo = data_new[data_new['rho']==r]['delta'].max()
        delta.append(maximo)
    coef_p = np.polyfit(delta, rho, degree_fit_polynomial)
    
    def polinomio(x):
        f = 0
        for i in list(np.arange(0,degree_fit_polynomial+1)):
            f=f+coef_p[i]*(x**(degree_fit_polynomial+1-(i+1)))
        return f
    plot_delta = list(np.arange(range_delta[0],range_delta[1],range_delta[2]))
    plt.plot(plot_delta,[polinomio(i) for i in plot_delta],
             label = rf'{frequency} $[kHz]$ - {current} $[\mu A/cm^2]$')
    #plt.show()
    return 

def approach_curve_HH(list_current, frequency, degree_fit_polynomial, range_delta, range_rho):
    """
    Parameters
    ----------
    
    Returns 
    -------

    """
    for current in list_current:
        generate_polynomial_HH(current, frequency, degree_fit_polynomial, range_delta)
    plt.xlim(range_delta[0], range_delta[1])
    plt.ylim(range_rho[0], range_rho[1])
    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$\rho$')
    plt.title('Approach Curve')
    plt.legend()
    if not os.path.exists("Figure"):
        os.makedirs("Figure")
    if not os.path.exists("Figure/HH"):
        os.makedirs("Figure/HH")    

    plt.savefig(f'Figure/HH/Approach_curve_{frequency}kHz_{degree_fit_polynomial}pol.png')
    plt.close()
    return


def main():
    # Execute experiment 1: amplitude of HF term vs slope for the HH model

    #parameter mesh
    range_delta = [0.0025,0.1275,0.0025]
    range_rho = [0.01,0.3,0.01]
    #range_delta = [0.0025,0.1275,0.01]
    #range_rho = [0.01,0.3,0.05]


    # Generate plots
    for current in [40, 80, 120]:
        for frequency in [3,5,10]:
            heatmap_delta_rho_HH(current, frequency) 
    for frequency in [3,5,10]:
        approach_curve_HH([40,80,120], frequency, degree_fit_polynomial = 6, range_delta= range_delta, range_rho = range_rho)
    
    return 
main()