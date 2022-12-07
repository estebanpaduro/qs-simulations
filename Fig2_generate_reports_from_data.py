import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import os

#Visualization for FHN, with some work the visualizations could be merged with HH
def heatmap_delta_rho_FHN(current, frequency):
    try: 
        data= pd.read_csv(f"DataFHN/f{frequency}_I{current}.csv")
    except:
        print("Missing simulation file" + f"DataFHN/f{frequency}_I{current}.csv")
        return
    #this is a patch to use old simulations, remove later
    data.columns = ["index0","index","delta","rho","num_peaks","frequency","current"]
    data['rho'] = data['rho'].round(2)
    data['delta'] = data['delta'].round(4)
    
    data=data.drop(columns=['frequency','current']).drop_duplicates()
    plt.figure(figsize = (10,6))
    data=data.pivot("rho", "delta", "num_peaks").astype(float).sort_values('rho',ascending=False)
    graf = sns.heatmap(data,cmap="GnBu",annot=False,vmin=0, vmax=5)
    plt.title(rf'$I_0$={current} $\omega$={frequency} [kHz]')
    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$\rho$')
    if not os.path.exists("Figure"):
        os.makedirs("Figure")
    if not os.path.exists("Figure/FHN"):
        os.makedirs("Figure/FHN")    
    plt.savefig(f'Figure/FHN/I{current}_F{frequency}.png')
    #plt.show()
    plt.close()
    
    return 




def generate_polynomial_FHN(current, frequency, degree_fit_polynomial, range_delta):

    try: 
        data= pd.read_csv(f"DataFHN/f{frequency}_I{current}.csv")
    except:
        print("Missing simulation file" + f"DataFHN/f{frequency}_I{current}.csv")
        return
    #this is a patch to use old simulations, remove later
    data.columns = ["index0","index","delta","rho","num_peaks","frequency","current"]
    data_new = data[data['num_peaks']==0]
    rho = data_new['rho'].unique()
    delta = []

    for r in rho:
        maximo = data_new[data_new['rho']==r]['delta'].max()
        delta.append(maximo)
    coef_p = np.polyfit(delta,rho,degree_fit_polynomial)
    
    def polinomio(x):
        f = 0
        for i in list(np.arange(0,degree_fit_polynomial+1)):
            f=f+coef_p[i]*(x**(degree_fit_polynomial+1-(i+1)))
        return f
    plot_delta = list(np.arange(range_delta[0],range_delta[1],range_delta[2]))
    return plt.plot(plot_delta,[polinomio(i) for i in plot_delta],
             label = rf'{frequency} $[kHz]$ - {current}')
     

def approach_curve_FHN(set_current, frequency, degree_fit_polynomial, range_delta, range_rho):
    for current in set_current:
        generate_polynomial_FHN(current, frequency, degree_fit_polynomial, range_delta)
    #There has to be a verification if all the files exists
    plt.xlim(range_delta[0], range_delta[1])
    plt.ylim(range_rho[0], range_rho[1])
    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$\rho$')
    plt.title('Approach Curve')
    plt.legend()
    plt.savefig(f'Figure/FHN/Approach_curve_{frequency}kHz_{degree_fit_polynomial}pol.png')
    plt.close()
    return 

def main():
    # Execute experiment 2: amplitude of HF term vs slope for the FHN model
  
    range_delta = [0.01,0.195,0.005]
    range_rho = [0.1,1.55,0.05]
    #range_delta = [0.01,0.2,0.05]
    #range_rho = [0.1,0.75,0.1]

    # Generate plots
    for frequency in [1,5,10,15,20]:
        for current in [0.15,0.3,0.5,0.5]:
            heatmap_delta_rho_FHN(current, frequency)
    
    #looks weird, maybe I broke something
    for frequency in [1,5,10,15,20]:   
        approach_curve_FHN([0.15,0.3,0.4,0.5], frequency, degree_fit_polynomial = 10, range_delta= range_delta, range_rho = range_rho)


main()