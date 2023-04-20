#Visualization for FHN, with some work the visualizations could be merged with HH
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import myokit
import gc


def heatmap_lambda_rho_FHN_exp1(frequency):
    try: 
        data= pd.read_csv(f"DataFHN_HF/f{frequency}.csv")
    except:
        print("Missing simulation file" + f"DataFHN_HF/f{frequency}.csv")
        return

    data['rho'] = data['rho'].round(4)
    data['lambda'] = data['lambda'].round(4)
    
    data=data.drop(columns=['frequency']).drop_duplicates()
    plt.figure(figsize = (10,6))
    data=data.pivot("lambda", "rho", "num_peaks").astype(float).sort_values('lambda',ascending=False)
    graf = sns.heatmap(data,cmap="GnBu",annot=False,vmin=0, vmax=5)
    plt.title(rf'$\omega$={frequency} [kHz]')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\rho$')
    plt.savefig(f'Figure/FHN_HF/F{frequency}.png')
    plt.close()
    return 


def generate_polynomial_FHN_exp1(frequency, degree_fit_polynomial, range_lambda, range_rho):
    try: 
        data= pd.read_csv(f"DataFHN_HF/f{frequency}.csv")
    except:
        print("Missing simulation file" + f"DataFHN_HF/f{frequency}.csv")
        return
    #this is a patch to use old simulations, remove later
    #data.columns = ["index0","index","delta","rho","num_peaks","frequency","current"]
    
    data_new = data[data['num_peaks']==0]
    lambda_ = data_new['lambda'].unique()
    rho = []
    for l in lambda_:
        maximo = data_new[data_new['lambda']==l]['rho'].max()
        #if maximo != range_rho[1]-range_rho[2]:
        rho.append(maximo)
    if len(rho)==len(lambda_):
        coef_p = np.polyfit(lambda_,rho,degree_fit_polynomial)
    #else:
     #   coef_p = np.polyfit(lambda_[len(lambda_)-len(rho):],rho,degree_fit_polynomial)
    def polinomio(x):
        f = 0
        for i in list(np.arange(0,degree_fit_polynomial+1)):
            f=f+coef_p[i]*(x**(degree_fit_polynomial+1-(i+1)))
        return f
    plot_lambda = list(np.arange(range_lambda[0],range_lambda[1],range_lambda[2]))
    plt.plot(plot_lambda,[polinomio(i) for i in plot_lambda],
             label = rf'{frequency} $[kHz]$')
    return 
def approach_curve_FHN_exp1(set_frequency, degree_fit_polynomial, range_lambda, range_rho):
    for frequency in set_frequency:
        generate_polynomial_FHN_exp1(frequency, degree_fit_polynomial, range_lambda,range_rho)
    plt.xlim(range_lambda[0], range_lambda[1])
    plt.ylim(range_rho[0], range_rho[1])
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\rho$')
    plt.title('Approach Curve')
    plt.legend()
    plt.savefig(f'Figure/FHN_HF/Approach_curve_{frequency}kHz_{degree_fit_polynomial}pol.png')
    plt.close()
    return 