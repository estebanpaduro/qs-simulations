import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import myokit
import gc
import FHN_exp1_generate_data_heatmap as FHN_exp1_generate
import FHN_exp2_generate_data_heatmap as FHN_exp2_generate
import HH_exp3_generate_data_heatmap as HH_exp3_generate
import HH_exp4_generate_data_heatmap as HH_exp4_generate
import FHN_exp1_report_from_data as FHN_exp1_report
import FHN_exp2_report_from_data as FHN_exp2_report
import HH_exp3_report_from_data as HH_exp3_report
import HH_exp4_report_from_data as HH_exp4_report



def main(experiment, list_frequency, list_current, range_rho, range_lambda=[],range_delta=[]):
    """
    Main function for run experiments
    Parameters
    ----------
    experiment      - Required : number of experiment (Int)
    list_frequency  - Required : list frequencies (Int) 
    list_current    - Required : list currents (Int)
    range_rho       - Required : Use for both experiments (Float)
                                list [star, stop, step]
    range_lambda    - Optional : Only use for experiment 1 (Float)
                                list [star, stop, step]        
    range_delta     - Optional : Only use for experiment 2 (Float)
                                list [star,stop,step]

    Returns
    ------
    Generate csv files for each current and frequency
    and generate figures for each current and frequency
    """
    if experiment==1:
        #Execute experiment 1: FHN model with I(t)=lambda*rho*omega*t*cos(omega*t)
        for frequency in list_frequency:
            # Generate simulation  
            output = FHN_exp1_generate.run_experiment1(frequency,range_rho,range_lambda) 
            # Generate plots
            FHN_exp1_report.heatmap_lambda_rho_FHN_exp1(frequency)
        ############## To fix the next function ############
        FHN_exp1_report.approach_curve_FHN_exp1(list_frequency, degree_fit_polynomial = 10, range_lambda= range_lambda,range_rho = range_rho)
    
    elif experiment==2:
        # Execute experiment 2: FHN model with I(t)=rho*omega*cos(omega*t)+delta*I_0 
        # Generate simulation  
        for frequency in list_frequency:
            for current in list_current:
                output = FHN_exp2_generate.run_experiment2(current,frequency,range_delta,range_rho) 
        # Generate plots
                FHN_exp2_report.heatmap_delta_rho_FHN_exp2(current,frequency)
            FHN_exp2_report.approach_curve_FHN_exp2(list_current, frequency, degree_fit_polynomial = 4, range_delta= range_delta, range_rho = range_rho)
    
    elif experiment==3:
        #Execute experiment 3: HH model with I(t)=lambda*rho*omega*t*cos(omega*t)
        for frequency in list_frequency:
            # Generate simulation  
            output = HH_exp3_generate.run_experiment3(frequency,range_rho,range_lambda) 
            # Generate plots
            HH_exp3_report.heatmap_lambda_rho_HH_exp3(frequency)
        ############## To fix the next function ############
        HH_exp3_report.approach_curve_HH_exp3(list_frequency, degree_fit_polynomial = 10, range_lambda= range_lambda,range_rho = range_rho)
    
    elif experiment==4:
        # Execute experiment 4: HH model with I(t)=rho*omega*cos(omega*t)+delta*I_0 
        # Generate simulation 
        for frequency in list_frequency:
            for current in list_current:
                output = HH_exp4_generate.run_experiment4(current,frequency,range_delta,range_rho) 
        # Generate plots
                HH_exp4_report.heatmap_delta_rho_HH_exp4(current,frequency)
            HH_exp4_report.approach_curve_HH_exp4(list_current, frequency, degree_fit_polynomial = 4, range_delta= range_delta, range_rho = range_rho)
    
    else:
        print(f"Error: Experiment number {exp} doesn't exists")

    return

#   If you want plot Action Potential
def plot_single_simulation(experiment, current=0,frequency=0,rho=0,delta=0,lambda_=0):
    if experiment == 1:
        FHN_exp1_generate.plot_single_simulation_FHN_exp1(frequency, rho, lambda_)
    if experiment == 2:
        FHN_exp2_generate.plot_single_simulation_FHN_exp2(current, frequency, delta , rho)
    if experiment == 3:
        HH_exp3_generate.plot_single_simulation_HH_exp3(frequency, rho, lambda_)
    if experiment == 4:
        HH_exp4_generate.plot_single_simulation_HH_exp4(current, frequency, delta , rho)
    else:
        print(f"Error: Experiment number {exp} doesn't exists")
    return


"""
Change the parameters
"""

exp=2

range_rho = [0.025,0.825,0.025]
#range_lambda = [0.025,1,0.025]
range_delta = [0.01,0.2,0.005] 

list_frequency = [5]
list_current = [0.3]

main(exp,list_frequency,list_current,range_rho =range_rho,range_delta=range_delta)


"""
plot single example
exp = 1
current = 0.3
frequency = 5
#delta = 0.1
lambda_ = 0.1
rho = 0.2
plot_single_simulation(experiment, current=current,frequency=frequency,rho=rho,lambda_=lambda_)
"""
