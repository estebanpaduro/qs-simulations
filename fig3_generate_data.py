# This code has been downloaded from: 
# https://github.com/estebanpaduro/qs-simulations
# and it is used for the generation of the simulation data and figures
# in the article
# The impact of high frequency-based stability on the onset of action 
# potentials in neuron models - E. Cerpa, N. Corrales, M. Courdurier, 
# L. E. Medina, E. Paduro

import numpy as np
import pandas as pd
import myokit
import os
from scipy.signal import find_peaks
from scipy.optimize import fsolve

## Experiment 2: Simulation for FHN model
def simulation_FHN_exp2(model, protocol, func, time, threshold , parameters_hf, parameters_dc, initial):
    """Detect blocking or not for a set parameters.
    Parameters
    ----------
    model           - Required : base model of myokit in .mmt
    protocol        - Required : in case the external stimulus needs a protocol
                                {[level, start, lenght, period, multiplier], None}
    func            - Required : list function for external stimulus
                                [oscillatory function, constant function]
    time            - Required : time simulation
    threshold       - Required : threshold action potential
    parameters_hf   - Required : parameters oscillatory function
                                {{'parameter':value},{}}, if the parameter depends on other 
                                value so format of value is str
    parameters_dc   - Required : parameters constant function
                                {{'parameter':value},{}}, if the parameter depends on other 
                                value so format of value is str
    initial         - Required : initial conditions 
                                [v,w]  
    Returns 
    -------
    df_times: data simulation in myokit
    list_peaks: position peaks
    """
    
    mod = myokit.load_model(model)
    prot = myokit.Protocol()
    if protocol is not None:
        prot.schedule(protocol[0],protocol[1],protocol[2],protocol[3],protocol[4])

    ##### Block Stimulus #####  
    EXT = mod.get('EXT')
    ## oscillatory function ##
    #--------Parameters------#
    for i in parameters_hf:
        if True == EXT.can_add_variable(i):
            var = EXT.add_variable(i)
            var.set_rhs(parameters_hf[i])
        else:
            var = mod.get('EXT.'+i)
            var.set_rhs(parameters_hf[i])
        
    stim_hf = EXT.add_variable('stim_hf')
    stim_hf.set_rhs(func[0])
    ## constant function ##
    #------Parameters-----#
    for i in parameters_dc:
        if True == EXT.can_add_variable(i):
            var = EXT.add_variable(i)
            var.set_rhs(parameters_dc[i])
        else:
            var = mod.get('EXT.'+i)
            var.set_rhs(parameters_dc[i])
    
    stim_dc = EXT.add_variable('stim_dc')
    stim_dc.set_rhs(func[1])
    
    ## final stimulus ##
    stim = mod.get('EXT.stim') #units [uA/cm^2]
    stim.set_rhs('stim_hf+stim_dc')
    ###########################
    membrane = mod.get('membrane')
    V_prom = membrane.add_variable('V_prom')
    V_prom.set_rhs('V-EXT.rho*sin(EXT.omega*environment.t)')
    
    ## Initial Conditions ##
    v = mod.get('membrane.V')
    v.set_initial_value(initial[0])
    w = mod.get('membrane.W')
    w.set_initial_value(initial[1])

    v_AV = mod.get('membrane.V_A1')
    v_AV.set_initial_value(initial[0])
    w_AV = mod.get('membrane.W_A1')
    w_AV.set_initial_value(initial[1])
    
    ### Validate model ###
    #print(m.warnings())
    ######################
    
    ## Simulation ##  
    s = myokit.Simulation(mod,prot)
    
    ### step default ###
    s.set_max_step_size(dtmax=0.005)
    s.set_min_step_size(dtmin=None)
    s.set_tolerance( abs_tol=1e-06 , rel_tol=0.0001 ) 
    ########################
     
    ## Run simulation ##
    df_times = s.run(time)
    val = 100
    for vm in df_times['environment.t']:
        if vm >=100:
            val=vm
            break
    ## find peaks HF
    list_peaks, _ = find_peaks(df_times['membrane.V_A1'][df_times['environment.t'].index(val)+1:], height = threshold,
                                prominence = (1,None))
    return df_times, list_peaks


def parameters_simulation_FHN_exp2(current,delta,rho,frequency=10):
    """Parameters simulation with base model FHN.mmt
    Parameters
    ----------
    delta       - Required : slope HF function
    rho         - Required : amplitude high frequency
    frequency   - Optional : frequency HF component
    current     - Required : current DC component
    """
    model = 'FHN_exp2.mmt'
    protocol = None
    parameters_hf = {'rho': rho,
                           'pi':np.pi,
                           'omega' : frequency}
    oscillatory_function = 'rho*omega*cos(omega*t)'
    parameters_dc = {'delta': delta,
                           'I':current,
                           't0':'(1/delta)+100',
                           'f1':0,
                           'f2': 'delta*I*(t-100)',
                           'f3':'piecewise(t<100,f1,f2)'}
    constant_function = 'piecewise(t<t0,f3,I)'
    functions = [oscillatory_function,constant_function]
    time = 500 
    
    threshold = 1
    ## this condition ensures starting from a stable state ##
    def system_P1(variables):
        v1, w1 = variables
        gamma = 0.5
        beta = 0.8
        eq_1 = (1-(rho**2)/2)*v1-(v1**3)/3-w1
        eq_2 = v1-gamma*w1+beta
        return [eq_1,eq_2]
    cond_init = [0,0]
    v1, w1 = fsolve(system_P1,cond_init)
    initial_conditions = [v1,w1]
    ## simulation ##
    return {"model": model
           , "protocol": protocol
           , "parameters_hf": parameters_hf
           , "oscillatory_function": oscillatory_function
           , "parameters_dc": parameters_dc
           , "constant_function": constant_function
           , "functions": functions
           , "time": time
           , "threshold": threshold
           , "initial_conditions": initial_conditions
           }

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()            
 

def run_experiment2(current, range_delta, range_rho):
    delta_col = []
    rho_col = []
    num_peaks_col = []
    
    #tracking progress 
    delta_range = np.arange(range_delta[0],range_delta[1],range_delta[2])
    rho_range = np.arange(range_rho[0],range_rho[1],range_rho[2])
    total = len(delta_range)*len(rho_range)
    progress = 0
    print(f'Simulation for FHN with current: {current}')
    printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for delta_log in delta_range:
        delta = 10**delta_log
        for rho_log in rho_range:
            rho = 10**rho_log
            progress = progress+1
            param= parameters_simulation_FHN_exp2(current,delta,rho)
            _, list_peaks = simulation_FHN_exp2(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["parameters_dc"]
                                       , param["initial_conditions"])
            delta_col.append(delta_log)
            rho_col.append(rho_log)
            num_peaks_col.append(len(list_peaks))
            printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    data = pd.DataFrame({'delta_ln': delta_col,'rho_ln': rho_col,'num_peaks': num_peaks_col})
    data['current'] = current
    data.to_csv(f'Data/exp2_FHN_I_{current}_log10.csv')
    return data


def main():
    if not os.path.exists("Data"):
        os.makedirs("Data")
    range_delta = [-2.2,0.05,0.05]
    range_rho = [-1.1,0.22,0.02]
    set_currents = [0.5, 0.475, 0.45, 0.425, 0.4,
                    0.375, 0.35, 0.325, 0.3,
                    0.275, 0.25, 0.225, 0.2,
                    0.175, 0.15, 0.125, 0.1]
    for current in set_currents:
        run_experiment2(current, range_delta, range_rho)

if __name__ == "__main__":
    main()