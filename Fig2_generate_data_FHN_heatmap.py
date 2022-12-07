import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import myokit
import os

## Experiment 2: Simulation for FHN model
def simulation_FHN(model, protocol, func, time, threshold , parameters_hf, parameters_dc, initial):
    """Detect blocking or not for a set parameters.
    Parameters
    ----------
    mod : base model of myokit in .mmt
    prot : in case the external stimulus needs a protocol
        {[level, start, lenght, period, multiplier], None}
    func : list function for external stimulus
        [oscillatory function, constant function]
    time : time simulation
    peaks : find peaks (action potential)
            [height, distance]
    max_peaks : integer
    par_o : parameters oscillatory function
        {{'parameter':value},{}}, if the parameter depends on other 
        value so format of value is str
    par_c : parameters constant function
        {{'parameter':value},{}}, if the parameter depends on other 
        value so format of value is str
    initial: initial conditions
            [v,w]
    value : bool
        plot data in matplotlib figure.
    
    Returns 
    -------
    {0,1} : 0 if there is blocking and 1 if there is not blocking
    d: data simulation in myokit
    num_peaks: peaks with only HF or HF+AMP
    """
    
    mod = myokit.load_model(model)
    prot = myokit.Protocol()
    if protocol is not None:
        prot.schedule(protocol[0],protocol[1],protocol[2],protocol[3],protocol[4])

    ##### Block Stimulus #####  
    ST = mod.get('ST')
    ## oscillatory function ##
    #--------Parameters------#
    for i in parameters_hf:
        var = ST.add_variable(i)
        var.set_rhs(parameters_hf[i])
        
    stim_hf = ST.add_variable('stim_hf')
    stim_hf.set_rhs(func[0])
    ## constant function ##
    #------Parameters-----#
    for i in parameters_dc:
        var = ST.add_variable(i)
        var.set_rhs(parameters_dc[i])
    
    stim_dc = ST.add_variable('stim_dc')
    stim_dc.set_rhs(func[1])
    
    ## final stimulus ##
    stim = mod.get('ST.stim') #units [uA/cm^2]
    stim.set_rhs('stim_hf+stim_dc')
    ###########################
    membrane = mod.get('membrane')
    V_prom = membrane.add_variable('V_prom')
    V_prom.set_rhs('V-ST.rho*sin(ST.omega*environment.t)')
    
    ## Initial Conditions ##
    v = mod.get('membrane.V')
    v.set_state_value(initial[0])
    w = mod.get('membrane.W')
    w.set_state_value(initial[1])
    
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
    dist=(1/(50/s.last_number_of_steps()))
    val = 50
    for vm in df_times['environment.t']:
        if vm >=50:
            val=vm
            break
    ## find peaks HF
    list_peaks, _ = find_peaks(df_times['membrane.V_prom'][df_times['environment.t'].index(val)+1:], height = threshold,
                               distance = dist)
  
    return df_times, list_peaks


def parameters_simulation_FHN(current,frequency,delta,rho):
    """parameters simulation with base model HH.mmt
    Parameters
    ----------
    delta : slope constant function
    rho : amplitude high frequency
    current : current with slope delta
    frequency: frequency oscillatory function
    value : bool
        plot data in matplotlib figure.
    
    Returns 
    -------
    {0,1} : 0 if there is blocking and 1 if there is not blocking
    d: data simulation in myokit
    num_peaks: peaks with only HF, only AMP or HF+AMP
    """
    model = 'FHN.mmt'
    protocol = None
    parameters_hf = {'rho': rho,
                           'pi':np.pi,
                           'omega':frequency}
    oscillatory_function = 'rho*omega*cos(omega*t)'
    parameters_dc = {'m': delta,
                           'C':current,
                           't0':'(1/m)+50',
                           'f1':0,
                           'f2': '(t-50)*m*C',
                           'f3':'piecewise(t<50,f1,f2)'}
    constant_function = 'piecewise(t<t0,f3,C)'
    functions = [oscillatory_function,constant_function]
    time = 200 # (current/delta)+150
    
    threshold = 1 
    ## this condition ensures starting from a stable state ##
    initial_conditions = [-1.12,-0.65]
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



    #return bloqueo, data, num_peaks, peaks

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
 

def run_experiment2(current, frequency, range_delta, range_rho):
    delta_col = []
    rho_col = []
    num_peaks_col = []
    
    #tracking progress 
    delta_range = np.arange(range_delta[0],range_delta[1],range_delta[2])
    rho_range = np.arange(range_rho[0],range_rho[1],range_rho[2])
    total = len(delta_range)*len(rho_range)
    progress = 0
    print(f'Simulation for FHN with current: {current}, and frequency: {frequency}')
    printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for delta in delta_range:
        for rho in rho_range:
            progress = progress+1
            param= parameters_simulation_FHN(current,frequency,delta,rho)
            _, list_peaks = simulation_FHN(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["parameters_dc"]
                                       , param["initial_conditions"])
            delta_col.append(delta)
            rho_col.append(rho)
            num_peaks_col.append(len(list_peaks))
            printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    data = pd.DataFrame({'delta': delta_col,'rho': rho_col,'num_peaks': num_peaks_col})
    data['frequency'] = frequency
    data['current'] = current
    if not os.path.exists("DataFHN"):
        os.makedirs("DataFHN")
    data.to_csv(f'DataFHN/f{frequency}_I{current}.csv')
    return data

def plot_single_simulation_FHN(current, frequency, delta, rho):

    param= parameters_simulation_FHN(current,frequency,delta,rho)
    df_times, _ = simulation_FHN(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["parameters_dc"]
                                       , param["initial_conditions"])
    marker_flag =0
    delay = [30]
    
    val = 50
    for vm in df_times['environment.t']:
        if vm >=50:
            val=vm
            break

    x = df_times['environment.t']
    w = df_times['membrane.W']
    y = df_times['membrane.V_prom']
    #plt.figure(figsize = (16,5))
    plt.plot(x,y,color='royalblue')

    #makes a mark at time  = val
    for i in delay:
        if marker_flag == 0:
            plt.plot(df_times['environment.t'][i+df_times['environment.t'].index(val)], round(df_times['membrane.V_prom'][i+df_times['environment.t'].index(val)],2),
                marker='x',label = str(round(df_times['membrane.V_prom'][i+df_times['environment.t'].index(val)],2))+' [mV]',color='green')
        else :
            plt.plot(df_times['environment.t'][i+df_times['environment.t'].index(val)], round(df_times['membrane.V_prom'][i+df_times['environment.t'].index(val)],2),
                marker='o',label = str(round(df_times['membrane.V_prom'][i+df_times['environment.t'].index(val)],2))+' [mV]', color='green')

    plt.title(rf'Action Potential - $\delta$:{delta} - $\rho$:{rho}')
    plt.xlabel('time [ms]')
    plt.ylabel('voltage [mV]')
    return

def main():
    # Execute experiment 2: amplitude of HF term vs slope for the FHN model
    
    #range_delta = [0.01,0.195,0.005]
    #range_rho = [0.1,1.55,0.05]
    range_delta = [0.01,0.2,0.05]
    range_rho = [0.1,0.75,0.1]

    # Generate simulation  
    for frequency in [1,5,10,15,20]:
        for current in [0.15,0.3,0.5,0.5]:
            output = run_experiment2(current,frequency,range_delta,range_rho) 
    #for debugging
    #plot_single_simulation_FHN(current, frequency, delta = 0.11, rho = 0.5)
    return

main()