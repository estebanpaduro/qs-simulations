import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import myokit
import os


#For the simulation

def simulation_HH(model, protocol, functions, time, threshold, 
               parameters_hf, parameters_dc, initial_conditions):
    """Detect blocking or not for a set parameters.
    Parameters
    ----------
    model : base model of myokit in .mmt
    protocol : in case the external stimulus needs a protocol
        {[level, start, lenght, period, multiplier], None}, 
        optional 
    functions : list function for external stimulus
        [oscillatory function, constant function]
    time : time simulation
    thresh : 
    max_peaks : integer
    parameters_hf : parameters oscillatory function
        {{'parameter':value},{}}, if the parameter depends on other 
        value so format of value is str
    parameters_dc : parameters constant function
        {{'parameter':value},{}}, if the parameter depends on other 
        value so format of value is str
    initial_conditions: initial conditions
            [m,n,h,voltage]
    
    Returns block, data, num_peaks
    -------
    block --> {0,1} : 0 if there is blocking and 1 if there is not blocking
    data --> data simulation 
    num_peaks --> number of peaks from this simulation
    """
    
    mod = myokit.load_model(model)
    
    prot = myokit.Protocol()
    if protocol is not None:
        prot.schedule(protocol[0], protocol[1], protocol[2], protocol[3], protocol[4])
        
    
    ##### Block Stimulus #####  
    ST = mod.get('ST')
    ## oscillatory function ##
    #--------Parameters------#
    for i in parameters_hf:
        var = ST.add_variable(i)
        var.set_rhs(parameters_hf[i])   
    stim_hf = ST.add_variable('stim_hf')
    stim_hf.set_rhs(functions[0])
    ## constant function ##
    #------Parameters-----#
    for i in parameters_dc:
        var = ST.add_variable(i)
        var.set_rhs(parameters_dc[i])
    
    stim_dc = ST.add_variable('stim_dc')
    stim_dc.set_rhs(functions[1])
    ## final stimulus ##
    stim = mod.get('ST.stim') #units [uA/cm^2]
    stim.set_rhs('stim_dc+stim_hf')
    ###########################
    
    ## Initial Conditions ##

    m = mod.get('INA.m')
    m.set_state_value(initial_conditions[0])
    n = mod.get('IK.n')
    n.set_state_value(initial_conditions[1])
    h = mod.get('INA.h')
    h.set_state_value(initial_conditions[2])
    V = mod.get('membrane.V')
    V.set_state_value(initial_conditions[3])
    
    ## Simulation ##  
    s = myokit.Simulation(mod,prot,path='my-sim.zip')
    
    ### step default ###
    s.set_max_step_size(dtmax=0.005)
    s.set_min_step_size(dtmin=None)
    s.set_tolerance( abs_tol=1e-06 , rel_tol=0.0001 ) 
    ######################## 

    ## Run simulation ##
    df_times = s.run(time)
    ## Find Peaks ##
    dist=(1/(50/s.last_number_of_steps()))
    list_peaks, _ = find_peaks(df_times['membrane.V'],
                               height = threshold, distance = dist)

    return df_times, list_peaks

def parameters_simulation_HH(current, frequency, delta, rho):
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
    data: data simulation in myokit
    num_peaks: peaks 
    """
    """Example
    parameters_simulation(5,0.15,80,4,True)
    """
    model = 'HH.mmt'
    protocol = None
    parameters_hf = {'rho': rho,
                           'f': frequency,
                           'pi':np.pi,
                           'omega_rad':'f*2*pi'}
    oscillatory_function = 'rho*f*1000*cos(omega_rad*t)'
    parameters_dc = {'m': delta,
                           'C':current,
                           't0':'(1/m)+50',
                           'f1':0,
                           'f2': 'C*(t-50)*m',
                           'f3':'piecewise(t<50,f1,f2)'}
    constant_function = 'piecewise(t<t0,f3,C)'
    
    functions = [oscillatory_function,constant_function]
    time = 200
    threshold = 0
    ## this condition ensures starting from a stable state ##
    initial_conditions = [0.0529,0.3176,0.5960,-65.005]
    
    return {"model": model
            , "protocol": protocol
            , "parameters_hf": parameters_hf
            , "oscillatory_function" : oscillatory_function
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
 
def run_experiment1(current, frequency,range_delta,range_rho):
    """
    Parameters
    ----------
    
    Returns 
    -------

    """
    
    delta_col=[]
    rho_col=[]
    num_peaks_col=[]
    
    #tracking progress 
    delta_range = np.arange(range_delta[0],range_delta[1],range_delta[2])
    rho_range = np.arange(range_rho[0],range_rho[1],range_rho[2])
    total = len(delta_range)*len(rho_range)
    progress = 0
    print(f'Simulation for HH with current: {current}, and frequency: {frequency}')
    printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for delta in delta_range:
        for rho in rho_range:
            progress = progress+1
            param = parameters_simulation_HH(current, frequency, delta, rho)
            _, list_peaks = simulation_HH(param["model"]
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

    data = pd.DataFrame({'delta':delta_col,'rho':rho_col,'num_peaks':num_peaks_col})
    data['frequency'] = frequency
    data['current'] = current
    if not os.path.exists("DataHH"):
        os.makedirs("DataHH")
    data.to_csv(f'DataHH/f{frequency}_I{current}.csv')
    return data

def plot_single_simulation_HH(current, frequency, delta, rho):

    param = parameters_simulation_HH(current, frequency, delta, rho)
    df_times, list_peaks = simulation_HH(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["parameters_dc"]
                                       , param["initial_conditions"])
    MAXPEAKS = 0
    #plot#
    x = df_times['environment.t']
    y = df_times['membrane.V']
    plt.figure(figsize = (16,5))
    plt.plot(x,y,color = 'green')
    for i in list_peaks:
        plt.plot(df_times['environment.t'][i], df_times['membrane.V'][i],
                 label = str(round(df_times['membrane.V'][i],2))+' [mV]', marker='o', color = 'red')
    plt.title('Action Potential')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.legend()
    plt.show()
    #if func[1] != '0':
    plot_data('environment.t','ST.stim_dc','time','[uA/cm^2]','stimulus','Constant Stimulus',df_times)
       #if func[0] != '0':
       #    plot_data('environment.t','ST.stim_hf','time','[uA/cm^2]','stimulus','Oscillatory Stimulus',data)
        
    return int(len(list_peaks) > MAXPEAKS), df_times, len(list_peaks)

def main():
    # Execute experiment 1: amplitude of HF term vs slope for the HH model

    #parameter mesh
    #range_delta = [0.0025,0.1275,0.0025]
    #range_rho = [0.01,0.3,0.01]
    
    #smaller range for testing
    range_delta = [0.0025,0.1275,0.01]
    range_rho = [0.01,0.3,0.05]

    for frequency in [3,5,10]:
        for current in [40, 80, 120]:
            output = run_experiment1(current,frequency,range_delta,range_rho) 
    
    #for debugging
    #current = 80
    #frequency = 5
    #plot_single_simulation_HH(current, frequency, delta = 0.05, rho = 0.06)
    return


main()