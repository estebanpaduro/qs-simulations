import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sym
import seaborn as sns
from scipy.signal import find_peaks
import myokit
import gc

## Experiment 1: Simulation for FHN model with input current of form ()
def simulation_FHN_exp1(model, protocol, func, time, threshold , parameters_hf, initial):
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
    ST = mod.get('ST')

    ## Oscillatory Function ##
    #--------Parameters------#
    for i in parameters_hf:
        var = ST.add_variable(i)
        var.set_rhs(parameters_hf[i])
        
    stim_hf = ST.add_variable('stim_hf')
    stim_hf.set_rhs(func[0])
    
    ## Final Stimulus ##
    stim = mod.get('ST.stim') #units [uA/cm^2]
    stim.set_rhs('stim_hf')

    ###########################
    membrane = mod.get('membrane')
    V_prom = membrane.add_variable('V_prom')
    V_prom.set_rhs('piecewise(environment.t<ST.t1,V-ST.lambda*environment.t*ST.rho*sin(ST.omega*environment.t),V-ST.lambda*ST.t1*ST.rho*sin(ST.omega*environment.t))')
    
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
     
    ## Run Simulation ##
    df_times = s.run(time)

    ## Find Peaks ##
    list_peaks, _ = find_peaks(df_times['membrane.V_prom'], height = threshold,
                               prominence = (1,None)) 
    return df_times, list_peaks

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
      
def parameters_simulation_FHN_exp1(frequency,rho,lambda_):
    """Parameters simulation with base model FHN.mmt
    Parameters
    ----------
    lambda_     - Required : slope HF function
    rho         - Required : amplitude high frequency
    frequency   - Required : frequency oscillatory function
    """
    
    model = 'FHN.mmt'   # In this line you can change base model
    protocol = None     # In this line you can change base protocol
    ### Stimulus ###
    parameters_hf = {'rho': rho,
                     'pi':np.pi,
                     'omega':frequency,
                     'lambda':lambda_,
                     'g1':'lambda*t*rho*omega*cos(omega*t)',
                     't1':'(1/lambda)',
                     'g2':'(lambda*t1)*rho*omega*cos(omega*t)'}
    oscillatory_function = 'piecewise(t<t1,g1,g2)'
    functions = [oscillatory_function]
    time = 300 
    threshold = 1   #This threshold is for FHN model
    ## this condition ensures starting from a stable state ##
    initial_conditions = [-1.12,-0.65]
    return {"model": model
           , "protocol": protocol
           , "parameters_hf": parameters_hf
           , "oscillatory_function": oscillatory_function
           , "functions": functions
           , "time": time
           , "threshold": threshold
           , "initial_conditions": initial_conditions
           }

def plot_single_simulation_FHN_exp1(frequency, rho,lambda_):

    param= parameters_simulation_FHN_exp1(frequency,rho,lambda_)
    df_times, list_peaks = simulation_FHN_exp1(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["initial_conditions"])

    x = df_times['environment.t']
    w = df_times['membrane.W']
    y = df_times['membrane.V_prom']
    z = df_times['ST.stim']

    fig = plt.figure(figsize = (16,10))
    fig, (ax1,ax2) =plt.subplots(2)
    ax1.plot(x,z,color='red')
    ax1.set_title('Stimulus')
    ax1.set_xlabel('time [ms]')
    ax1.set_ylabel('[uA]')
    ax1.set_xlim(0,param['time'])
    ax2.plot(x,y,color='royalblue')
    if len(list_peaks)!=0:
        for i in list_peaks:
            ax2.plot(x[i],y[i], marker='o',color='green', label = str(round(y[i],2)))
    ax2.set_title(rf'Action Potential - $\lambda$:{lambda_} - $\rho$:{rho}')
    ax2.set_xlabel('time [ms]')
    ax2.set_ylabel('voltage [mV]')
    ax2.set_xlim(0,param['time'])
    ax2.set_ylim(-3,3)
    ax2.legend()
    fig.tight_layout()
    fig.savefig(f'Figure/FHN_HF/{frequency}kHz_{lambda_}.png')
    plt.close(fig)
    
    return

def run_experiment1(frequency, range_rho,range_lambda):
    rho_col = []
    lambda_col =[]
    num_peaks_col = []
    
    #tracking progress 
    lambda_range = np.arange(range_lambda[0],range_lambda[1],range_lambda[2])
    rho_range = np.arange(range_rho[0],range_rho[1],range_rho[2])
    total = len(lambda_range)*len(rho_range)
    progress = 0
    print(f'Simulation for FHN_HF with frequency: {frequency}')
    printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for lamb in lambda_range:
        for rho in rho_range:
            #print('{:.1%}'.format(1.0*progress/total))
            progress = progress+1
            param= parameters_simulation_FHN_exp1(frequency,lamb,rho)
            _, list_peaks = simulation_FHN_exp1(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters_hf"]
                                       , param["initial_conditions"])

            lambda_col.append(lamb)
            rho_col.append(rho)
            num_peaks_col.append(len(list_peaks))
            printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    data = pd.DataFrame({'lambda': lambda_col,'rho': rho_col,'num_peaks': num_peaks_col})
    data['frequency'] = frequency
    data.to_csv(f'DataFHN_HF/f{frequency}.csv')
    return data