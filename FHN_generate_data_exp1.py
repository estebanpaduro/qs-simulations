import numpy as np
import pandas as pd
import myokit
from scipy.signal import find_peaks
from scipy.optimize import fsolve


## Experiment 1: Simulation for FHN model with input current of form (3)
def simulation_FHN_exp1(model, protocol, func, time, threshold , parameters, initial, beta):
    """Detect action potential or not for a set parameters.

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
    EXT = mod.get('EXT')

    ## High Frequency term ##
    #--------Parameters------#
    for i in parameters:
        if True == EXT.can_add_variable(i):
            var = EXT.add_variable(i)
            var.set_rhs(parameters[i])
        else:
            var = mod.get('EXT.'+i)
            var.set_rhs(parameters[i])
        
    stim_hf = EXT.add_variable('stim_hf')
    stim_hf.set_rhs(func[0])
    
    ## Final Stimulus ##
    stim = mod.get('EXT.stim') #units [uA/cm^2]
    stim.set_rhs('stim_hf')

    b = mod.get('membrane.beta')
    b.set_rhs(beta)

    ###########################
    membrane = mod.get('membrane')
    V_prom = membrane.add_variable('V_prom')
    V_prom.set_rhs('piecewise(t<(1/EXT.lambda),V-EXT.lambda*t*EXT.rho*sin(EXT.omega*t),V-EXT.rho*sin(EXT.omega*t))')
    
    ## Initial Conditions ##
    v = mod.get('membrane.V')
    v.set_state_value(initial[0])
    w = mod.get('membrane.W')
    w.set_state_value(initial[1])

    v_AV = mod.get('membrane.V_A1')
    v_AV.set_state_value(initial[0])
    w_AV = mod.get('membrane.W_A1')
    w_AV.set_state_value(initial[1])
    
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
    list_peaks, _ = find_peaks(df_times['membrane.V_A1'], height = threshold,
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
      
def parameters_simulation_FHN_exp1(beta,rho,lambda_,frequency=10):
    """Parameters simulation with base model FHN.mmt
    Parameters
    ----------
    beta        - Required : value of beta
    lambda_     - Required : slope HF function
    rho         - Required : amplitude high frequency
    frequency   - Optional : frequency oscillatory function
    """
    
    model = 'FHN_exp1.mmt'   # In this line you can change base model
    protocol = None     # In this line you can change base protocol
    ### Stimulus ###
    parameters = {'rho': rho,
                     'pi':np.pi,
                     'omega':frequency,
                     'lambda':lambda_,
                     'g1':'lambda*t*rho*omega*cos(omega*t)',
                     't1':'(1/lambda)',
                     'g2':'rho*omega*cos(omega*t)'}
    high_freq_term = 'piecewise(t<t1,g1,g2)'
    functions = [high_freq_term]
    time = 200 
    threshold = 1   #This threshold is for FHN model
    ## this condition ensures starting from a stable state ##
    def system_P0(variables):
        v0, w0 = variables
        gamma = 0.5
        bet = beta
        eq_1 = v0-v0**3/3-w0
        eq_2 = v0-gamma*w0+bet
        return [eq_1,eq_2]
    cond_init = [0,0]
    v0, w0 = fsolve(system_P0,cond_init)
    initial_conditions = [v0,w0]
    return {"model": model
           , "protocol": protocol
           , "parameters": parameters
           , "oscillatory_function": high_freq_term
           , "functions": functions
           , "time": time
           , "threshold": threshold
           , "initial_conditions": initial_conditions
           , 'beta':beta
           }

def run_experiment1(beta,range_rho,range_lambda):
    rho_col = []
    lambda_col =[]
    num_peaks_col = []
    
    #tracking progress 
    lambda_range = np.arange(range_lambda[0],range_lambda[1],range_lambda[2])
    rho_range = np.arange(range_rho[0],range_rho[1],range_rho[2])
    total = len(lambda_range)*len(rho_range)
    progress = 0
    print(f'Simulation for FHN_HF with beta {beta}:')
    printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    for lamb_log in lambda_range:
        lamb = 10**lamb_log
        for rho_log in rho_range:
            rho = 10**rho_log
            #rho = ((r.round(4))/(lamb.round(4))).round(4)

            #print('{:.1%}'.format(1.0*progress/total))
            progress = progress+1
            param= parameters_simulation_FHN_exp1(beta,rho,lamb)
            _, list_peaks = simulation_FHN_exp1(param["model"]
                                       , param["protocol"]
                                       , param["functions"]
                                       , param["time"]
                                       , param["threshold"]
                                       , param["parameters"]
                                       , param["initial_conditions"]
                                       , param["beta"])

            lambda_col.append(lamb_log)
            rho_col.append(rho_log)
            num_peaks_col.append(len(list_peaks))
            printProgressBar(progress, total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    data = pd.DataFrame({'lambda_ln': lambda_col,'rho_ln': rho_col,'num_peaks': num_peaks_col})
    data.to_csv(f'exp1_FHN_beta_{beta}_log10.csv')
    return data
