import FHN_generate_data_exp2 as generate_exp2
import FHN_generate_data_exp1 as generate_exp1

### This block generate data experiment 1 ###
range_rho = [-0.9,-0.16,0.02]
range_lambda = [-1.8,0.06,0.03]
set_betas = [0.65,0.675,0.7,0.725,0.75,
             0.775,0.8,0.825,0.85,0.875,0.9]
for beta in set_betas:
    generate_exp1.run_experiment1(beta,range_rho,range_lambda)
#############################################

### This block generate data experiment 2 ###
range_delta = [-2.2,0.05,0.05]
range_rho = [-0.9,0.22,0.02]
set_currents = [0.5, 0.475, 0.45, 0.425, 0.4,
                0.375, 0.35, 0.325, 0.3,
                0.275, 0.25, 0.225, 0.2,
                0.175, 0.15, 0.125, 0.1]
for current in set_currents:
    generate_exp2.run_experiment2(current, range_delta, range_rho)
#############################################

