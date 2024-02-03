import numpy as np
import matplotlib.pyplot as plt
import FHN_report_from_data_exp1 as report_exp1
import FHN_report_from_data_exp2 as report_exp2

############################
""" Figure 1: Example I(t)"""

lambda_val = 0.05
rho = 1
omega = 1.5
delta = 0.05
I_0 = 0.5

def slope(t):
    if t <=1 /lambda_val:
        return t*lambda_val*rho*omega
    else:
        return
    
def dc_component(t):
    if t <= 40:
        return 0
    elif 40 < t <= 40+1/delta:
        return delta*(t-40)*I_0
    else:
        return I_0

def I_function(t):
    if t <= 1/lambda_val:
        return t*lambda_val*omega*rho*np.cos(omega*t)
    elif 1/lambda_val < t <= 40:
        return rho*omega*np.cos(omega*t)
    elif 40 < t <=40+ 1/delta:
        return rho*omega*np.cos(omega*t)+(delta*(t-40)*I_0)
    else:
        return I_0+rho*omega*np.cos(omega*t)


t_values = np.linspace(0, 100, 10000)
t_values_2 = np.linspace(0, 40, 10000)

I_values = [I_function(t) for t in t_values]
dc_values = [dc_component(t) for t in t_values]
slope_values = [slope(t) for t in t_values_2]

plt.plot(t_values, I_values, label=r'$I(t)$',color='royalblue')
plt.plot(t_values, dc_values, label=r'dc Component', color='red')
plt.plot(t_values_2, slope_values, color='black', linestyle='--')
plt.xlabel('t')
plt.ylabel('y(t)')
annotate_box = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.7)
plt.text(8, 0.98, r'$\lambda\rho\omega$', ha='center', va='center', fontsize=14, color='black', bbox=annotate_box)
plt.text(45, 0.3, r'$\delta I_0$', ha='center', va='center', fontsize=14, color='black', bbox=annotate_box)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Figures/figure1.png', bbox_inches='tight', pad_inches=0.5)
plt.close()
############################

""" Experiment 1"""

betas = [0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9]
markers = ['*','o','P','s']

### This block generate figure 2.a and 2.b using beta=0.75 ###
list_rho_plot = [0.4,0.6]
list_lambda_plot = [0.04,0.9]
report_exp1.region_lambda_rho_FHN_exp1(0.75, list_rho_plot,list_lambda_plot,markers,'figure_2a.png')
report_exp1.plot_simulations_FHN_exp1(0.75,list_rho_plot, list_lambda_plot,markers,'figure_2b.png')
##########################################################

### This block generate figure 2.c ###
report_exp1.curve(betas,'figure_2c.png')
######################################

""" Experiment 2"""

currents = [0.5,0.475,0.45,0.4,0.375,0.35,0.325,0.3,0.275,0.25,0.225,0.2,0.175,0.15,0.125,0.1]
markers = ['*','o','P','s']

### This block generate figure 3.a and 3.b using I_0 = 0.2###
list_rho_plot = [0.5,1.25]
list_delta_plot = [0.02,0.3]
report_exp2.region_delta_rho_FHN_exp2(0.2, list_rho_plot,list_delta_plot,markers,'figure_3a.png')
report_exp2.plot_simulations_FHN_exp2(0.2, list_rho_plot, list_delta_plot,markers,'figure_3b.png')
#############################################################

### This block generate figure 3.c and 3.d using I_0 = 0.4###
list_rho_plot = [0.5,1.25]
list_delta_plot = [0.3,0.02]
report_exp2.region_delta_rho_FHN_exp2(0.4, list_rho_plot,list_delta_plot,markers,'figure_3c.png')
report_exp2.plot_simulations_FHN_exp2(0.4, list_rho_plot, list_delta_plot,markers,'figure_3d.png')
#############################################################

### This block generate figure 3.e ###
report_exp2.curve(currents,name_fig='figure_3e.png')
######################################

