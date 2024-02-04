
import numpy as np
import matplotlib.pyplot as plt
import os

#check for output directory
if not os.path.exists("Figures"):
    os.makedirs("Figures")

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
plt.savefig('Figures/Fig1.png', bbox_inches='tight', pad_inches=0.5)
plt.close()
############################