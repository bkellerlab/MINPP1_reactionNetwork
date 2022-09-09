# This file contains the fit functions and the analytic derivatives of the fit functions for all progress
# curves of the InsP5[2OH] dephosphorylation pathway.

import numpy as np

########################################################################################################
#### D A T A   I N T E R P O L A T I O N   #############################################################
########################################################################################################
   
# --- all fit function -------------------------------------------------------------

def fit_101111(t, k1=0.100, a0=182.834):
    return a0*np.exp(-k1*t)

def fit_101101(t, a=0.734, b=0.459, k=0.0008):
    return a*t**b*np.exp(-t*k)

def fit_100111(t, a=7.508, b=1.322, k=0.048):
    return a*t**b*np.exp(-t*k)

def fit_100110(t, a=0.00684, b=0.000262, c0=167.378):
    return ((a*c0)/(b-a))*(np.exp(-a*t) - np.exp(-b*t))

def fit_100101(t, a=0.181, b=0.937, k=0.003):
    return a*t**b*np.exp(-t*k)

def fit_100100(t, a=-0.0529, b=0.000317, S=70.446):
    return S - (S - a)*np.exp(-b*t)

def fit_100010(t, a=-0.289, b=0.000305, S=17.523):
    return S - (S - a)*np.exp(-b*t)

def fit_000110(t, a=-0.582, b=0.000222, S=71.352):
    return S - (S - a)*np.exp(-b*t)

def fit_100000(t, a=0.0000000789):
    return a*t**2

def fit_000100(t, a=0.0000000706):
    return a*t**2



# --- first derivative of all fitfunctions ----------------------------------------

def fit_derivative_101111(t, k1=0.100, a0=182.834):
    return -a0*k1*np.exp(-k1*t)

def fit_derivative_101101(t, a=0.734, b=0.459, k=0.0008):
    return a*np.exp(-k*t)*(b*t**(b-1) - k*t**b)

def fit_derivative_100111(t, a=7.508, b=1.322, k=0.048):
    return a*np.exp(-k*t)*(b*t**(b-1) - k*t**b)

def fit_derivative_100110(t, a=0.00684, b=0.000262, c0=167.378):
    return ((a*c0)/(b-a))*(-a*np.exp(-a*t) + b*np.exp(-b*t))

def fit_derivative_100101(t, a=0.181, b=0.937, k=0.003):
    return a*np.exp(-k*t)*(b*t**(b-1) - k*t**b)

def fit_derivative_100100(t, a=-0.0529, b=0.000317, S=70.446):
    return (S - a)*b*np.exp(-b*t)

def fit_derivative_100010(t, a=-0.289, b=0.000305, S=17.523):
    return (S - a)*b*np.exp(-b*t)

def fit_derivative_000110(t, a=-0.582, b=0.000222, S=71.352):
    return (S - a)*b*np.exp(-b*t)

def fit_derivative_100000(t, a=0.0000000789):
    return 2*a*t

def fit_derivative_000100(t, a=0.0000000706):
    return 2*a*t


# --- combine time evolution of all species in a matrix ----------------------------

def calculate_exp_data_interpolation(maxstep, tau):
    # order of species in concentration array
    header = ["101111", "101101", "100111", "100110", "100101", "100100", "100010", "000110", "100000", "000100"]
    # index of species in concentration array
    index = {}
    for i in range(0,len(header)):
        index[header[i]] = i
     
    # time evolution
    steps = np.arange(1, maxstep, tau)
        
    # interpolated concentration of all species
    conc = np.zeros((10,len(steps)))
    conc[0,:] = fit_101111(steps)
    conc[1,:] = fit_101101(steps)
    conc[2,:] = fit_100111(steps)
    conc[3,:] = fit_100110(steps)
    conc[4,:] = fit_100101(steps)
    conc[5,:] = fit_100100(steps)
    conc[6,:] = fit_100010(steps)
    conc[7,:] = fit_000110(steps)
    conc[8,:] = fit_100000(steps)
    conc[9,:] = fit_000100(steps)
    # set all negative entries to zero
    conc[conc < 0] = 0
    return index, steps, conc

def calculate_derivative_exp_data_interpolation(maxstep, tau):
    # time evolution
    steps = np.arange(1, maxstep, tau)
        
    # interpolated concentration of all species
    conc = np.zeros((10,len(steps)))
    conc[0,:] = fit_derivative_101111(steps)
    conc[1,:] = fit_derivative_101101(steps)
    conc[2,:] = fit_derivative_100111(steps)
    conc[3,:] = fit_derivative_100110(steps)
    conc[4,:] = fit_derivative_100101(steps)
    conc[5,:] = fit_derivative_100100(steps)
    conc[6,:] = fit_derivative_100010(steps)
    conc[7,:] = fit_derivative_000110(steps)
    conc[8,:] = fit_derivative_100000(steps)
    conc[9,:] = fit_derivative_000100(steps)
    return steps, conc   

def return_plot_color():
    # color definitions
    cyan = "#1f77b4"
    orange = "#ff7f0e"
    green = "#2ca02c"
    red = "#d62728"
    violet = "#9467bd"
    brown = "#8c564b"
    rose = "#e377c2"
    gray = "#7f7f7f"
    blue = 'b'
    black = 'k'
    colr = np.array([cyan, black, orange, red, green, violet, blue, brown, rose, gray])
    return colr


#------------------------------------------------------------------------------------


