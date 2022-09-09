# This file contains the fit functions and the analytic derivatives of the fit functions for all progress
# curves of the InsP6 dephosphorylation pathway.

import numpy as np

########################################################################################################
#### D A T A   I N T E R P O L A T I O N   #############################################################
########################################################################################################
   
# --- all fit function -------------------------------------------------------------

def fit_111111(t, k=0.00092049, a0=175):
    return a0*np.exp(-k*t)

def fit_111011(t, a=6.2e-18, b=6.697, k=0.00460):
    return a*t**b*np.exp(-k*t)

def fit_110111(t, k1=0.000741, k2=0.00093705, c0=175):
    return ((k1*c0)/(k2-k1))*(np.exp(-k1*t) - np.exp(-k2*t))

def fit_111001(t, a=45.149, mean=2057.099, var=709.875):
    return a*np.exp(-(t-mean)**2/var**2)

def fit_110011(t, a=0.00000211, b=2.603, k=0.00272):
    return a*t**b*np.exp(-k*t) 

def fit_110110(t, a=55.0451, mean=2849.4626, var=1045.3938):
    return a*np.exp(-(t-mean)**2/var**2)

def fit_111000(t, a=61.897, b=4.578, c=0.000768, d=3133.404):
    return a*np.exp(-b*(1 - np.exp(-c*(t-d)))**2)

def fit_110001(t, a=0.0973, b=45.655, c=0.00249, d=7.480, e=0.000143):
    return d/(a + b*np.exp(-c*t))*np.exp(-e*t)

def fit_110010(t, a=9.9097, mean=3338.0448, var=626.4720):
    return a*np.exp(-(t-mean)**2/var**2)

def fit_110000(t, a=0.00004617, b=1394.199, c=0.004425, d=0.005908):
    return d/(a + b*np.exp(-c*t))

def fit_010010(t, a=11.970, b=0.0000240, c=-0.00497, d=0.000468, e=-0.00486):
    return d/(a + b*np.exp(-c*t))*np.exp(-e*t)

def fit_010000(t, a=-8.741, b=-3609.113, c=0.000978, d=-168.344):
    return d/(a + b*np.exp(-c*t))



# --- first derivative of all fitfunctions ----------------------------------------

def fit_derivative_111111(t, k=0.00092049, a0=175):
    return -k*a0*np.exp(-k*t)

def fit_derivative_111011(t, a=6.2e-18, b=6.697, k=0.00460):
    return a*np.exp(-k*t)*(b*t**(b-1) - k*t**b)

def fit_derivative_110111(t, k1=0.000741, k2=0.00093705, c0=175):
    return ((k1*c0)/(k2-k1))*(-k1*np.exp(-k1*t) + k2*np.exp(-k2*t))

def fit_derivative_111001(t, a=45.149, mean=2057.099, var=709.875):
    return -a*np.exp(-(t-mean)**2/var**2)*((2*(t-mean))/var**2)

def fit_derivative_110011(t, a=0.00000211, b=2.603, k=0.00272):
    return a*np.exp(-k*t)*(b*t**(b-1) - k*t**b)

def fit_derivative_110110(t, a=55.0451, mean=2849.4626, var=1045.3938):
    return -a*np.exp(-(t-mean)**2/var**2)*((2*(t-mean))/var**2)

def fit_derivative_111000(t, a=61.897, b=4.578, c=0.000768, d=3133.404):
    return -2*a*b*c*np.exp(-b*(1 - np.exp(-c*(t-d)))**2)*(1 - np.exp(-c*(t-d)))*np.exp(-c*(t-d))

def fit_derivative_110001(t, a=0.0973, b=45.655, c=0.00249, d=7.480, e=0.000143):
    return ((b*c*d*np.exp(-c*t))/(a + b*np.exp(-c*t)) - e*d)*(np.exp(-e*t)/(a + b*np.exp(-c*t)))

def fit_derivative_110010(t, a=9.9097, mean=3338.0448, var=626.4720):
    return -a*np.exp(-(t-mean)**2/var**2)*((2*(t-mean))/var**2)

def fit_derivative_110000(t, a=0.00004617, b=1394.199, c=0.004425, d=0.005908):
    return (b*c*d*np.exp(-c*t))/(a + b*np.exp(-c*t))**2

def fit_derivative_010010(t, a=11.970, b=0.0000240, c=-0.00497, d=0.000468, e=-0.00486):
    return ((b*c*d*np.exp(-c*t))/(a + b*np.exp(-c*t)) - e*d)*(np.exp(-e*t)/(a + b*np.exp(-c*t)))

def fit_derivative_010000(t, a=-8.741, b=-3609.113, c=0.000978, d=-168.344):
    return (b*c*d*np.exp(-c*t))/(a + b*np.exp(-c*t))**2


# --- combine time evolution of all species in a matrix ----------------------------

def calculate_exp_data_interpolation(maxstep, tau):
    # order of species in concentration array
    header = ["111111", "111011", "110111", "111001", "110011", "110110", "111000", "110001", "110010", "110000", "010010", "010000"]
    # index of species in concentration array
    index = {}
    for i in range(0,len(header)):
        index[header[i]] = i
     
    # time evolution
    steps = np.arange(1, maxstep, tau)
        
    # interpolated concentration of all species
    conc = np.zeros((12,len(steps)))
    conc[0,:] = fit_111111(steps)
    conc[1,:] = fit_111011(steps)
    conc[2,:] = fit_110111(steps)
    conc[3,:] = fit_111001(steps)
    conc[4,:] = fit_110011(steps)
    conc[5,:] = fit_110110(steps)
    conc[6,:] = fit_111000(steps)
    conc[7,:] = fit_110001(steps)
    conc[8,:] = fit_110010(steps)
    conc[9,:] = fit_110000(steps)
    conc[10,:] = fit_010010(steps)
    conc[11,:] = fit_010000(steps)

    # convert all nan to 0
    conc = np.nan_to_num(conc, 0)
    # set all negative entries to zero
    conc[conc < 0] = 0
    return index, steps, conc

def calculate_derivative_exp_data_interpolation(maxstep, tau):
    # time evolution
    steps = np.arange(1, maxstep, tau)
        
    # interpolated concentration of all species
    conc = np.zeros((12,len(steps)))
    conc[0,:] = fit_derivative_111111(steps)
    conc[1,:] = fit_derivative_111011(steps)
    conc[2,:] = fit_derivative_110111(steps)
    conc[3,:] = fit_derivative_111001(steps)
    conc[4,:] = fit_derivative_110011(steps)
    conc[5,:] = fit_derivative_110110(steps)
    conc[6,:] = fit_derivative_111000(steps)
    conc[7,:] = fit_derivative_110001(steps)
    conc[8,:] = fit_derivative_110010(steps)
    conc[9,:] = fit_derivative_110000(steps)
    conc[10,:] = fit_derivative_010010(steps)
    conc[11,:] = fit_derivative_010000(steps)
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
    yellow = "#e3d200"
    bright_green = "#66ff66"
    colr = np.array([cyan, black, orange, green, blue, red, violet, rose, brown, gray, bright_green, yellow])
    return colr


#------------------------------------------------------------------------------------


