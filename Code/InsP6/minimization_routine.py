# This file contains the full code for the minimization routine for the InsP6 dephosphorylation.
# The analysis is based on the fitted experimental data contained in fitted_exp_data.py

import matplotlib.pyplot as plt
import numpy as np
import fitted_exp_data as func
from scipy import linalg
from scipy import optimize



#####################################################################
### Error function that has to be minimized #########################
#####################################################################

def calculate_error_for_one_tau(k01, k02, k13, k14, k24, k25, k36, k37, k47, k48, k58, k69, k79, k89, k810, k911, k1011, a):
    ''' Calculate Delta(n tau) as in equation B.9
    '''
    # concentration for time increment a (fits to experimental data)
    x = phi[:,a]
    # time derivative of concentration for time increment a (calculated analytically from fits)
    y = phidot[:,a]
  
    # assign experimental values for concentration
    x0 = x[0]
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    x4 = x[4]
    x5 = x[5]
    x6 = x[6]
    x7 = x[7]
    x8 = x[8]
    x9 = x[9]
    x10 = x[10]

    # experimental values for time derivative of concentration
    yexp0 = y[0]
    yexp1 = y[1]
    yexp2 = y[2]
    yexp3 = y[3]
    yexp4 = y[4]
    yexp5 = y[5]
    yexp6 = y[6]
    yexp7 = y[7]
    yexp8 = y[8]
    yexp9 = y[9]
    yexp10 = y[10]
    yexp11 = y[11]
    
    # compute diagonal elements as
    k00 = - (k01 + k02)
    k11 = - (k13 + k14)
    k22 = - (k24 + k25)
    k33 = - (k36 + k37)
    k44 = - (k47 + k48)
    k55 = - (k58)
    k66 = - (k69)
    k77 = - (k79)
    k88 = - (k89 + k810)
    k99 = - (k911)
    k1010 = - (k1011)
        
    # system of differential equations
    y0 = k00*x0
    y1 = k01*x0 + k11*x1
    y2 = k02*x0 + k22*x2
    y3 = k13*x1 + k33*x3
    y4 = k14*x1 + k24*x2 + k44*x4
    y5 = k25*x2 + k55*x5
    y6 = k36*x3 + k66*x6
    y7 = k37*x3 + k47*x4 + k77*x7
    y8 = k48*x4 + k58*x5 + k88*x8
    y9 = k69*x6 + k79*x7 + k89*x8 + k99*x9
    y10 = k810*x8 + k1010*x10
    y11 = k911*x9 + k1011*x10
           
    # define error function as square (equation B.9)
    Delta = (y0 - yexp0)**2 + (y1 - yexp1)**2 + (y2 - yexp2)**2 + (y3 - yexp3)**2 + (y4 - yexp4)**2 + (y5 - yexp5)**2 + (y6 - yexp6)**2 + (y7 - yexp7)**2 + (y8 - yexp8)**2 + (y9 - yexp9)**2 + (y10 - yexp10)**2 + (y11 - yexp11)**2    
    return Delta

def calculate_overall_error(k):
    ''' Calculate Delta = Delta(n tau) + Delta(m tau) + ... + Delta(k tau)
    '''
    
    # assign array elements to crresponding rates
    k02 = 7.4e-04  # constraint
    k01 = 0.000932 - k02  # constraint
    k13 = k[0]
    k14 = k[1]
    k24 = k[2] 
    k25 = 0.00093 - k24  # constraint
    k36 = k[3]
    k37 = k[4]
    k47 = k[5]
    k48 = k[6]
    k58 = k[7]
    k69 = k[8]
    k79 = k[9]
    k89 = k[10]
    k810 = k[11]
    k911 = k[12]
    k1011 = k[13]
    
    
    # define overall error function
    time = np.arange(0, maxstep, tau)
    result = 0
    for i in range(0, len(time)-1):
        result = result + calculate_error_for_one_tau(k01, k02, k13, k14, k24, k25, k36, k37, k47, k48, k58, k69, k79, k89, k810, k911, k1011, time[i])
    return result



####################################################################
### Analysis #######################################################
####################################################################

# maximal time value in min and time increment tau in min
maxstep = 10081
tau = 1

# load fitted progress curves and time derivative of fitted progress curves
index, steps, phi = func.calculate_exp_data_interpolation(maxstep, tau)
steps, phidot = func.calculate_derivative_exp_data_interpolation(maxstep, tau)

# initial guess for all rates
initial_guess = np.array([2.9e-03, 1.0e-05, 3.0e-04, 5.0e-04, 3.7e-04, 1.4e-03, 1.5e-04, 4.2e-04, 2.0e-04, 2.0e-04, 3.9e-03, 1.0e-05, 1.0e-03, 1.0e-03])

# upper and lower bound for each rate
boundss = ((0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.001,1), (0.00001,1), (0.0001,1))

# minimize error functions
minimized = optimize.minimize(fun=calculate_overall_error, x0=initial_guess, bounds = boundss)



####################################################################
### Backcompute progress curves from computed rates ################
####################################################################

# --- compute time evolution of density using previously computed rates -------

# prepare array for concentrations
result = np.zeros((12,len(steps)))
result[0,0] = 175

# set-up rate matrix
K = np.zeros((12, 12))
K[0,2] = 7.4e-04
K[0,1] = 0.000932 - K[0,2]
K[1,3] = minimized.x[0]
K[1,4] = minimized.x[1]
K[2,4] = minimized.x[2]
K[2,5] = 0.00093 - K[2,4]
K[3,6] = minimized.x[3]
K[3,7] = minimized.x[4]
K[4,7] = minimized.x[5]
K[4,8] = minimized.x[6]
K[5,8] = minimized.x[7]
K[6,9] = minimized.x[8]
K[7,9] = minimized.x[9]
K[8,9] = minimized.x[10]
K[8,10] = minimized.x[11]
K[9,11] = minimized.x[12]
K[10,11] = minimized.x[13]

# compute diagonal elements of rate matrix
for i in range(0, 12):
    K[i,i] = -np.sum(K[i,:])

# compute propagator
P = linalg.expm(K*tau)

# propagate concentrations in time
for i in range(0,len(steps)-1):
    result[:,i+1] = np.dot(result[:,i],P)
    
# get individual plot color for each species
colr = func.return_plot_color()
# get labels for all species labels
labs = list(index.keys())



##################################################################
### compare results to original fit functions ####################
##################################################################


fig = plt.figure()
for i in range(0,12):
    plt.plot(steps, phi[i,:], color=colr[i], label=labs[i])
    plt.plot(steps, result[i,:], '--', color=colr[i])

plt.xlabel('time in min')
plt.ylabel('concentration in $\\mu$M')
plt.legend(fontsize=16, framealpha=1, loc=1, bbox_to_anchor=(1.4,1.13))







