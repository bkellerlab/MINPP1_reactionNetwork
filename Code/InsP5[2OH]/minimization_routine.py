# This file contains the full code for the minimization routine for the InsP5[2OH] dephosphorylation.
# The analysis is based on the fitted experimental data contained in fitted_exp_data.py


import matplotlib.pyplot as plt
import numpy as np
import fitted_exp_data as func
from scipy import linalg
from scipy import optimize

#####################################################################
### Error function that has to be minimized #########################
#####################################################################

def calculate_error_for_one_tau(k01, k02, k14, k23, k24, k35, k36, k37, k45, k58, k59, k68, k79, a):
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
    
    # compute diagonal elements as
    k00 = - (k01 + k02)
    k11 = - (k14)
    k22 = - (k23 + k24)
    k33 = - (k35 + k36 + k37)
    k44 = - (k45)
    k55 = - (k58 + k59)
    k66 = - (k68)
    k77 = - (k79)
        
    # system of differential equations
    y0 = k00*x0
    y1 = k01*x0 + k11*x1
    y2 = k02*x0 + k22*x2
    y3 = k23*x2 + k33*x3
    y4 = k14*x1 + k24*x2 + k44*x4
    y5 = k35*x3 + k45*x4 + k55*x5
    y6 = k36*x3 + k66*x6
    y7 = k37*x3 + k77*x7
    y8 = k58*x5 + k68*x6
    y9 = k59*x5 + k79*x7
           
    # define error function as square (equation B.9)
    Delta = (y0 - yexp0)**2 + (y1 - yexp1)**2 + (y2 - yexp2)**2 + (y3 - yexp3)**2 + (y4 - yexp4)**2 + (y5 - yexp5)**2 + (y6 - yexp6)**2 + (y7 - yexp7)**2 + (y8 - yexp8)**2 + (y9 - yexp9)**2    
    return Delta

def calculate_overall_error(k):
    ''' Calculate Delta = Delta(n tau) + Delta(m tau) + ... + Delta(k tau)
    '''
    
    # assign array elements to crresponding rates
    k01 = 0      # constraint
    k02 = k[0]
    k14 = 0      # constraint
    k23 = k[1]
    k24 = 0.001  # constraint
    k35 = k[2]
    k36 = k[3] 
    k37 = 0.000262 - k35 - k36  # constraint
    k45 = k[4]
    k58 = k[5]
    k59 = k[6]
    k68 = k[7]
    k79 = k[8]
    
    # define overall error function
    time = np.arange(0, maxstep, tau)
    result = 0
    for i in range(0, len(time)-1):
        result = result + calculate_error_for_one_tau(k01, k02, k14, k23, k24, k35, k36, k37, k45, k58, k59, k68, k79, time[i])
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
initial_guess = np.array([0.0976, 0.00684, 0.000109, 0.0000416, 0.000611, 0.00001, 0.0000108, 0.0000577, 0.0000024])

# upper and lower bound for each rate
boundss = ((0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1), (0.000001,1))

# minimize error functions
minimized = optimize.minimize(fun=calculate_overall_error, x0=initial_guess, bounds = boundss)



####################################################################
### Backcompute progress curves from computed rates ################
####################################################################

# prepare array for concentrations
result = np.zeros((10,len(steps)))
result[0,0] = 175

# set-up rate matrix
K = np.zeros((10, 10))
K[0,1] = 0
K[0,2] = minimized.x[0]
K[1,4] = 0
K[2,3] = minimized.x[1]
K[2,4] = 0.001
K[3,5] = minimized.x[2]
K[3,6] = minimized.x[3]
K[3,7] = 0.000262 - K[3,5] - K[3,6]
K[4,5] = minimized.x[4]
K[5,8] = minimized.x[5]
K[5,9] = minimized.x[6]
K[6,8] = minimized.x[7]
K[7,9] = minimized.x[8]

# compute diagonal elements of rate matrix
for i in range(0, 10):
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
for i in range(0,10):
    plt.plot(steps, phi[i,:], color=colr[i], label=labs[i])
    plt.plot(steps, result[i,:], '--', color=colr[i])

plt.xlabel('time in min')
plt.ylabel('concentration in $\\mu$M')
plt.legend(fontsize=16, framealpha=1, loc=1, bbox_to_anchor=(1.4,1.13))




















