# -*- coding: utf-8 -*-
"""
This programm solves the transcendent equation appearing to determine
an estimation value of the mean lifetime of a myon after applying the
max-log-likelihood-method on an exponential decaying law.

Input: N measured samples {t_1,...,t_N}
Output: Estimation value for the mean lifetime
"""

import numpy as np
from matplotlib import pyplot as plt

        
def onMouseClick(event):
    """ Start the dynamic plot of the location-distribution."""    
    mode = plt.get_current_fig_manager().toolbar.mode
    # plot should just be started by clicking the upper left window
    if  event.button == 1 and event.inaxes == axHist and mode == '':            
        global x_i   
        print "x_i = ", x_i
        x_ip1 =  getNextIteration(x_i, t, T, N)
        print "x_i+1 = ", x_ip1
        axHist.plot(x_i,0, ls='', marker='o')
        plt.draw()
        x_i = x_ip1
        

def f(x_i, t, T, N):
    """
        contractive function which fixpoint will be determined
    """
    funcVal = 0.25 + np.e**(-1./x_i)/(1- np.e**(-1./x_i)) - x_i
    print "f(x_i) = ", funcVal
    return  funcVal
    
def Df(x_i):
    """
        derivative of f with respect to x
    """
    funcVal =  np.e**(1./x_i)/(x_i**2 * (np.e**(1./x_i) - 1)**2) - 1
    return np.e**(1./x_i)/(x_i**2 * (np.e**(1./x_i) - 1)**2) - 1

def getNextIteration(x_i, t, T, N):
    """ This function calculates the next fixpoint-iteration."""
    return x_i - f(x_i, t, T, N)/Df(x_i)

    
T = 4.167                           # maximal measurable time in us
t = np.linspace(1,T,10)             # measured samples for lifetime
N = len(t)                          # sample size
x_0 = 0.5                          # start value of the iteration 
global x_i
x_i =  x_0  
x = np.linspace(10e-9,1,100)            # sample points for plotting f
print 1./N * np.sum(t/T)

fig = plt.figure("Pressure measurement", figsize=(15,8))
axHist = fig.add_subplot(111)
axHist.set_xlabel(r'$x$')
axHist.set_ylabel(r'$f(x)$')
axHist.set_xlim([0, 1])
axHist.set_title('Newton-Method')
axHist.plot(x, f(x, t, T, N))
axHist.plot(x, np.zeros(len(x)))

plt.connect('button_press_event', onMouseClick)
plt.show()

