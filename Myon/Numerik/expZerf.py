# -*- coding: utf-8 -*-
"""
This programm solves the transcendent equation appearing to determine
an estimation value of the mean lifetime of a myon after applying the
max-log-likelihood-method on an exponential decaying law.

Input: N measured samples {t_1,...,t_N}
Output: Estimation value for the mean lifetime
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

def readFile(fn, startRow):
    """Read data from file fn. Returns an array with the data"""
    data = np.loadtxt(fn, skiprows = startRow)
    return data
    
        
def onMouseClick(event):
    """ Start the dynamic plot of the location-distribution."""    
    mode = plt.get_current_fig_manager().toolbar.mode
    # plot should just be started by clicking the upper left window
    if  event.button == 1 and event.inaxes == axHist and mode == '':            
        global x_i   
        print "x_i = ", x_i
        print "tau_i =", x_i * T, " us"
        x_ip1 =  getNextIteration(x_i, t, counts, T, N, cMin, cMax)
        axHist.plot(x_i,0, ls='', marker='o')
        axHist.set_title(r'Mittlere Lebensdauer von Myonen mit Exponentialfunktion $\tau =$' +
                         str(round(x_i * T,2)) + r'$\pm$' + str(round(sigma_tau,2)) + r' $\mu s$')
        plt.draw()
        x_i = x_ip1
        
        

def f(x_i, t, counts, T, N, cMin, cMax):
    """
        contractive function which fixpoint will be determined
    """
        
    funcVal = 1./(N*T) * np.sum(counts[cMin-1:cMax] * t[cMin-1:cMax]) \
              + np.e**(-1./x_i)/(1- np.e**(-1./x_i)) - x_i
    return  funcVal
    
def Df(x_i):
    """
        derivative of f with respect to x
    """
    funcVal =  np.e**(1./x_i)/(x_i**2 * (np.e**(1./x_i) - 1)**2) - 1
    return np.e**(1./x_i)/(x_i**2 * (np.e**(1./x_i) - 1)**2) - 1

def getNextIteration(x_i, t, counts, T, N, cMin, cMax):
    """ This function calculates the next fixpoint-iteration."""
    return x_i - f(x_i, t, counts, T, N, cMin, cMax)/Df(x_i)

    

fn = "mdat.txt"                     # file name for measuredata
startRow = 85                       # start row
data = readFile(fn, startRow)       # array with measuredata
cNr = data[:,0]                     # channelnumbers
cWdth = 1./24                       # channelwidth in us
t = cNr * cWdth                     # measureable lifetimes
cMin = 20                         # lower bond for channels
cMax = 175                          # upper bond for channels
counts =  data[:,1]                 # counts per channel
countErr = np.sqrt(counts)          # error on the counts

print t[8]



N = np.sum(counts[cMin-1:cMax])     # total number of counted events
print "N = ", N
T = t[cMax-1]                       # maximal measurable time in us
print "T = ", T, " us"
sigma_tau = 1./N * np.sqrt(np.sum(counts[cMin-1:cMax] * t[cMin-1:cMax]**2))
print sigma_tau

x_0 = 0.2                           # start value of the iteration 
global x_i
x_i =  x_0  
x = np.linspace(10e-9,1,100)        # sample points for plotting f

fig = plt.figure("Mittlere Lebensdauer von Myonen", figsize=(15,8))
axHist = fig.add_subplot(111)
axHist.set_xlabel(r'$x = \frac{\tau}{T}$')
axHist.set_ylabel(r'$f(x)$')
axHist.set_xlim([0, 1])
axHist.set_title(r'Mittlere Lebensdauer von Myonen mit Exponentialfunktion')
axHist.plot(x, f(x, t, counts, T, N, cMin, cMax))
axHist.plot(x, np.zeros(len(x)))
matplotlib.rcParams.update({'font.size': 20})

plt.connect('button_press_event', onMouseClick)
plt.show()

