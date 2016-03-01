# -*- coding: utf-8 -*-
"""
This programm determines the mean lifetime of myons
with the maximum-likelihood-method. The distribution
is assumed as poisson-distribution.

Input: N measured samples {t_1,...,t_N}
Output: Estimation value for the mean lifetime
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def readFile(fn, startRow):
    """Read data from file fn. Returns an array with the data"""
    data = np.loadtxt(fn, skiprows = startRow)
    return data

def getN0(tau, N, t, cWdth, cMin, cMax):
    """Calculates the norming constant N0"""
    N0 = N/(np.e**(-t[cMin-1]/tau)-np.e**(-(t[cMax-1]+cWdth)/tau))
    return N0

def getfi(tau, N0, cWdth, ti):
    """Calculate the expectation values f for the channels"""
    fi = N0 * cWdth * np.e**(-(ti + cWdth/2.)/tau) /tau
    return fi

def getLikelihood(tau, t, N0, cWdth, counts, N):
    funcVal = np.zeros(len(tau))
    f = np.zeros((len(tau), len(t)))
    for i in range(len(t)):
        f[:,i] = getfi(tau, N0, cWdth, t[i])
        funcVal +=    -2 * counts[i] * np.log(f[:,i])
    return funcVal

           
def main():
    fn = "mdat.txt"                     # file name for measuredata
    startRow = 85                       # start row
    data = readFile(fn, startRow)       # array with measuredata
    cNr = data[:,0]                     # channelnumbers
    cWdth = 1./24                       # channelwidth
    t = cNr * cWdth                     # measureable lifetimes
    cMin = 20                            # lower bond for channels
    cMax = 175                          # upper bond for channels
    counts =  data[:,1]                 # counts per channel
    countErr = np.sqrt(counts)          # error on the counts
    N = np.sum(counts[cMin-1:cMax])     # total number of counted events
    tau = np.linspace(t[0], t[cMax-1],1000)# possible mean lifetimes in us
    N0 = getN0(tau, N, t, cWdth, cMin, cMax)  # norming constant
    # effective likelihood-function
    L_eff = getLikelihood(tau, t[cMin-1:cMax],
                      N0, cWdth, 
                      counts[cMin-1:cMax], N)   
    tau_ew = tau[np.argmin(L_eff)]
    print "Mean Lifetime tau_ew = ",tau_ew, " us"
    L_min = np.min(L_eff)
    tau_sigma =  tau[np.abs(L_eff - L_min) < 1]
    print tau_sigma
    dTau = np.max(np.abs(tau_ew - tau_sigma))
    print "With standarddeviation sigma_tau = ", dTau, " us"

    
    
    fig = plt.figure("Mittlere Lebensdauer von Myonen", figsize=(15,8))
    axHist = fig.add_subplot(111)
    axHist.set_xlabel(r'$\tau\ /\mu s$')
    axHist.set_ylabel(r'$-2\ \ln{(L)}$')
    axHist.set_xlim([t[0], t[cMax-1]])
    axHist.set_title(r'Mittlere Lebendauer von Myonen mit'+
		     ' Poissonverteilung $\tau = $' +
		     str(round(tau_ew,3)) + r'$\pm$'+ str(round(dTau,3))
			+ r' $\mu s$')
    axHist.plot(tau, L_eff)
    axHist.plot(tau_ew, L_min, ls='', marker='o')
    matplotlib.rcParams.update({'font.size': 20})
    plt.show()

if __name__ == '__main__':
    main()    




