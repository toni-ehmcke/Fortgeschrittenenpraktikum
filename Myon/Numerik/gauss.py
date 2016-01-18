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

def getLikelihood_gauss(tau, t, N0, cWdth, counts, N):
    funcVal = np.zeros(len(tau))
    f = np.zeros((len(tau), len(t)))
    for i in range(len(t)):
        f[:,i] = getfi(tau, N0, cWdth, t[i])
        funcVal +=    (counts[i]-f[:,i])**2 /counts[i]
    return funcVal

           
def main():
    fn = "mdat.txt"                     # file name for measuredata
    startRow = 85                       # start row
    data = readFile(fn, startRow)       # array with measuredata
    
    #--------------------values without combining channels--------------------    
    cNr = data[:,0]                     # channelnumbers
    cWdth = 2./48
    t = cNr * cWdth                         # measureable lifetimes
    cMin = 20                               # lower bond for channels
    cMax = 175                              # upper bond for channels
    counts =  data[:,1]                     # counts per channel
    N = np.sum(counts[cMin-1:               # total number of counted events
                      cMax])
    tau = np.linspace(t[0],                 # possible mean lifetimes in us
                      t[cMax-1],1000) 
    # norming constant
    N0 = getN0(tau, N, t, cWdth,
                     cMin, cMax)
    # effective likelihood-function
    L_eff = getLikelihood_gauss(tau, t[cMin-1:cMax],
                      N0, cWdth, 
                      counts[cMin-1:cMax], N) 
    
    tau_ew = tau[np.argmin(L_eff)]
    print "Mean Lifetime tau_ew without combining channels = ",tau_ew, " us"
    L_min = np.min(L_eff)
    tau_sigma =  tau[np.abs(L_eff - L_min) < 1]
    dTau = np.max(np.abs(tau_ew - tau_sigma))
    print "With standarddeviation sigma_tau = ", dTau, " us"
    
    #---------------------special values for gaussian distribution------------
    counts_gauss = np.zeros(len(counts) # combine neigboured channel-counts
                            /2) 
    for i in range(len(counts_gauss)):
        counts_gauss[i] = counts[2*i] + counts[2*i+1]
    cWdth_gauss = 4./48                                 # channelwidth
    t_gauss = cWdth_gauss * cNr[0:len(cNr)/2] - 1./48   # measureable lifetimes
    cMin_gauss = np.ceil(cMin / 2.)                     # lower bond for channels
    cMax_gauss = np.ceil(cMax / 2.)                     # upper bond for channels
    N_gauss = np.sum(counts[cMin_gauss-1:               # total number of counted events
                      cMax_gauss])
    tau_gauss = np.linspace(t_gauss[0],                 # possible mean lifetimes in us
                      t_gauss[cMax_gauss-1],1000) 
    # norming constant
    N0_gauss = getN0(tau_gauss, N_gauss, t_gauss, cWdth_gauss,
                     cMin_gauss, cMax_gauss)
    # effective likelihood-function
    L_eff_gauss = getLikelihood_gauss(tau_gauss, t_gauss[cMin_gauss-1:cMax_gauss],
                      N0_gauss, cWdth_gauss, 
                      counts_gauss[cMin_gauss-1:cMax_gauss], N_gauss)  
   
    tau_ew_gauss = tau_gauss[np.argmin(L_eff_gauss)]
    print "Mean Lifetime tau_ew  with combined channels= ",tau_ew_gauss, " us"
    L_min_gauss = np.min(L_eff_gauss)
    tau_sigma_gauss =  tau_gauss[np.abs(L_eff_gauss - L_min_gauss) < 1]
    dTau_gauss = np.max(np.abs(tau_ew_gauss - tau_sigma_gauss))
    print "With standarddeviation sigma_tau = ", dTau_gauss, " us"

    
    
    fig = plt.figure("Mittlere Lebensdauer von Myonen", figsize=(15,8))
    axHist = fig.add_subplot(111)
    axHist.set_xlabel(r'$\tau\ /\mu s$')
    axHist.set_ylabel(r'$-2 \ln{L}$')
    axHist.set_xlim([t_gauss[0], t_gauss[cMax_gauss-1]])
    axHist.set_title(r'Mittlere Lebensdauer von Myonen mit Gaussverteilung')
    axHist.plot(tau_gauss, L_eff_gauss)
    axHist.plot(tau, L_eff)
    axHist.text(3, 200000, r'Mittlere Lebensdauer mit kombinierten Kanaelen $\tau_k = $'
                            + str(round(tau_ew_gauss,2)) + r'$\pm$'+ str(round(dTau_gauss,2)))
    axHist.text(3, 220000, r'Mittlere Lebensdauer ohne kombinierte Kanaele $\tau = $' 
                        + str(round(tau_ew,2)) + r'$\pm$'+ str(round(dTau,2)))
    plt.show()

if __name__ == '__main__':
    main()    




