# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 09:32:26 2015

@author: toni
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm

N = 5                                                  # Dimension
BIN = 1                                                # bin width of filter
inp = np.zeros(shape=(N,N))                            # input matrix
inp[0,2] = 1
inp[2,0] = 3
inp[3,3] = 7

filt = np.zeros(2*BIN+1)
filt[0] = -.1
filt[1] = .25
filt[2] = -.1

WProj = np.zeros(shape=(9,5))         # weightingmatrix for projection
WProj[0,0] = .14
WProj[1,0] = .88
WProj[2,0] = .39
WProj[1,1] = .01
WProj[2,1] = .61
WProj[3,1] = .75
WProj[4,1] = .04
WProj[3,2] = .25
WProj[4,2] = .91
WProj[5,2] = .25
WProj[4,3] = .04
WProj[5,3] = .75
WProj[6,3] = .61
WProj[7,3] = .01
WProj[6,4] = .39
WProj[7,4] = .88
WProj[8,4] = .14

WInvProj = np.zeros(shape=(5,9))      # weightingmatrix for inverse projection
WInvProj[0,0] = .27
WInvProj[0,1] = 1.
WInvProj[0,2] = .73
WInvProj[1,1] = .01
WInvProj[1,2] = .93
WInvProj[1,3] = .99
WInvProj[1,4] = .07
WInvProj[2,3] = .5
WInvProj[2,4] = 1.
WInvProj[2,5] = .5
WInvProj[3,4] = .07
WInvProj[3,5] = .99
WInvProj[3,6] = .93
WInvProj[3,7] = .01
WInvProj[4,6] = .73
WInvProj[4,7] = 1.
WInvProj[4,8] = .27

#--------------------------- Fill the projection matrix------------------------
Proj45 = np.mat(np.zeros(shape=(1,2*N-1)))
Proj135 = np.mat(np.zeros(shape=(1,2*N-1)))

for i in np.arange(-N+1,N):
    Proj45[0,i+N-1] = np.sum(np.diag(inp, k = i)) 
    Proj135[0,-i-N] = np.sum(np.diag(inp[::-1], k=i))
    # mirror matrix with respect to the center row: secdiag --> maindiag
    # then fill the Proj135 array from last to first position

Proj0 = np.zeros(N)
Proj90 = np.zeros(N)

for i in np.arange(N):
    Proj0[i] = np.sum(inp[:,i])
    Proj90[-i-1] = np.sum(inp[i,:])

Proj45Wei = Proj45 * WProj                      # weighted projection 45 deg
Proj135Wei = Proj135 * WProj                    # weighted projection 135 deg

Proj = np.zeros(shape=(4,N+2))                      # Projection (unfiltered)
Proj[0,1:N+1] = Proj0
Proj[1,1:N+1] = Proj45Wei
Proj[2,1:N+1] = Proj90
Proj[3,1:N+1] = Proj135Wei

#---------------------------Fill filtered projection---------------------------
ProjFilt = np.zeros(shape=(4,N))                    # Projection (filtered)
for i in range(4):
     for j in range(N):
         ProjFilt[i,j] = np.dot(Proj[i,j:j+3],filt) 

#-------------------------Backprojection---------------------------------------
# unfiltered summands of the backprojection
Backproj0 = np.zeros(shape=(N,N))
Backproj45 = np.zeros(shape=(N,N))
Backproj90 = np.zeros(shape=(N,N))
Backproj135 = np.zeros(shape=(N,N))

# filtered summands of the backprojection
Backproj0Filt = np.zeros(shape=(N,N))
Backproj45Filt = np.zeros(shape=(N,N))
Backproj90Filt = np.zeros(shape=(N,N))
Backproj135Filt = np.zeros(shape=(N,N))
  
for i in range(N):
    Backproj0[i,:] = Proj0  
    Backproj90[::-1][:,i] = Proj90
    Backproj0Filt[i,:] = ProjFilt[0,:]
    Backproj90Filt[::-1][:,i] = ProjFilt[2,:]
    

Backproj45Wei = Proj45Wei*WInvProj              # weighted backprojection
Backproj135Wei = Proj135Wei*WInvProj
Backproj45WeiFilt = np.mat(ProjFilt[1,:]) * WInvProj # weighted backproj. filtered
Backproj135WeiFilt = np.mat(ProjFilt[3,:]) * WInvProj


for i in np.arange(-N+1,N):
    # iterate through the elements x of the weighted matrices of the backproj. and create 
    # matrices with x on the i main (45 deg) / secondary diagonal(135 deg)
    Backproj45 += np.diag(Backproj45Wei[0,i+N-1] * np.ones(N - np.abs(i)),k=i)
    Backproj135[::-1] += np.diag(Backproj135Wei[0,i+N-1] * np.ones(N - np.abs(i)),k=-i)
    Backproj45Filt += np.diag(Backproj45WeiFilt[0,i+N-1] * np.ones(N - np.abs(i)),k=i)
    Backproj135Filt[::-1] += np.diag(Backproj135WeiFilt[0,i+N-1] * np.ones(N - np.abs(i)),k=-i)
        
outUnfilt = Backproj0 + 0.5 * Backproj45 + Backproj90 \
                + 0.5 * Backproj135 # unfiltered output
outFilt = Backproj0Filt + 0.5 * Backproj45Filt \
                + Backproj90Filt + 0.5 * Backproj135Filt # unfiltered output

print outFilt

fig = plt.figure("Radon-Transformation", figsize=(15,9))

axInp = fig.add_subplot(221)
axInp.set_xlabel("location x")
axInp.set_ylabel("location y")
axInp.set_title("Ausgangsverteilung")

axOutUnfilt = fig.add_subplot(223)
axOutUnfilt.set_xlabel("location x")
axOutUnfilt.set_ylabel("location y")
axOutUnfilt.set_title("Ungefilterte Rueckprojektion")

axOutFilt = fig.add_subplot(224)
axOutFilt.set_xlabel("location x")
axOutFilt.set_ylabel("location y")
axOutFilt.set_title("Gefilterte Rueckprojektion")

print inp
plotInp = axInp.imshow(inp,interpolation='None', cmap = cm.Greys,
                       extent = [0,N,0,N])
fig.colorbar(plotInp, ax=axInp)

plotOutUnfilt = axOutUnfilt.imshow(outUnfilt,interpolation='None', cmap = cm.Greys,
                       extent = [0,N,0,N])
fig.colorbar(plotOutUnfilt, ax=axOutUnfilt)

plotOutFilt = axOutFilt.imshow(outFilt,interpolation='None', cmap = cm.Greys,
                       extent = [0,N,0,N])
fig.colorbar(plotOutFilt, ax=axOutFilt)

plt.show()