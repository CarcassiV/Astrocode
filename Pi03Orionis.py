from astropy.io import fits
import oifits
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
import scipy.integrate as integrate
import scipy.stats as stats
import sympy as sp
import random
from pi03OrionisParameters import uniformDiskFitErrorBarTest, limbdarkenedThetaTest

hdulist = fits.open('MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')
hdulist['OI_ARRAY'].header['OI_REVN'] = 1
oifitsobj = oifits.open(hdulist)

#To Do:
#   check, Fit uniform disk theta, should get around 1. something milliarc seconds
#        -with 100 trials using the error bar method, I received a theta of 1.4702000000000017 milliarc seconds
#   Fit limb darkened theta
#        -with 1 trial, angular diameter of 1.47 milliarc seconds and alpha of 0.01. One trial took 43 mins T-T. How can I make it more efficient????
#   check, Graph visibilities
#   Graph closure phases, why are there more closure phases than spatial frequencies?
#   Conduct a literature review and compile a list of stellar parameters for pi03Ori
#   Solve for my own parameters
#   Figure out how to use Candid to look for a binary, keep good data about the results

def flatten(ndarr): #assumes each row is of the same length
    rows = len(ndarr)
    cols = len(ndarr[0])
    flatarr = []
    for i in range(0,rows):
        for j in range(0,cols):
            flatarr.append(ndarr[i][j])
    return flatarr

twoDVisibilities = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDVisibilities.append(np.ma.getdata(oifitsobj.vis2[i].vis2data))
    i += 1

twoDVisibilitiesNotSquared = []
i = 0
while i < np.size(oifitsobj.vis):
    twoDVisibilitiesNotSquared.append(np.ma.getdata(oifitsobj.vis[i].visamp))
    i += 1

twoDVisibilitiesError = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDVisibilitiesError.append(np.ma.getdata(oifitsobj.vis2[i].vis2err))
    i += 1

twoDVisibilitiesNotSquaredErr = []
i = 0
while i < np.size(oifitsobj.vis):
    twoDVisibilitiesNotSquaredErr.append(np.ma.getdata(oifitsobj.vis[i].visamperr))
    i += 1

twoDClosurePhases = []
i = 0
while i < np.size(oifitsobj.t3):
    twoDClosurePhases.append(np.ma.getdata(oifitsobj.t3[i].t3phi))
    i += 1

twoDClosurePhasesErrors = []
i = 0
while i < np.size(oifitsobj.t3):
    twoDClosurePhasesErrors.append(np.ma.getdata(oifitsobj.t3[i].t3phierr))
    i += 1

twoDSpatialFrequency = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDSpatialFrequency.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    i += 1

visibilitiesSquared = flatten(twoDVisibilities)
print(len(visibilitiesSquared))
visibilitiesNotSquared = flatten(twoDVisibilitiesNotSquared)
visibilitiesSquaredErr = flatten(twoDVisibilitiesError)
visibilitiesNotSquaredErr = flatten(twoDVisibilitiesNotSquaredErr)
spatialFrequencies = flatten(twoDSpatialFrequency)
closurePhases = flatten(twoDClosurePhases)
closurePhasesErr = flatten(twoDClosurePhasesErrors)
print(len(closurePhases))

#uniformDiskTheta, uniformDiskError, visibilityAtCenter, visibilityAtCenterError = uniformDiskFitErrorBarTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 1000)
#print("Uniform Disk Theta:", uniformDiskTheta, " Error:", uniformDiskError, "Visibility at Center:", visibilityAtCenter, "Error:", visibilityAtCenterError)

limbdarkenedTheta, alpha = limbdarkenedThetaTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 1)
print("Limb darkened theta", limbdarkenedTheta/((1/1000)*(1/60)*(1/60)*(np.pi/180)))
print("Limb darkened coefficient", alpha)

uniformDiskTheta = 1.4845949999999937*((1/1000)*(1/60)*(1/60)*(np.pi/180))
limbdarkenedTheta = 1.4*((1/1000)*(1/60)*(1/60)*(np.pi/180))
alpha = 0.05

x = np.arange(10, 225, .2) #for the visibility squared curve

#Visibility Squared plot
fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].plot(spatialFrequencies, visibilitiesSquared, '.')
ax[0].errorbar(spatialFrequencies, visibilitiesSquared, yerr=visibilitiesSquaredErr, fmt = '.')
ax[1].set_xlabel('Spatial Frequency (Mλ)(baseline/wavelength)')
ax[0].set_ylabel('Visibilities Squared')
#ax[0].set_yscale('log', base=10)

#Closure Phases plot
"""ax[1].plot(spatialFrequencies, closurePhases, '.')
ax[1].errorbar(spatialFrequencies, closurePhases, yerr=closurePhasesErr, fmt = '.')
ax[1].set_ylabel('Closure Phases (degrees)')"""

#Different alpha value plots
mew = np.arange(0, 1, .0001)
r = np.arange(0, 1, 0.0001)
limbDarkeningCoefficient = [0, 0.2, 0.5, 1, 1.5, 3, 7]

figTwo, axTwo = plt.subplots(2, 1)

for i in range(0, np.size(limbDarkeningCoefficient)):
    axTwo[0].plot(mew, mew**limbDarkeningCoefficient[i])
    axTwo[1].plot(r, (1-r**2)**(limbDarkeningCoefficient[i]/2))
axTwo[0].set_xlabel('μ')
axTwo[1].set_xlabel('r')
axTwo[0].set_ylabel('Intensity, I(r)')
axTwo[1].set_ylabel('Intensity, I(r)')

#Plot models
ax[0].plot(x, ((2*jv(1, np.pi*uniformDiskTheta*x*1e6))/(np.pi*uniformDiskTheta*x*1e6))**2, label='Uniform Disk Model') #uniform disk model

limbDarkenedValues = []
i = 0
while i < np.size(x):
    function = lambda r : ((1-r**2)**(alpha/2))*jv(0, np.pi*limbdarkenedTheta*r*x[i]*1e6)*r
    integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
    limbDarkenedValues.append((alpha+2)*np.abs(integral))
    i+=1

ax[0].plot(x, limbDarkenedValues, label='Limb Darkened Model')
ax[0].legend()

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()