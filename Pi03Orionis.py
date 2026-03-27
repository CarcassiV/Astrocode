from astropy.io import fits
import oifits
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, j0, j1
import scipy.integrate as integrate
import scipy.stats as stats
import sympy as sp
import random
from oifitsTools import openFile, limbdarkenedThetaTest, uniformDiskFitBootstrap, uniformDiskFitErrorBarTest, uniformDiskFitErrorBarTestWithCenterVisibility


visibilitiesSquaredSigGem, visibilitiesSquaredErrSigGem, spatialFrequenciesSigGem, spatialFrequenciesClosurePhasesSigGem,closurePhasesSigGem,closurePhasesErrSigGem = openFile('CHARADataSigGem/2011Dec07.17ms.sigGem.oifits')

visibilitiesSquaredCHARA, visibilitiesSquaredErrCHARA, spatialFrequenciesCHARA, spatialFrequenciesClosurePhasesCHARA, closurePhasesCHARA, closurePhasesErrCHARA = openFile('CHARADataPiOri/MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')

visibilitiesSquaredVLTI, visibilitiesSquaredErrVLTI, spatialFrequenciesVLTI, spatialFrequenciesClosurePhasesVLTI, closurePhasesVLTI, closurePhasesErrVLTI = openFile('VLTIDataPiOri/PIONI.2019-11-28T04-54-51.340_oidataCalibrated.fits')
print(visibilitiesSquaredVLTI)

print(len(visibilitiesSquaredVLTI))

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-28T05-13-08.623_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-29T05-18-37.271_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-29T05-34-15.601_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T08p16p19.858_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T08p44p52.527_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T09p12p58.922_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T09p45p15.562_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)


#uniformDiskTheta, uniformDiskError, chiSquareValues = uniformDiskFitErrorBarTest(visibilitiesSquaredSigGem, visibilitiesSquaredErrSigGem, spatialFrequenciesSigGem, 1)
#print("Uniform Disk Theta:", uniformDiskTheta, " Error:", uniformDiskError)

#uniformDiskTheta, visibility = uniformDiskFitErrorBarTestWithCenterVisibility(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 3)
#print("Uniform Disk Theta:", uniformDiskTheta, "Visibility", visibility)

#limbdarkenedTheta, alpha = limbdarkenedThetaTest(visibilitiesSquaredSigGem, visibilitiesSquaredErrSigGem, spatialFrequenciesSigGem, 1)
#print("Limb darkened theta", limbdarkenedTheta/((1/1000)*(1/60)*(1/60)*(np.pi/180)))
#print("Limb darkened coefficient", alpha)

uniformDiskTheta = 1.412*((1/1000)*(1/60)*(1/60)*(np.pi/180)) # VLTI results
limbdarkenedTheta = 1.41*((1/1000)*(1/60)*(1/60)*(np.pi/180)) 
alpha = 0

# Sigma Gem Results
#uniformDiskTheta = 2.256*((1/1000)*(1/60)*(1/60)*(np.pi/180))
#limbdarkenedTheta = 2.26*((1/1000)*(1/60)*(1/60)*(np.pi/180)) 
#alpha = 0.0

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True, figsize=(11,8))
print("about to plot")

# CHARA Data
ax[0].errorbar(spatialFrequenciesCHARA, visibilitiesSquaredCHARA, yerr=visibilitiesSquaredErrCHARA, fmt = '.', color='grey')
ax[1].set_xlabel('Spatial Frequency (Mλ)(baseline/wavelength)')
ax[0].set_ylabel('Visibilities Squared')
ax[0].set_yscale('log', base=10)
ax[1].errorbar(spatialFrequenciesClosurePhasesCHARA, closurePhasesCHARA, yerr=closurePhasesErrCHARA, fmt = '.')
ax[1].set_ylabel('Closure Phases (degrees)')

# VLTI Data
ax[0].errorbar(spatialFrequenciesVLTI, visibilitiesSquaredVLTI, yerr=visibilitiesSquaredErrVLTI, fmt = '*')
ax[0].set_yscale('log', base=10)
ax[1].errorbar(spatialFrequenciesClosurePhasesVLTI, closurePhasesVLTI, yerr=closurePhasesErrVLTI, fmt = '*')

#Plot models
ax[0].plot( #uniform disk model
    x, 
    ((2*j1(np.pi*uniformDiskTheta*x*1e6))/(np.pi*uniformDiskTheta*x*1e6))**2, 
    label='Uniform Disk Model',
    color="red",
    linewidth=2
) 

#uniform disk model
limbDarkenedValues = []
i = 0
while i < np.size(x):
    function = lambda r : ((1-r**2)**(alpha/2))*j0(np.pi*limbdarkenedTheta*r*x[i]*1e6)*r
    integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
    limbDarkenedValues.append(((alpha+2)*np.abs(integral))**2)
    i+=1

ax[0].plot(
    x, 
    limbDarkenedValues, 
    label='Limb Darkened Model',
    color="green",
    linewidth=2
)
ax[0].legend()

# Sigma Gem
#ax[0].errorbar(spatialFrequenciesSigGem, visibilitiesSquaredSigGem, yerr=visibilitiesSquaredErrSigGem, fmt = '*')
#ax[1].errorbar(spatialFrequenciesClosurePhasesSigGem, closurePhasesSigGem, yerr=closurePhasesErrSigGem, fmt = '*')

# Define font sizes
SIZE_DEFAULT = 17

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = SIZE_DEFAULT
plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.05)
plt.legend()
plt.show()


#uniform disk chisquared values plot
"""chiSquareValuesList = list(chiSquareValues.values())
i = 0
while i < np.size(chiSquareValuesList):
    if(chiSquareValuesList[i] > 500):
        chiSquareValuesList[i] = np.nan
    i += 1

figThree, axThree = plt.subplots()

axThree.plot(list(chiSquareValues.keys()), chiSquareValuesList)
axThree.set_xlabel("Angular Diameter in Milliarcseconds")
axThree.set_ylabel("Chi Squared Value")
axThree.set_title("Uniform Disk Model Chi Squared Values")"""

#Different alpha value plots
"""mew = np.arange(0, 1, .0001)
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
"""