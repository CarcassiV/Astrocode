from astropy.io import fits
import oifits
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
import scipy.integrate as integrate
import random
from dataclasses import dataclass
from typing import Self

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
#   Conduct a literature review and compile a list of stellar paramaters for pi03Ori
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

@dataclass
class chiSquareTestResult:
    angularDiameter: np.float64
    alpha: np.float64
    chiSquare: np.float64

    def compareTo(self, other: Self) -> Self:
        if(other.chiSquare < self.chiSquare):
            return other
        else:
            return self

twoDVisibilities = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDVisibilities.append(np.ma.getdata(oifitsobj.vis2[i].vis2data))
    i += 1

twoDVisibilitiesError = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDVisibilitiesError.append(np.ma.getdata(oifitsobj.vis2[i].vis2err))
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
twoDSpatialFrequencyFiveNights = []
i = 0
while i < np.size(oifitsobj.vis2):
    twoDSpatialFrequency.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    if(i<5):
        twoDSpatialFrequencyFiveNights.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    i += 1 

visibilitiesSquared = flatten(twoDVisibilities)
visibilitiesSquaredErr = flatten(twoDVisibilitiesError)
spatialFrequencies = flatten(twoDSpatialFrequency)
#spatialFrequenciesFiveNights = flatten(twoDSpatialFrequencyFiveNights, 5, 8)
closurePhases = flatten(twoDClosurePhases)
closurePhasesErr = flatten(twoDClosurePhasesErrors)


def uniformDiskFitErrorBarTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    thetas = []
    #for each visibility, randomly sample a point on the error bar
    for i in range(0, numberOfTrials):
        sampleVisibilitiesSquared = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i]) #should be normal distribution 
            sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)
        #with these points, for each theta from 1.7milliarc second to 2.7 milliarc second, calculate the 
        #chi square value and find the optimal theta for that dataset
        chiSquareValues = {}

        thetaMilliArcSeconds = np.arange(1, 1.6, 0.01)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians):
            chiSquare = 0
            j = 0
            while j < np.size(sampleVisibilitiesSquared):
                observed = sampleVisibilitiesSquared[j]
                expected = ((2*jv(1, np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))/(np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))**2
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                j += 1
            print('Theta:', thetaRadians[i], ', Chi Squared Value:', chiSquare)
            chiSquareValues.update({thetaRadians[i]: chiSquare})
            i+=1
        thetas.append((min(chiSquareValues, key=chiSquareValues.get)))
        print(min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180))) 
    #repeat 500 times
    #take the average of all the theta
    print("done")
    i = 0
    sum = 0
    while i < np.size(thetas):
        sum += thetas[i]
        i += 1
    theta = sum/np.size(thetas)
    print(theta/((1/1000)*(1/60)*(1/60)*(np.pi/180)))
    return theta

#uniformDiskTheta = uniformDiskFitErrorBarTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 100)

def limbdarkenedThetaTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    chiSquareTestValues = []
    minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e30)
    #limb-darkened model angular diameter
    for a in range(0, numberOfTrials):
        sampleVisibilitiesSquared = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i]) #should be normal distribution 
            sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)
        sampleVisibilities = sampleVisibilitiesSquared
        
        thetaMilliArcSeconds = np.arange(1.2, 1.55, 0.01)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians): 
            j = 0
            minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e10)
            alpha = np.arange(0, .3, 0.01)
            while j < np.size(alpha):
                k = 0
                chiSquare = 0
                alphaValuesForOneTheta = {}
                while k < np.size(sampleVisibilities):
                    observed = sampleVisibilities[k]
                    function = lambda r : ((1-r**2)**(alpha[j]/2))*jv(0, np.pi*thetaRadians[i]*r*spatialFrequencies[k]*1e6)*r
                    integral, err = integrate.quad(function, 0, 1)
                    expected = ((alpha[j]+2)*integral)**2
                    chiSquareValue = ((observed-expected)**2)/expected
                    if not np.isnan(chiSquareValue):
                        chiSquare += chiSquareValue
                    k += 1
                print("Theta value:", thetaRadians[i], "alpha value:", alpha[j], "chisquare:", chiSquare)
                if (chiSquare < minChiSquareTestResultForATheta.chiSquare):
                    minChiSquareTestResultForATheta = chiSquareTestResult(thetaRadians[i], alpha[j], chiSquare)
                    #print(minChiSquareTestResultForATheta)
                    #print('switched value')
                j += 1
            print(minChiSquareTestResultForATheta)
            chiSquareTestValues.append(minChiSquareTestResultForATheta)
            i += 1
        a+=1
    print('done')
    i=0
    min = chiSquareTestResult(0,0,1e30)
    while i < np.size(chiSquareTestValues):
        if chiSquareTestValues[i].chiSquare < min.chiSquare:
            min = chiSquareTestValues[i]
            print('computing')
        i += 1
    print(min)
    theta = min.angularDiameter
    alpha = min.alpha
    return theta, alpha

limbdarkenedTheta, alpha = limbdarkenedThetaTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 1)
uniformDiskTheta = 1.48*((1/1000)*(1/60)*(1/60)*(np.pi/180))
#limbdarkenedTheta = 1.4*((1/1000)*(1/60)*(1/60)*(np.pi/180))
#alpha = 0.05

print(len(spatialFrequencies))
print(len(closurePhases))

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].plot(spatialFrequencies, visibilitiesSquared, '.')
ax[0].plot(x, ((2*jv(1, np.pi*uniformDiskTheta*x*1e6))/(np.pi*uniformDiskTheta*x*1e6))**2)
ax[0].errorbar(spatialFrequencies, visibilitiesSquared, yerr=visibilitiesSquaredErr, fmt = '.')
ax[0].set_ylabel('Visibilities Squared')
ax[0].set_yscale('log', base=10)

limbDarkenedValues = []
i = 0
while i < np.size(x):
    function = lambda r : ((1-r**2)**(alpha/2))*jv(0, np.pi*limbdarkenedTheta*r*x[i]*1e6)*r
    integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
    limbDarkenedValues.append((alpha+2)*np.abs(integral))
    i+=1

ax[0].plot(x, limbDarkenedValues)

"""ax[1].plot(spatialFrequencies, closurePhases, '.')
ax[1].errorbar(spatialFrequencies, closurePhases, yerr=closurePhasesErr, fmt = '.')
ax[1].set_xlabel('Spatial Frequency (Mλ)')
ax[1].set_ylabel('Closure Phases (degrees)')"""

mew = np.arange(0, 1, .0001)
r = np.arange(0, 1, 0.0001)
limbDarkeningCoefficient = [0, 0.2, 0.5, 1, 1.5, 3, 7]

"""figTwo, axTwo = plt.subplots(2, 1)

for i in range(0, np.size(limbDarkeningCoefficient)):
    axTwo[0].plot(mew, mew**limbDarkeningCoefficient[i])
    axTwo[1].plot(r, (1-r**2)**(limbDarkeningCoefficient[i]/2))
axTwo[0].set_xlabel('μ')
axTwo[1].set_xlabel('r')
axTwo[0].set_ylabel('Intensity, I(r)')
axTwo[1].set_ylabel('Intensity, I(r)')"""

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()