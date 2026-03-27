import numpy as np
from astropy.io import fits
import oifits
import matplotlib.pyplot as plt
from scipy.special import jv
import scipy.integrate as integrate
import random
from dataclasses import dataclass
from typing import Self

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

def flatten(ndarr, rows, cols):
    flatarr = []
    for i in range(0,rows):
        for j in range(0,cols):
            flatarr.append(ndarr[i][j])
    return flatarr

def uniformDiskFitErrorBar(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
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

        thetaMilliArcSeconds = np.arange(1.7, 2.7, 0.0001)
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
    print(theta)
    return theta

def uniformDiskFitBootstrap(visibilitiesSquared, spatialFrequencies, numberOfTrials):
    #makes a dictionary where each visibility will correspond with a certain spatial frequency
    visDict = {}
    for i in range(0, np.size(spatialFrequencies)):
        if not np.isnan(visibilitiesSquared[i]):
            visDict.update({visibilitiesSquared[i]:spatialFrequencies[i]})

    thetas = []
    #for each visibility, randomly sample a point on the error bar
    for i in range(0, numberOfTrials):
        sampleVisibilities = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySample = visibilitiesSquared[random.randint(0, np.size(visibilitiesSquared)-1)]
            if not np.isnan(randomVisibilitySample):
                sampleVisibilities.append(randomVisibilitySample)
        #with these points, for each theta from 1.7milliarc second to 2.7 milliarc second, calculate the 
        #chi square value and find the optimal theta for that dataset
        chiSquareValues = {}

        thetaMilliArcSeconds = np.arange(1.7, 2.7, 0.0001)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians):
            chiSquare = 0
            j = 0
            while j < np.size(sampleVisibilities):
                observed = sampleVisibilities[j]
                expected = ((2*jv(1, np.pi*thetaRadians[i]*visDict.get(sampleVisibilities[j])*1e6))/(np.pi*thetaRadians[i]*visDict.get(sampleVisibilities[j])*1e6))**2
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                j += 1
            print('Theta:', thetaRadians[i], ', Chi Squared Value:', chiSquare)
            chiSquareValues.update({thetaRadians[i]: chiSquare})
            i+=1
        thetas.append((min(chiSquareValues, key=chiSquareValues.get)))
        print(min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180))) 
        #repeat number of trials amount of times
    #take the average of all the theta
    print("done")
    i = 0
    sum = 0
    while i < np.size(thetas):
        sum += thetas[i]
        i += 1
    theta = sum/np.size(thetas)
    print(theta)
    return theta

oifitsSigGem = oifits.open('CHARADataSigGem/2011Dec07.17ms.sigGem.oifits')

twoDVisibilitiesSigGem = []
i = 0
while i < np.size(oifitsSigGem.vis2):
    twoDVisibilitiesSigGem.append(np.ma.getdata(oifitsSigGem.vis2[i].vis2data))
    i += 1

twoDVisibilitiesErrorSigGem = []
i = 0
while i < np.size(oifitsSigGem.vis2):
    twoDVisibilitiesErrorSigGem.append(np.ma.getdata(oifitsSigGem.vis2[i].vis2err))
    i += 1

twoDClosurePhasesSigGem = []
i = 0
while i < np.size(oifitsSigGem.t3):
    twoDClosurePhasesSigGem.append(np.ma.getdata(oifitsSigGem.t3[i].t3phi))
    i += 1

twoDClosurePhasesErrorsSigGem = []
i = 0
while i < np.size(oifitsSigGem.t3):
    twoDClosurePhasesErrorsSigGem.append(np.ma.getdata(oifitsSigGem.t3[i].t3phierr))
    i += 1

twoDSpatialFrequencySigGem = []
i = 0
while i < np.size(oifitsSigGem.vis2):
    twoDSpatialFrequencySigGem.append(np.sqrt((oifitsSigGem.vis2[i].ucoord)**2 + (oifitsSigGem.vis2[i].vcoord)**2)/oifitsSigGem.vis2[i].wavelength.eff_wave/1e6)
    i += 1 

twoDSpatialFrequencyClosurePhasesSigGem = []
i = 0
while i < np.size(oifitsSigGem.t3):
    twoDSpatialFrequencyClosurePhasesSigGem.append(np.sqrt((oifitsSigGem.t3[i].u1coord)**2 + (oifitsSigGem.t3[i].v1coord)**2)/oifitsSigGem.t3[i].wavelength.eff_wave/1e6)
    i += 1

visibilitiesSquared = flatten(twoDVisibilitiesSigGem, 8, 8)
visibilitiesSquaredErr = flatten(twoDVisibilitiesErrorSigGem, 8, 8)
spatialFrequencies = flatten(twoDSpatialFrequencySigGem, 8, 8)
closurePhases = flatten(twoDClosurePhasesSigGem, 5, 8)
closurePhasesErr = flatten(twoDClosurePhasesErrorsSigGem, 5, 8)
spatialFrequenciesClosurePhases = flatten(twoDSpatialFrequencyClosurePhasesSigGem, 5, 8)


hdulist = fits.open('CHARADataPiOri/MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')
hdulist['OI_ARRAY'].header['OI_REVN'] = 1
oifitsobj = oifits.open(hdulist)

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

twoDSpatialFrequencyClosurePhases = []
i = 0
while i < np.size(oifitsobj.t3):
    twoDSpatialFrequencyClosurePhases.append(np.sqrt((oifitsobj.t3[i].u1coord)**2 + (oifitsobj.t3[i].v1coord)**2)/oifitsobj.t3[i].wavelength.eff_wave/1e6)
    i += 1

visibilitiesSquaredPiOri = flatten(twoDVisibilities, 8, 8)
visibilitiesSquaredErrPiOri = flattenvisibilitiesNotSquaredErr = flatten(twoDVisibilitiesNotSquaredErr, 8, 8)
spatialFrequenciesPiOri = flatten(twoDSpatialFrequency, 8, 8)
closurePhasesPiOri = flatten(twoDClosurePhases, 8, 8)
closurePhasesErrPiOri = flatten(twoDClosurePhasesErrors, 8, 8)
spatialFrequenciesClosurePhasesPiOri = flatten(twoDSpatialFrequencyClosurePhases, 8, 8)

#uniformDiskTheta = uniformDiskFitErrorBar(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, 500)

chiSquareTestValues = []
minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e30)
#limb-darkened model angular diameter
for a in range(0, 1):
    sampleVisibilitiesSquared = []
    for i in range(0, np.size(visibilitiesSquared)):
        randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i]) #should be normal distribution 
        sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)
    sampleVisibilities = sampleVisibilitiesSquared

    sampleVisibilities = visibilitiesSquared
    
    thetaMilliArcSeconds = np.arange(2.2, 2.45, 0.01)
    thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
    i = 0
    while i < np.size(thetaRadians): 
        j = 0
        minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e30)
        alpha = np.arange(0, .4, 0.01)
        while j < np.size(alpha):
            k = 0
            chiSquare = 0
            alphaValuesForOneTheta = {}
            while k < np.size(sampleVisibilities):
                observed = sampleVisibilities[k]
                #print(alpha[j])
                function = lambda r : ((1-r**2)**(alpha[j]/2))*jv(0, np.pi*thetaRadians[i]*r*spatialFrequencies[k]*1e6)*r
                integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
                expected = ((alpha[j]+2)*integral)**2
                #if(alpha[j] == 0 and k == 1):
                    #print(sampleVisibilities[k], spatialFrequencies[k], thetaRadians[i])
                    #print(expected)
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                k += 1
            #print(chiSquare)
            print("Theta value:", thetaRadians[i], "alpha value:", alpha[j], "chisquare:", chiSquare)
            if (chiSquare < minChiSquareTestResultForATheta.chiSquare): #the less than comparison isnt currently working properly it seems
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

thetaSigGem = 2.335*((1/1000)*(1/60)*(1/60)*(np.pi/180))
limbDarkenedThetaSigGem = 2.4*((1/1000)*(1/60)*(1/60)*(np.pi/180))
alphaSigGem = 0.25

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
#ax[0].plot(spatialFrequenciesSigGem, visibilitiesSquaredSigGem, '.')
#ax[0].plot(x, ((2*jv(1, np.pi*thetaSigGem*x*1e6))/(np.pi*thetaSigGem*x*1e6))**2)
ax[0].errorbar(spatialFrequenciesSigGem, visibilitiesSquaredSigGem, yerr=visibilitiesSquaredErrSigGem, fmt = '*')
ax[0].set_ylabel('Visibilities Squared')
#ax[0].set_yscale('log', base=10)

#ax[0].plot(spatialFrequenciesPiOri, visibilitiesSquaredPiOri, '.')
ax[0].errorbar(spatialFrequenciesPiOri, visibilitiesSquaredPiOri, yerr=visibilitiesSquaredErrPiOri, fmt = '.')

limbDarkenedValues = []
i = 0
while i < np.size(x):
    function = lambda r : ((1-r**2)**(alphaSigGem/2))*jv(0, np.pi*limbDarkenedThetaSigGem*r*x[i]*1e6)*r
    integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
    limbDarkenedValues.append((alphaSigGem+2)*np.abs(integral))
    i+=1

#ax[0].plot(x, limbDarkenedValues)

ax[1].errorbar(spatialFrequenciesClosurePhasesSigGem, closurePhasesSigGem, yerr=closurePhasesErrSigGem, fmt = '*')
ax[1].set_xlabel('Spatial Frequency (Mλ)')
ax[1].set_ylabel('Closure Phases (degrees)')

ax[1].errorbar(spatialFrequenciesClosurePhasesPiOri, closurePhasesPiOri, yerr=closurePhasesErrPiOri, fmt = '.')

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

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()