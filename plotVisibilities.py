import numpy as np
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
    alphaValues: np.float64
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
            #print('Theta:', thetaRadians[i], ', Chi Squared Value:', chiSquare)
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
            #print('Theta:', thetaRadians[i], ', Chi Squared Value:', chiSquare)
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

oifitsobj = oifits.open('2011Dec07.17ms.sigGem.oifits')

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

visibilitiesSquared = flatten(twoDVisibilities, 8, 8)
visibilitiesSquaredErr = flatten(twoDVisibilitiesError, 8, 8)
spatialFrequencies = flatten(twoDSpatialFrequency, 8, 8)
spatialFrequenciesFiveNights = flatten(twoDSpatialFrequencyFiveNights, 5, 8)
closurePhases = flatten(twoDClosurePhases, 5, 8)
closurePhasesErr = flatten(twoDClosurePhasesErrors, 5, 8)

chiSquareTestValues = []
minChiSquareTestResultForATheta = chiSquareTestResult(0,0,1e30)
#limb-darkened model angular diameter
for a in range(0, 1):
    sampleVisibilitiesSquared = []
    for i in range(0, np.size(visibilitiesSquared)):
        randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i]) #should be normal distribution 
        sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)
    sampleVisibilities = np.sqrt(sampleVisibilitiesSquared)
    
    thetaMilliArcSeconds = np.arange(1.7, 2.7, 0.01)
    thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
    i = 0
    while i < np.size(thetaRadians): 
        j = 0
        alpha = np.arange(0, 1, 0.01)
        while j < np.size(alpha):
            k = 0
            chiSquare = 0
            alphaValuesForOneTheta = {}
            while k < np.size(sampleVisibilities):
                observed = np.sqrt(sampleVisibilities[k])
                #print(observed)
                function = lambda r : ((1-r**2)**(alpha[j]/2))*jv(0, np.pi*thetaRadians[i]*r*spatialFrequencies[k])*r
                integral, err = integrate.quad(function, 0, 1)
                expected = (alpha[j]+2)*np.abs(integral)
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                k += 1
            #print("Theta value:", thetaRadians[i], "alpha value:", alpha[j], "chisquare:", chiSquare)
            if j == 0 or chiSquare < minChiSquareTestResultForATheta.chiSquare: #the less than comparison isnt currently working properly it seems
                minChiSquareTestResultForATheta = chiSquareTestResult(thetaRadians[i], alpha[j], chiSquare)
            j += 1
        print(minChiSquareTestResultForATheta)
        chiSquareTestValues.append(minChiSquareTestResultForATheta)
        i += 1
    a+=1
#while i < np.size(chiSquareTestValues):


theta = 2.335*((1/1000)*(1/60)*(1/60)*(np.pi/180))

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].plot(spatialFrequencies, visibilitiesSquared, '.')
ax[0].plot(x, ((2*jv(1, np.pi*theta*x*1e6))/(np.pi*theta*x*1e6))**2)
ax[0].errorbar(spatialFrequencies, visibilitiesSquared, yerr=visibilitiesSquaredErr, fmt = '.')
ax[0].set_ylabel('Visibilities Squared')
ax[0].set_yscale('log', base=10)

ax[1].plot(spatialFrequenciesFiveNights, closurePhases, '.')
ax[1].errorbar(spatialFrequenciesFiveNights, closurePhases, yerr=closurePhasesErr, fmt = '.')
ax[1].set_xlabel('Spatial Frequency (Mλ)')
ax[1].set_ylabel('Closure Phases (degrees)')

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