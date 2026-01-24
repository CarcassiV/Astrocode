import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, j0, j1
import random
import scipy.integrate as integrate
from dataclasses import dataclass
from typing import Self
import sympy as sp

def convertMilliArcSecToRadian(num):
    return num*((1/1000)*(1/60)*(1/60)*(np.pi/180))

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

#will take about 3hr minutes to run 1000 trials
#uniform disk model angular diameter
def uniformDiskFitErrorBarTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    thetas = []
    visibilitiesAtCenter = []
    #for each visibility, randomly sample a point on the error bar
    for n in range(0, numberOfTrials):
        sampleVisibilities = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i])
            sampleVisibilities.append(randomVisibilitySquaredSample)
        chiSquareValues = {}

        bottomThetaRange = 1.4
        topThetaRange = 1.55

        thetaMilliArcSeconds = np.arange(bottomThetaRange, topThetaRange, 0.001)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians):
            chiSquare = 0
            j = 0
            while j < np.size(sampleVisibilities):
                observed = sampleVisibilities[j]
                expected = ((2*j1(np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))/(np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))**2
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                j += 1
            print('Theta:', thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), ', Chi Squared Value:', chiSquare)
            chiSquareValues.update({thetaRadians[i]: chiSquare})
            i += 1
        thetas.append(min(chiSquareValues, key=chiSquareValues.get))
        print("Trial Number:", n, "Theta:", min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180))) 

        """x = sp.symbols('x')
        equation = ((2*sp.besselj(np.pi*min(chiSquareValues, key=chiSquareValues.get)*x*1e6, 1))/(np.pi*min(chiSquareValues, key=chiSquareValues.get)*x*1e6))**2
        visibilitiesAtCenter.append(sp.limit(equation, x, 0))"""
        
        #Test if this ever hits the extreme values of theta, if so expand the range. Try to make range as small as possible to save time
        if(min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180)) == bottomThetaRange or min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180)) == topThetaRange):
            print("Expand theta range, stopping loop")
            break
        n += 1
    #take the average of all the theta
    print("done")
    i = 0
    sum = 0
    while i < np.size(thetas):
        sum += thetas[i]
        i += 1
    theta = sum/np.size(thetas)
    print(theta/((1/1000)*(1/60)*(1/60)*(np.pi/180)))

    #The error is calculated from the standard deviation of the trials
    thetaError = np.std(thetas)

    """i = 0
    sum = 0
    while i < np.size(visibilitiesAtCenter):
        sum += visibilitiesAtCenter[i]
        i += 1
    visibilityAtCenter = sum/np.size(visibilitiesAtCenter)
    print(visibilityAtCenter)

    visibilityAtCenterError = np.std(visibilitiesAtCenter)
    print(visibilityAtCenterError)"""

    visibilityAtCenter = 0
    visibilityAtCenterError = 0

    return theta, thetaError, chiSquareValues, visibilityAtCenter, visibilityAtCenterError


#limb-darkened model angular diameter
def limbdarkenedThetaTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    chiSquareTestValues = []
    minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e10)
    
    for n in range(0, numberOfTrials):
        sampleVisibilitiesSquared = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i])
            sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)

        bottomThetaRange = 1.3
        topThetaRange = 2.1
        
        thetaMilliArcSeconds = np.arange(bottomThetaRange, topThetaRange, 0.1)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians): #are there libraries that already run chi square tests? More efficient than mine perhaps?
            j = 0
            minChiSquareTestResultForATheta = chiSquareTestResult(0,.5,1e10)

            bottomAlphaRange = 0.0
            topAlphaRange = 0.22

            alpha = np.arange(bottomAlphaRange, topAlphaRange, 0.01)
            while j < np.size(alpha):
                k = 0
                chiSquare = 0
                alphaValuesForOneTheta = {}
                while k < np.size(sampleVisibilitiesSquared):
                    observed = sampleVisibilitiesSquared[k]
                    function = lambda r : ((1-r**2)**(alpha[j]/2))*j0(0, np.pi*thetaRadians[i]*r*spatialFrequencies[k]*1e6)*r
                    integral, err = integrate.quad(function, 0, 1) #is quad an efficient function? Does python have more efficient integration?
                    expected = ((alpha[j]+2)*integral)**2
                    chiSquareValue = ((observed-expected)**2)/expected #if I do the whole thing with visibility not square more efficient?
                    if not np.isnan(chiSquareValue):
                        chiSquare += chiSquareValue
                    k += 1
                print("Theta value:", thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), "alpha value:", alpha[j], "chisquare:", chiSquare)
                if (chiSquare < minChiSquareTestResultForATheta.chiSquare): #I compare every theta and alpha pair, is it more efficient to store them all and sort/compare at end to find min?
                    minChiSquareTestResultForATheta = chiSquareTestResult(thetaRadians[i], alpha[j], chiSquare)
                j += 1
            print(minChiSquareTestResultForATheta)
            chiSquareTestValues.append(minChiSquareTestResultForATheta)
            print(chiSquareTestValues)
            """if(minChiSquareTestResultForATheta.alpha == bottomAlphaRange or minChiSquareTestResultForATheta.alpha == topAlphaRange):
                print("Expand alpha range, stopping loop")
                break"""
            i += 1
        #Test if this every hits the extreme values of theta, if so expand the range
        """if(minChiSquareTestResultForATheta.angularDiameter/((1/1000)*(1/60)*(1/60)*(np.pi/180)) == bottomThetaRange or minChiSquareTestResultForATheta.angularDiameter/((1/1000)*(1/60)*(1/60)*(np.pi/180)) == topThetaRange):
                print("Expand theta range, stopping loop")
                break"""
        n += 1
    print('done')
    i=0
    min = chiSquareTestResult(0,0,1e10)
    while i < np.size(chiSquareTestValues):
        if chiSquareTestValues[i].chiSquare < min.chiSquare:
            min = chiSquareTestValues[i]
            print('computing')
        i += 1
    print(min)

    

    theta = min.angularDiameter
    alpha = min.alpha
    return theta, alpha


thetaMilliarcSeconds = 0 #get this once I run the limbdarkened fit