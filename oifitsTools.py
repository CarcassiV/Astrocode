import numpy as np
from astropy.io import fits
import oifits
import matplotlib.pyplot as plt
from scipy.special import jv, j0, j1
from scipy.stats import chisquare
import random
import scipy.integrate as integrate
from dataclasses import dataclass
from typing import Self
import sympy as sp
import csv

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

def openFile(fileName):

    def flatten(ndarr): #assumes each row is of the same length
        rows = len(ndarr)
        cols = len(ndarr[0])
        flatarr = []
        for i in range(0,rows):
            for j in range(0,cols):
                flatarr.append(ndarr[i][j])
        return flatarr

    hdulist = fits.open(fileName)
    hdulist['OI_ARRAY'].header['OI_REVN'] = 1
    oifitsobj = oifits.open(hdulist)

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
    i = 0
    while i < np.size(oifitsobj.vis2):
        twoDSpatialFrequency.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
        i += 1

    twoDSpatialFrequencyClosurePhases = []
    i = 0
    while i < np.size(oifitsobj.t3):
        twoDSpatialFrequencyClosurePhases.append(np.sqrt((oifitsobj.t3[i].u1coord)**2 + (oifitsobj.t3[i].v1coord)**2)/oifitsobj.t3[i].wavelength.eff_wave/1e6)
        i += 1
    
    visibilitiesSquared = flatten(twoDVisibilities)
    visibilitiesSquaredErr = flatten(twoDVisibilitiesError)
    spatialFrequencies = flatten(twoDSpatialFrequency)
    closurePhases = flatten(twoDClosurePhases)
    closurePhasesErr = flatten(twoDClosurePhasesErrors)
    spatialFrequenciesClosurePhases = flatten(twoDSpatialFrequencyClosurePhases)
    
    return visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, spatialFrequenciesClosurePhases, closurePhases, closurePhasesErr


#31.48s for 10 trials with only calculating up to 3 chi squared
#1min 18secs for 10 trials from before
#Nearly 2.5 times faster
#uniform disk model angular diameter
def uniformDiskFitErrorBarTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    # Make CSV file to store trial results, that way even if code crashes, we have some results.
    with open('uniformDiskFit.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Trial', 'Angular Diameter', 'Chi-Squared Value'])

    thetas = []
    #for each visibility, randomly sample a point on the error bar
    for n in range(0, numberOfTrials):
        """sampleVisibilities = []
        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i])
            sampleVisibilities.append(randomVisibilitySquaredSample)
        chiSquareValues = {}

        sampleVisibilities = visibilitiesSquared"""

        bottomThetaRange = 1.8
        topThetaRange = 2.8

        thetaMilliArcSeconds = np.arange(bottomThetaRange, topThetaRange, 0.001)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
        i = 0
        while i < np.size(thetaRadians):
            chiSquare = 0
            j = 0
            while j < np.size(spatialFrequencies):
                observed = sampleVisibilities[j]
                expected = ((2*j1(np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))/(np.pi*thetaRadians[i]*spatialFrequencies[j]*1e6))**2
                chiSquareValue = ((observed-expected)**2)/expected
                if not np.isnan(chiSquareValue):
                    chiSquare += chiSquareValue
                """if chiSquare > 3:
                    break"""
                j += 1
            print('Theta:', thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), ', Chi Squared Value:', chiSquare)
            chiSquareValues.update({thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)): chiSquare})
            i += 1
        thetas.append(min(chiSquareValues, key=chiSquareValues.get))
        print("Trial Number:", n, "Theta:", min(chiSquareValues, key=chiSquareValues.get)) 
        
        #Test if this ever hits the extreme values of theta, if so expand the range. Try to make range as small as possible to save time
        if(min(chiSquareValues, key=chiSquareValues.get) == bottomThetaRange or min(chiSquareValues, key=chiSquareValues.get) == topThetaRange):
            print("Expand theta range, stopping loop")
            break

        with open('uniformDiskFit.csv', 'a', newline='') as file: #Could open this less often, add every 5 trials for example, to speed it up
            writer = csv.writer(file)
            writer.writerow([n, min(chiSquareValues), min(chiSquareValues)])

        n += 1
    #take the average of all the theta
    print("done")
    i = 0
    sum = 0
    while i < np.size(thetas):
        sum += thetas[i]
        i += 1
    theta = sum/np.size(thetas)
    print(theta)

    #The error is calculated from the standard deviation of the trials
    thetaError = np.std(thetas)

    return theta, thetaError, chiSquareValues


"""
    Ideas to make the code more efficient:
    - Write to the csv file less often?
    - Explore chi squared libraries? scipy.stats.chisquare? For the goodness of fit test? https://www.statology.org/how-to-conduct-chi-square-tests-scipy/ 
    - Write down all the steps and see if there are unnecessary arrays or the like

"""


def uniformDiskFitErrorBarTestWithCenterVisibility(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    # Make CSV file to store trial results, that way even if code crashes, we have some results.
    with open('trialData.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Trial', 'Angular Diameter', 'Visibility Value', 'Chi-Squared Value'])

    minChiSquareTestResults = []
    
    for n in range(0, numberOfTrials):
        sampleVisibilitiesSquared = []
        minChiSquareTestResultForATrial = chiSquareTestResult(0,.5,1e10)

        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i])
            sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)


        bottomThetaRange = 1.4
        topThetaRange = 1.55
        
        # Add another parameter V0, which varies from .95-1.05 (check if hits edge). expected = V0*(current function)
        
        thetaMilliArcSeconds = np.arange(bottomThetaRange, topThetaRange, 0.001)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))

        minChiSquareTestResultForATheta = chiSquareTestResult(0.0,.5,1e10)

        bottomVisibilityRange = .95
        topVisibilityRange = 1.05

        visibilityRange = np.arange(bottomVisibilityRange, topVisibilityRange, 0.01)

        i = 0
        while i < np.size(thetaRadians): #are there libraries that already run chi square tests? More efficient than mine perhaps?
            j = 0
            while j < np.size(visibilityRange):
                k = 0
                chiSquare = 0
                while k < np.size(sampleVisibilitiesSquared):
                    observed = sampleVisibilitiesSquared[k]
                    expected = visibilityRange[j]*((2*j1(np.pi*thetaRadians[i]*spatialFrequencies[k]*1e6))/(np.pi*thetaRadians[i]*spatialFrequencies[k]*1e6))**2
                    chiSquareValue = ((observed-expected)**2)/expected
                    if not np.isnan(chiSquareValue):
                        chiSquare += chiSquareValue
                    k += 1
                #print("Theta value:", thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), "visibility value:", visibilityRange[j], "chisquare:", chiSquare)
                if (chiSquare < minChiSquareTestResultForATheta.chiSquare):
                    minChiSquareTestResultForATheta = chiSquareTestResult(thetaRadians[i], visibilityRange[j], chiSquare)
                j += 1
            if(minChiSquareTestResultForATheta.chiSquare < minChiSquareTestResultForATrial.chiSquare):
                minChiSquareTestResultForATrial = minChiSquareTestResultForATheta
            i += 1
        #print(minChiSquareTestResultForATrial)
        minChiSquareTestResults.append(minChiSquareTestResultForATrial)

        with open('trialData.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([n, minChiSquareTestResultForATrial.angularDiameter/((1/1000)*(1/60)*(1/60)*(np.pi/180)), minChiSquareTestResultForATrial.alpha, minChiSquareTestResultForATrial.chiSquare])

        n += 1

    print('done')
    i=0
    min = chiSquareTestResult(0,0,1e10)
    while i < np.size(minChiSquareTestResults):
        if minChiSquareTestResults[i].chiSquare < min.chiSquare:
            min = minChiSquareTestResults[i]
            print('computing')
        i += 1
    print(min)

    theta = min.angularDiameter
    visibility = min.alpha
    return theta, visibility


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


# table with spatial frequency, visibilities squared, expected, chi square

#limb-darkened model angular diameter
def limbdarkenedThetaTest(visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, numberOfTrials):
    # Make CSV file to store trial results, that way even if code crashes, we have some results.
    with open('limbdarkened.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Trial', 'Angular Diameter', 'Alpha Value', 'Chi-Squared Value'])

    """with open('chiSquare.csv', 'w', newline='') as fileOne:
        writer = csv.writer(fileOne)
        writer.writerow(['Theta', 'Alpha', 'Spatial Frequency', 'Visibilities Squared', 'Expected Value', 'Chi-Squared Value'])"""

    minChiSquareTestResults = []
    
    for n in range(0, numberOfTrials):
        sampleVisibilitiesSquared = []
        minChiSquareTestResultForATrial = chiSquareTestResult(0,.5,1e10)

        for i in range(0, np.size(visibilitiesSquared)):
            randomVisibilitySquaredSample = random.gauss(visibilitiesSquared[i], visibilitiesSquaredErr[i])
            sampleVisibilitiesSquared.append(randomVisibilitySquaredSample)

        sampleVisibilitiesSquared = visibilitiesSquared

        # calculate expected and put into an array. Don't need to recalculate each time.

        bottomThetaRange = 1.6
        topThetaRange = 2.6
        
        # Add another parameter V0, which varies from .95-1.05 (check if hits edge). expected = V0*(current function)
        
        thetaMilliArcSeconds = np.arange(bottomThetaRange, topThetaRange, 0.01)
        thetaRadians = thetaMilliArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))

        minChiSquareTestResultForATheta = chiSquareTestResult(0.0,.5,1e10)

        bottomAlphaRange = 0.00
        topAlphaRange = 0.4

        alpha = np.arange(bottomAlphaRange, topAlphaRange, 0.01)

        i = 0
        while i < np.size(thetaRadians): #are there libraries that already run chi square tests? More efficient than mine perhaps?
            j = 0
            while j < np.size(alpha):
                k = 0
                chiSquare = 0
                while k < np.size(sampleVisibilitiesSquared):
                    #if(spatialFrequencies[k] < 160):
                    observed = sampleVisibilitiesSquared[k]
                    function = lambda r : ((1-r**2)**(alpha[j]/2))*j0(np.pi*thetaRadians[i]*r*spatialFrequencies[k]*1e6)*r
                    integral, err = integrate.quad(function, 0, 1) 
                    expected = ((alpha[j]+2)*integral)**2
                    chiSquareValue = ((observed-expected)**2)/expected
                    """with open('chiSquare.csv', 'a', newline='') as fileOne:
                        writer = csv.writer(fileOne)
                        writer.writerow([thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)) , alpha[j], spatialFrequencies[k], observed, expected, chiSquareValue])
                    """
                    if not np.isnan(chiSquareValue):
                        chiSquare += chiSquareValue
                    #if(chiSquare > 5 or chiSquare > minChiSquareTestResultForATheta.chiSquare):
                        #break
                    k += 1
                print("Theta value:", thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), "alpha value:", alpha[j], "chisquare:", chiSquare)
                with open('limbdarkened.csv', 'a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow([n, thetaRadians[i]/((1/1000)*(1/60)*(1/60)*(np.pi/180)), alpha[j], chiSquare])
                if (chiSquare < minChiSquareTestResultForATheta.chiSquare): #I compare every theta and alpha pair, is it more efficient to store them all and sort/compare at end to find min?
                    minChiSquareTestResultForATheta = chiSquareTestResult(thetaRadians[i], alpha[j], chiSquare)
                j += 1
            if(minChiSquareTestResultForATheta.chiSquare < minChiSquareTestResultForATrial.chiSquare):
                minChiSquareTestResultForATrial = minChiSquareTestResultForATheta
            i += 1
        print(minChiSquareTestResultForATrial)
        minChiSquareTestResults.append(minChiSquareTestResultForATrial)

        """with open('limbdarkened.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([n, minChiSquareTestResultForATrial.angularDiameter/((1/1000)*(1/60)*(1/60)*(np.pi/180)), minChiSquareTestResultForATrial.alpha, minChiSquareTestResultForATrial.chiSquare])
"""
        n += 1

    print('done')
    i=0
    min = chiSquareTestResult(0,0,1e10)
    while i < np.size(minChiSquareTestResults):
        if minChiSquareTestResults[i].chiSquare < min.chiSquare:
            min = minChiSquareTestResults[i]
            print('computing')
        i += 1
    print(min)    

    theta = min.angularDiameter
    alpha = min.alpha
    return theta, alpha