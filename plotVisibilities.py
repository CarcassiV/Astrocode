import numpy as np
import oifits
import matplotlib.pyplot as plt
from scipy.special import jv
import random

#To work on:
#How do I access the closure phases array?
#What is effective wavelength vs. wavelength?
#Are the squared visibilities corresponding to the correct wavelengths? How can I tell?
#Is the angular diameter separate from the wavelength over the baseline? But the angular diameter is in the units of wavelength over baseline?
#Why does the angular diameter in the summer school paper so much smaller than Korolik, when both are baseline/wavelength?
def flatten(ndarr, rows, cols):
    flatarr = []
    for i in range(0,rows):
        for j in range(0,cols):
            flatarr.append(ndarr[i][j])
    return flatarr

oifitsobj = oifits.open('2011Dec07.17ms.sigGem.oifits')

visibilities = []
i = 0
while i < np.size(oifitsobj.vis2):
    visibilities.append(np.ma.getdata(oifitsobj.vis2[i].vis2data))
    i += 1

visibilitiesError = []
i = 0
while i < np.size(oifitsobj.vis2):
    visibilitiesError.append(np.ma.getdata(oifitsobj.vis2[i].vis2err))
    i += 1

closurePhases = []
i = 0
while i < np.size(oifitsobj.t3):
    closurePhases.append(np.ma.getdata(oifitsobj.t3[i].t3phi))
    i += 1

closurePhasesErrors = []
i = 0
while i < np.size(oifitsobj.t3):
    closurePhasesErrors.append(np.ma.getdata(oifitsobj.t3[i].t3phierr))
    i += 1

spatialFrequency = []
spatialFrequencyFiveNights = []
i = 0
while i < np.size(oifitsobj.vis2):
    spatialFrequency.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    if(i<5):
        spatialFrequencyFiveNights.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    i += 1 
print(spatialFrequency)

oneDVis = flatten(visibilities, 8, 8)
oneDVisErr = flatten(visibilitiesError, 8, 8)
oneDSpatial = flatten(spatialFrequency, 8, 8)
oneDClose = flatten(closurePhases, 5, 8)
oneDCloseErr = flatten(closurePhasesErrors, 5, 8)
oneDSpatialFive = flatten(spatialFrequencyFiveNights, 5, 8)
"""
thetas = []
#for each visibility, randomly sample a point on the error bar
for i in range(0, 1000):
    sampleVisibilities = []
    for i in range(0, np.size(oneDVis)):
        randomVisibilitySample = random.uniform(oneDVis[i]-oneDVisErr[i], oneDVis[i]+oneDVisErr[i]) #should be normal distribution 
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
            expected = ((2*jv(1, np.pi*thetaRadians[i]*oneDSpatial[j]*1e6))/(np.pi*thetaRadians[i]*oneDSpatial[j]*1e6))**2
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
#repeat the same process for the limb-darkened angular diameter?
#okay idk how to use that equation...

#make a dictionary where each visibility will correspond with a certain spatial frequency
visDict = {}
for i in range(0, np.size(oneDSpatial)):
    if not np.isnan(oneDVis[i]):
        visDict.update({oneDVis[i]:oneDSpatial[i]})

thetas = []
#for each visibility, randomly sample a point on the error bar
for i in range(0, 1000):
    sampleVisibilities = []
    for i in range(0, np.size(oneDVis)):
        randomVisibilitySample = oneDVis[random.randint(0, np.size(oneDVis)-1)]
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

#repeat 500 times
#take the average of all the theta
print("done")
i = 0
sum = 0
while i < np.size(thetas):
    sum += thetas[i]
    i += 1
theta = sum/np.size(thetas)
print(theta) """

theta = 2.335*((1/1000)*(1/60)*(1/60)*(np.pi/180))

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].plot(oneDSpatial, oneDVis, '.')
ax[0].plot(x, ((2*jv(1, np.pi*theta*x*1e6))/(np.pi*theta*x*1e6))**2)
ax[0].errorbar(oneDSpatial, oneDVis, yerr=oneDVisErr, fmt = '.')
ax[0].set_ylabel('Visibilities Squared')
#ax[0].set_yscale('log', base=10)

ax[1].plot(oneDSpatialFive, oneDClose, '.')
ax[1].errorbar(oneDSpatialFive, oneDClose, yerr=oneDCloseErr, fmt = '.')
ax[1].set_xlabel('Spatial Frequency (Mλ)')
ax[1].set_ylabel('Closure Phases (degrees)')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()