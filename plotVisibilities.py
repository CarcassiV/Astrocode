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

spatialFrequency = []
spatialFrequencyFiveNights = []
i = 0
while i < np.size(oifitsobj.vis2):
    spatialFrequency.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    if(i<5):
        spatialFrequencyFiveNights.append(np.sqrt((oifitsobj.vis2[i].ucoord)**2 + (oifitsobj.vis2[i].vcoord)**2)/oifitsobj.vis2[i].wavelength.eff_wave/1e6)
    i += 1 
#print(spatialFrequency)


"""chiSquareValues = {}

thetaMicroArcSeconds = np.arange(1.7, 2.7, 0.01)
thetaRadians = thetaMicroArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))  #(1 arcsecs / 1000 millisecs)*(1 arcmin / 60 arcsecs)*(1 degree / 60 arcmin)*(pi radians / 180 degrees)
#print(thetaRadians)
while i < np.size(thetaRadians):
    chiSquare = 0
    j = 0
    while j < np.size(spatialFrequency,0):
        k = 0
        while k < np.size(spatialFrequency,1):
            observed = visibilities[j][k]
            expected = ((2*jv(1, np.pi*thetaRadians[i]*spatialFrequency[j][k]*1e6))/(np.pi*thetaRadians[i]*spatialFrequency[j][k]*1e6))**2
            print(np.pi*thetaRadians[i]*spatialFrequency[j][k]*1e6)
            chiSquareValue = ((observed-expected)**2)/expected
            if not np.isnan(chiSquareValue):
                chiSquare += chiSquareValue
            k += 1
        j += 1
    #print('Theta:', thetaRadians[i], ', Chi Squared Value:', chiSquare)
    chiSquareValues.update({thetaRadians[i]: chiSquare})
    i+=1

#print(chiSquareValues) 
print(min(chiSquareValues, key=chiSquareValues.get)/((1/1000)*(1/60)*(1/60)*(np.pi/180))) 
theta = min(chiSquareValues, key=chiSquareValues.get) #assigns theta to the theta that produced the smallest chi squared value
#2.335 milliarc seconds should be the correct angular diameter
# = theta * 1e9"""

flatVis = flatten(visibilities, 8, 8)
flatErr = flatten(visibilitiesError, 8, 8)
flatSpatial = flatten(spatialFrequency, 8, 8)
flatClose = flatten(closurePhases, 5, 8)
flatSpatialFive = flatten(spatialFrequencyFiveNights, 5, 8)


thetas = []
#for each visibility, randomly sample a point on the error bar
for i in range(0, 100):
    sampleVisibilities = []
    for i in range(0, np.size(flatVis)):
        randomVisibilitySample = random.uniform(flatVis[i]-flatErr[i], flatVis[i]+flatErr[i])
        sampleVisibilities.append(randomVisibilitySample)
    #with these points, for each theta from 1.7milliarc second to 2.7 milliarc second, calculate the 
    #chi square value and find the optimal theta for that dataset
    chiSquareValues = {}

    thetaMicroArcSeconds = np.arange(1.7, 2.7, 0.0001)
    thetaRadians = thetaMicroArcSeconds*((1/1000)*(1/60)*(1/60)*(np.pi/180))
    i = 0
    while i < np.size(thetaRadians):
        chiSquare = 0
        j = 0
        while j < np.size(sampleVisibilities):
            observed = sampleVisibilities[j]
            expected = ((2*jv(1, np.pi*thetaRadians[i]*flatSpatial[j]*1e6))/(np.pi*thetaRadians[i]*flatSpatial[j]*1e6))**2
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

x = np.arange(10, 225, .2) #for the visibility squared curve

fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].plot(flatSpatial, flatVis, '.')
ax[0].plot(x, ((2*jv(1, np.pi*theta*x*1e6))/(np.pi*theta*x*1e6))**2)
ax[0].errorbar(flatSpatial, flatVis, yerr=flatErr, fmt = '.')
ax[0].set_ylabel('Visibilities Squared')
#ax[0].set_yscale('log', base=10)

ax[1].plot(flatSpatialFive, flatClose, '.')
ax[1].set_xlabel('Spatial Frequency')
ax[1].set_ylabel('Closure Phases')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()