import numpy as np
import oifits
import matplotlib.pyplot as plt
from scipy.special import jv

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

chiSquareValues = {}

theta = 1.7*4.8481368110954e-3 #theta in microradians (is this the right unit?????)
thetaMax = 2.7*4.8481368110954e-3
dtheta = 1e-6
while theta <= thetaMax:
    chiSquare = 0
    i = 0
    while i < np.size(spatialFrequency,0):
        j = 0
        while j < np.size(spatialFrequency,1):
            O = visibilities[i][j]
            E = ((2*jv(1, np.pi*theta*spatialFrequency[i][j]))/(np.pi*theta*spatialFrequency[i][j]))**2
            chiSquareValue = ((O-E)**2)/E
            if not np.isnan(chiSquareValue):
                chiSquare += chiSquareValue
            j += 1
        i += 1
    print('Theta:', theta, ', Chi Squared Value:', chiSquare)
    chiSquareValues.update({theta: chiSquare})
    theta += dtheta

#print(chiSquareValues) 
print(min(chiSquareValues, key=chiSquareValues.get)/(4.8481368110954e-3)) #prints, in microrad, the theta that produced the smallest chi squared value
theta = min(chiSquareValues, key=chiSquareValues.get) #assigns theta to the theta that produced the smallest chi squared value
#2.335 milliarc seconds should be the correct angular diameter

x = np.arange(10, 225, .2) #for the visibility squared curve

flatVis = flatten(visibilities, 8, 8)
flatErr = flatten(visibilitiesError, 8, 8)
flatSpatial = flatten(spatialFrequency, 8, 8)
flatClose = flatten(closurePhases, 5, 8)
flatSpatialFive = flatten(spatialFrequencyFiveNights, 5, 8)

fig, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(flatSpatial, flatVis, '.')
ax[0].plot(x, ((2*jv(1, np.pi*theta*x))/(np.pi*theta*x))**2)
ax[0].errorbar(flatSpatial, flatVis, yerr=flatErr, fmt = '.')
ax[0].set_ylabel('Visibilities Squared')
#ax[0].set_yscale('log', base=10)

ax[1].plot(flatSpatialFive, flatClose, '.')
ax[1].set_xlabel('Spatial Frequency')
ax[1].set_ylabel('Closure Phases')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()