import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

def convertMilliArcSecToRadian(num):
    return num*((1/1000)*(1/60)*(1/60)*(np.pi/180))

angularDiameter = 2.335 #milliarcseconds


# How to find the distance if you have the parallax angle
# distance = (radius of the two places you're comparing the parallax from ie. earth at one point in the orbit and the opposite side) 
#      / tan(parallax angle)
# distance in parsecs = 1 / P, where P is the parallax angle in arc seconds
# Find this from Gaia data directly
# The parallex of σ GEMINORUM is 27.121043168963787 milliarcseconds according to Gaia data
distance = (1 / ((27.121043168963787)*(1/1000))) #distance from earth in parsecs
print(distance, 'parsecs')

# How to find the radius from the distance and angular diameter:
# R = d tan (a/2), where R is the radius of the star in km, d is the distance to the star in km and a is the angular size of the
#   star in arcseconds
distanceKm = 30856775812800*distance #distance in km
radiusKm = distanceKm * np.tan(convertMilliArcSecToRadian(angularDiameter)/2) #angular diameter needs to be in radians for python...
radiusS = radiusKm*(1/695700) #convert to solar radii (1 solar radius/695700km)
print(radiusS, 'solar radii')
#radius should be 10.1 solar radii

# How to calculate effective temperature with angular diameter and the bolometric flux
# Teff = ((4Fbol)/(sigma*(thetaLD)^2))^(1/4), where sigma is the Stefan-Boltzmann constant
# need to find the bolometric flux from another paper, units should be erg s^-1 cm^-2
# can't find this online...
