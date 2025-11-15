import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
import random

def convertMilliArcSecToRadian(num):
    return num*((1/1000)*(1/60)*(1/60)*(np.pi/180))

angularDiameter = 2.335 #milliarcseconds
angularDiameterErr = .007 #milliarcSeconds
angularDiameterRadians = convertMilliArcSecToRadian(angularDiameter)
angularDiameterErrRadians = convertMilliArcSecToRadian(angularDiameterErr)

# How to find the distance if you have the parallax angle
# distance = (radius of the two places you're comparing the parallax from ie. earth at one point in the orbit and the opposite side) 
#      / tan(parallax angle)
# distance in parsecs = 1 / P, where P is the parallax angle in arc seconds
# Find this from Gaia data directly
# The parallex of σ GEMINORUM is 27.121043168963787 milliarcseconds according to Gaia data
distance = (1 / ((27.121043168963787)*(1/1000))) #distance from earth in parsecs
print("Distance from Earth:", distance, 'parsecs')

# How to find the radius from the distance and angular diameter:
# R = d tan (a/2), where R is the radius of the star in km, d is the distance to the star in km and a is the angular size of the
#   star in arcseconds
distanceKm = 30856775812800*distance #distance in km
radii = []
for i in range(0, 100):
    sampleAngularDiameter = random.gauss(angularDiameterRadians, angularDiameterErrRadians)
    radii.append(distanceKm * np.tan(convertMilliArcSecToRadian(angularDiameter)/2))
radiusKm = sum(radii) / len(radii)
radiusS = radiusKm*(1/695700) #convert to solar radii (1 solar radius/695700km)
print("Radius:", radiusS, 'solar radii')
#radius should be 10.1 solar radii

# How to calculate effective temperature with angular diameter and the bolometric flux
# Teff = ((4Fbol)/(sigma*(thetaLD)^2))^(1/4), where sigma is the Stefan-Boltzmann constant
# L = 4piR^2Fbol, L/(4piR^2)=Fbol
# L = 39 ± 2 (Roettenbacher et al. 2015)
luminosity = 39 # watts
luminosityErr = 2
radiusM = radiusKm * (1000)

bolometricFluxes = []
for i in range(0, 100): #will have to account for the error in radius but didn't invest enough skill points into stats yet
    sampleLuminosity = random.gauss(luminosity, luminosityErr)
    bolometricFluxes.append(sampleLuminosity / (4*np.pi*radiusM**2))
bolometricFluxM = sum(bolometricFluxes) / len(bolometricFluxes)
print("Bolometrix Flux:", bolometricFluxM, "watts m^-2")

effectiveTemperatures = []
bolometricFluxCm = bolometricFluxM * (10000) #(1000cm^2/1m^2)
stefanBoltzmannConstant = 5.670367e-8 #W m^−2 K^−4
for i in range(0,100):
    sampleAngularDiameter = random.gauss(angularDiameterRadians, angularDiameterErrRadians)
    effectiveTemperatures.append(((4*bolometricFluxCm)/(stefanBoltzmannConstant*(sampleAngularDiameter**2)))**(1/4))
effectiveTemperature = sum(effectiveTemperatures) / len(effectiveTemperatures)
print("Effective Temperature:", effectiveTemperature, "K")

# surface gravity = G*M/R^2, where G is the gravitational constant, M is the mass of the star in Kg and R is the radius of the star in Km
# Mass = 1.28 ± 0.07 solar masses (Roettenbacher et al. 2015)
massSolarMass = 1.28 
massSolarMassErr = .07
massKg = 1.28 * 1.98847E+30 #1.98847E+30kg in one solar mass
massKgErr = .07 * 1.98847E+30

surfaceGravities = []
for i in range(0, 100):
    sampleMass = random.gauss(massKg, massKgErr)
    surfaceGravities.append()
