import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sympy as sp
import oifitsTools
import random

#angularDiameter=np.float64(7.3255347215650886e-09), alpha=np.float64(0.1), chiSquare=np.float64(2.0436841620921675)

""" From Gaia EDR3, data about pi03 Ori to check results with
        Distance: 8.016 parsecs
        Teff: 6154 [6150, 6158] K
        log(g): 4.031 [4.025, 4.044] log(cm s-2)

        From SPORES catalogue, v2.1.0
        Distance = 8.0684 parsecs
        Parallax = 123.94 mas from 2007 tho from https://ui.adsabs.harvard.edu/link_gateway/2007A&A...474..653V/doi:10.1051/0004-6361:20078357
        ParallaxErr = .17 mas
        Teff = 6443 K from 2022
        TeffErr = 14
        Radius = 1.321 solar radiuses
        Theta = 1.523 mas (no )

        PMOIRed, Antoine Me'rand
"""
# Solar parameters (for conversion), from internet
solarLuminosityWatts = 3.828e26 # Watts
solarRadiusKm = 695700 # Km
solarMassKg = 1.988e30 # Kg

angularDiameterRadian = oifitsTools.convertMilliArcSecToRadian(1.49)
angularDiameterRadianErr = .00000000001 # have to still calculate this, filler number

# From Gaia data, from Gaia FGK benchmark stars paper, doi: 10.1051/0004-6361/202347136 
parallax = 124.62
parallaxErr = 0.22
bolometricFlux = 139.928e-11 # Watt per square meter
bolometricFluxErr = 2.943e-11
logg = 4.31 #cm/s^2
loggErr = 0.01 #cm/s^2


luminositySolar = 2.816 # solar luminosity
luminosityWatts = luminositySolar*solarLuminosityWatts
#massSolar = 1.280 # solar mass
#massSolarKg = massSolar*solarMassKg

distances = []
for n in range(0,5000):
        parallaxSample = random.gauss(parallax, parallaxErr)
        distances.append((1 / ((parallaxSample)*(1/1000))))
distanceParsecs = np.mean(distances)
distanceParsecsErr = np.std(distances)

print("Distance from Earth:", distanceParsecs, '±', distanceParsecsErr, 'parsecs')
distanceKm = 30856775812800*distanceParsecs
distanceKmErr = 30856775812800*distanceParsecsErr

# How to find the radius from the distance and angular diameter:
# R = d tan (a/2), where R is the radius of the star in km, d is the distance to the star in km and a is the angular size of the
#   star in radians
radii = []
for n in range(0,5000):
        distanceKmSample = random.gauss(distanceKm, distanceKmErr)
        angularDiameterRadianSample = random.gauss(angularDiameterRadian, angularDiameterRadianErr)
        radiusKmSample = distanceKmSample * np.tan(angularDiameterRadianSample/2)
        radii.append(radiusKmSample*(1/solarRadiusKm))
radiusS = np.mean(radii)
radiusSErr = np.std(radii)
print("Radius:", radiusS, '±', radiusSErr, 'solar radii')
radiusKm = radiusS*solarRadiusKm
radiusKmErr = radiusSErr*solarRadiusKm
#print("Radius:", radiusKm, '±', radiusKmErr, 'km')

# How to calculate effective temperature with angular diameter and the bolometric flux
# effT = ((4Fbol)/(sigma*(thetaLD)^2))^(1/4), where sigma is the Stefan-Boltzmann constant
stefanBoltzmannConstant = 5.670374419e-8 # W m^2 K^-4
temps = []
for n in range(0, 5000):
        bolometricFluxSample = random.gauss(bolometricFlux, bolometricFluxErr)
        angularDiameterRadianSample = random.gauss(angularDiameterRadian, angularDiameterRadianErr)
        effTempSample = ((4*bolometricFluxSample)/(stefanBoltzmannConstant*((angularDiameterRadianSample)**2)))**(1/4)
        temps.append(effTempSample)
effTemp = np.mean(temps)
effTempErr = np.std(temps)
print("Effective Temperature: ", effTemp, '±', effTempErr, " Kelvin")
# effT = (L/(4piR^2*sigma))^(-1/4) where sigma is the Stefan-Boltzmann constant
#effTemp = (luminosityWatts/(4*np.pi*((radiusKm*1000)**2)*stefanBoltzmannConstant))**(1/4)
#print("Effective Temperature, ", effTemp, "Kelvin")


gravitationalConstant = 6.67408e-11 #m^3 kg^-1 s^-2
masses = []
for n in range(0,5000):
        surfaceGravitySample = random.gauss(logg, loggErr)
        radiusKmSample = random.gauss(radiusKm, radiusKmErr)
        massSample = ((10**(surfaceGravitySample)/100)*(radiusKm*1000)**2)/gravitationalConstant
        masses.append(massSample)
massKg = np.mean(masses)
massKgErr = np.std(masses)
massS = massKg * (1/solarMassKg)
massSErr = massKgErr * (1/solarMassKg)
print('Mass: ', massS, '±', massSErr, ' solar mass') 


# Upper limits on rotation period
vsini = 17.3 # km / s
vsiniErr = 0.2

periods = []
for n in range(0,5000):
        vsiniSample = random.gauss(vsini, vsiniErr)
        radiusSample = random.gauss(radiusKm, radiusKmErr)
        periods.append(((2*np.pi*radiusKmSample)/vsiniSample) / (60*60*24))
ProtSini = np.mean(periods)
ProtSiniErr = np.std(periods)
print("Rotational Period: ", ProtSini, '±', ProtSiniErr, ' days')
print(ProtSini/np.sin(np.pi/10000))
print(ProtSini/np.sin(np.pi/4))
print(ProtSini/np.sin(np.pi/3))
print(ProtSini/np.sin(np.pi/2))

# Age of star, log t = 1/n(log(P_rot)-log(a)-b*log(B-V)-0.4)
a = 0.7725
aErr = .011
b = 0.601 
bErr = 0.024
n = 0.5189
nErr = 0.0070
B = 3.630
BErr = .007
V = 3.190
VErr = 0.009

ages = []
for i in range(0,5000):
        aSamp = random.gauss(a, aErr)
        bSamp = random.gauss(b, bErr)
        nSamp = random.gauss(n, nErr)
        Bsamp = random.gauss(B, BErr)
        Vsamp = random.gauss(V, VErr)
        ProtSiniSamp = random.gauss(ProtSini, ProtSiniErr)
        ages.append(10**((1/n)*(np.log(ProtSiniSamp)-np.log(aSamp)-bSamp*np.log(Bsamp-Vsamp)-0.4)))
minAge = np.mean(ages)
minAgeErr = np.std(ages)
print("Min Age:", minAge, '±', minAgeErr, " ")
