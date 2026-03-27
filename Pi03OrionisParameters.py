import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sympy as sp
import oifitsTools

#angularDiameter=np.float64(7.3255347215650886e-09), alpha=np.float64(0.1), chiSquare=np.float64(2.0436841620921675)

""" From Gaia EDR3, data about pi03 Ori to check results with
        Distance: 8.016 parsecs
        Teff: 6154 [6150, 6158] K
        log(g): 4.031 [4.025, 4.044] log(cm s-2)

        From SPORES catalogue, v2.1.0
        Distance = 8.0684 parsecs
        Parallax = 123.94 mas from 2007 tho
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

angularDiameterRadian = oifitsTools.convertMilliArcSecToRadian(1.52)
angularDiameterRadianErr = .0001 # have to still calculate this, filler number

# From SPORES Catalogue
parallaxRadian = oifitsTools.convertMilliArcSecToRadian(123.94)
parallaxRadianErr = 0.17

# From Gaia data, from Gaia FGK benchmark stars paper, doi: 10.1051/0004-6361/202347136 
bolometricFlux = 139.928e-11 # Watt per square meter
luminositySolar = 2.816 # solar luminosity
luminosityWatts = luminositySolar*solarLuminosityWatts
massSolar = 1.280 # solar mass
massKg = 1.280*solarMassKg

# to get standard deviation of distance, do parallax + sigma to get distance1 and parallax-sigma to get distance 2. One sigma for distance is the difference of the two
distanceParsecs = (1 / ((124.62)*(1/1000))) #distance from earth in parsecs
print("Distance from Earth:", distanceParsecs, 'parsecs')
distanceKm = 30856775812800*distanceParsecs

# How to find the radius from the distance and angular diameter:
# R = d tan (a/2), where R is the radius of the star in km, d is the distance to the star in km and a is the angular size of the
#   star in radians
radiusKm = distanceKm * np.tan(angularDiameterRadian/2)
radiusS = radiusKm*(1/solarRadiusKm) #convert to solar radii 
print("Radius:", radiusS, 'solar radii')

# How to calculate effective temperature with angular diameter and the bolometric flux
# effT = ((4Fbol)/(sigma*(thetaLD)^2))^(1/4), where sigma is the Stefan-Boltzmann constant
stefanBoltzmannConstant = 5.670374419e-8 # W m^2 K^-4
effTemp = ((4*bolometricFlux)/(stefanBoltzmannConstant*((angularDiameterRadian)**2)))**(1/4) # try erg cm-2 s-1, milliarcseconds
#Teff = 2341(Fbol/angularDiameterMilliarcSeconds**2)**(0.25)
print("Effective Temperature: ", effTemp, " Kelvin")
# effT = (L/(4piR^2*sigma))^(-1/4) where sigma is the Stefan-Boltzmann constant
effTemp = (luminosityWatts/(4*np.pi*((radiusKm*1000)**2)*stefanBoltzmannConstant))**(1/4)
print("Effective Temperature, ", effTemp, "Kelvin")

gravitationalConstant = 6.67408e-11 #m^3 kg^-1 s^-2

surfaceGravityM = (gravitationalConstant*massKg)/((radiusKm*1000)**2)
surfaceGravityCm = surfaceGravityM*(100)
print("Surface Gravity log(g):", np.log10(surfaceGravityCm), "cm s^-2")

# How to find rotational period for my star?

