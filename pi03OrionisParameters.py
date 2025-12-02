import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
import random

def convertMilliArcSecToRadian(num):
    return num*((1/1000)*(1/60)*(1/60)*(np.pi/180))

# uniform disk fit: with 100 trials using the error bar method, I received a theta of 1.4702000000000017 milliarc seconds

thetaMilliarcSeconds = 0 #get this once I run the limbdarkened fit