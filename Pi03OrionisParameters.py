import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, j0, j1
import scipy.integrate as integrate
import scipy.stats as stats
import sympy as sp
import oifitsTools

angularDiameterRadian = oifitsTools.convertMilliArcSecToRadian(1.45)
angularDiameterRadianError = .0001

# From Gaia EDR3 data
