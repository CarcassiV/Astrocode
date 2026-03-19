import numpy as np
from astropy.io import fits
import oifits
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, j0, j1
import scipy.integrate as integrate
import scipy.stats as stats
import sympy as sp
import random
from oifitsTools import openFile


visibilitiesSquaredCHARA, visibilitiesSquaredErrCHARA, spatialFrequenciesCHARA, spatialFrequenciesClosurePhasesCHARA, closurePhasesCHARA, closurePhasesErrCHARA = openFile('CHARADataPiOri/MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')
print(visibilitiesSquared)

visibilitiesSquaredVLTI, visibilitiesSquaredVLTI, spatialFrequenciesVLTI, spatialFrequenciesVLTI, spatialFrequenciesClosurePhasesVLTI, closurePhasesVLTI, closurePhasesErrVLTI = openFile()
