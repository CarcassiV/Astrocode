import numpy as np
from astropy.io import fits
import oifits
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import random
from oifitsTools import openFile

visibilitiesSquaredCHARA, visibilitiesSquaredErrCHARA, spatialFrequenciesCHARA, spatialFrequenciesClosurePhasesCHARA, closurePhasesCHARA, closurePhasesErrCHARA = openFile('CHARADataPiOri/MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')
print(visibilitiesSquaredCHARA)

visibilitiesSquaredVLTI, visibilitiesSquaredErrVLTI, spatialFrequenciesVLTI, spatialFrequenciesClosurePhasesVLTI, closurePhasesVLTI, closurePhasesErrVLTI = openFile('VLTIDataPiOri/PIONI.2019-11-28T04-54-51.340_oidataCalibrated.fits')
print(visibilitiesSquaredVLTI)

print(len(visibilitiesSquaredVLTI))

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-28T05-13-08.623_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-29T05-18-37.271_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONI.2019-11-29T05-34-15.601_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T08p16p19.858_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T08p44p52.527_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T09p12p58.922_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

v, ve, sf, sfcp, cp, cpe = openFile('VLTIDataPiOri/PIONIER.2011-09-27T09p45p15.562_oidataCalibrated.fits')

visibilitiesSquaredVLTI.extend(v)
visibilitiesSquaredErrVLTI.extend(ve)
spatialFrequenciesVLTI.extend(sf)
spatialFrequenciesClosurePhasesVLTI.extend(sfcp)
closurePhasesVLTI.extend(cp)
closurePhasesErrVLTI.extend(cpe)

print(len(visibilitiesSquaredVLTI))
print(len(visibilitiesSquaredCHARA))

#Visibility Squared plot
fig, ax = plt.subplots(2,1, sharex=True)
print("about to plot")
ax[0].errorbar(spatialFrequenciesCHARA, visibilitiesSquaredCHARA, yerr=visibilitiesSquaredErrCHARA, fmt = '.')
ax[0].errorbar(spatialFrequenciesVLTI, visibilitiesSquaredVLTI, yerr=visibilitiesSquaredErrVLTI, fmt = '*')
ax[1].set_xlabel('Spatial Frequency (Mλ)(baseline/wavelength)')
ax[0].set_ylabel('Visibilities Squared')
ax[0].set_yscale('log', base=10)

#Closure Phases plot
ax[1].errorbar(spatialFrequenciesClosurePhasesCHARA, closurePhasesCHARA, yerr=closurePhasesErrCHARA, fmt = '.')
ax[1].errorbar(spatialFrequenciesClosurePhasesVLTI, closurePhasesVLTI, yerr=closurePhasesErrVLTI, fmt = '*')
ax[1].set_ylabel('Closure Phases (degrees)')

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=.085)
plt.show()