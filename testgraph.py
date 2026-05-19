import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import jv, j0, j1
import scipy.integrate as integrate
from oifitsTools import openFile
from matplotlib import rc

# Following guide from https://jonchar.net/notebooks/matplotlib-styling/

def stylize_axes(ax, title, xlabel, ylabel):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    


visibilitiesSquared, visibilitiesSquaredErr, spatialFrequencies, spatialFrequenciesClosurePhases, closurePhases, closurePhasesErr = openFile('CHARADataPiOri/MIRC_L2.2025Sep21.pi03_Ori.MIRCX_IDL.RMR_deepedge.AVG5m.oifits')

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':12})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11,11))

#xticks = np.linspace(0, 230, 50)
#yticks = np.linspace(0, 1.1, .2)
stylize_axes(ax[0], 'Pi03 Orionis Visibility', 'Spatial Frequency', 'Visibility Squared')
stylize_axes(ax[0], 'Pi03 Orionis Closure Phases', 'Spatial Frequency', 'Closure Phases (degrees)')

#ax[0].set_ylim(0, 1.05)
ax[1].set_ylim(-180, 180)

# Plot visibility data
ax[0].errorbar(
    spatialFrequencies, 
    visibilitiesSquared, 
    visibilitiesSquaredErr, 
    color='red', 
    ls = 'none', 
    marker='o',
    capsize = 1, 
    capthick = 1, 
    ecolor = 'red',
    zorder=1
)

# Plot closure phases
ax[1].errorbar(spatialFrequenciesClosurePhases, 
    closurePhases, 
    closurePhasesErr, 
    color='red', 
    ls = 'none', 
    marker='.',
    capsize = 1, 
    capthick = 1,
    ecolor = 'red',
)

# Plot Models
uniformDiskTheta = 1.485*((1/1000)*(1/60)*(1/60)*(np.pi/180)) 
limbdarkenedTheta = 1.5*((1/1000)*(1/60)*(1/60)*(np.pi/180)) 
alpha = 0.11

x = np.arange(10, 215, .5) #for the visibility squared curve

ax[0].plot( #uniform disk model
    x, 
    ((2*j1(np.pi*uniformDiskTheta*x*1e6))/(np.pi*uniformDiskTheta*x*1e6))**2, 
    label='Uniform Disk Model',
    color="blue",
    linewidth=3,
    zorder=2
) 

#uniform disk model
limbDarkenedValues = []
i = 0
while i < np.size(x):
    function = lambda r : ((1-r**2)**(alpha/2))*j0(np.pi*limbdarkenedTheta*r*x[i]*1e6)*r
    integral, err = integrate.quad(function, 0, 1, epsabs=1e-14)
    limbDarkenedValues.append(((alpha+2)*np.abs(integral))**2)
    i+=1

ax[0].plot(
    x, 
    limbDarkenedValues, 
    label='Limb Darkened Model',
    color="green",
    linestyle = '--',
    linewidth=3,
    zorder=3
)

ax[0].legend()
ax[0].set_yscale('log', base=10)
plt.show()

