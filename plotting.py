import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import importlib

def plotting(wave, smooth_csre_spec, errs, name, numtran = 1, chisq = None, numsig = None, mu = None):
    """
    Makes some nice plots.
    """
    
    # loading parameter file for the exoplanet system
    params = importlib.import_module(f'{name}_params')
    
    # loading the Rp and Rs from CSRE params
    Rp, Rs = params.csre_params['Rp'], params.csre_params['Rs']
    
    # fixed transit depth for comparison
    # includes no absorption effects from planet OR star!
    fixed_depth = (Rp / Rs)**2 * 1e6 # ppm
    
    fig, ax = plt.subplots(figsize = (9,4))
    ax.plot(wave, smooth_csre_spec, c = 'seagreen', label = 'CSRE', alpha = 0.75)
    ax.axhline(fixed_depth, c = 'black', ls = 'dotted', label = 'fixed transit depth')
    
    ### Errors ###
    x,y,yerr = errs
    # plots multiple errorbars corresponding for different number of transits
    if isinstance(numtran, list):
        numtran.sort() # makes sure it plots smaller errorbars on top of larger ones
        alphas = np.linspace(0.5,1,len(numtran))
        for i in range(len(numtran)):
            ax.errorbar(x, y, yerr/np.sqrt(numtran[i]), None, 'o', ms = 2, c = 'dimgray', label = f'{int(numtran[i])} transits', capsize = 2.5, alpha = alphas[i]) 
    else:
        ax.errorbar(x, y, yerr/np.sqrt(numtran), None, 'o', ms = 2, c = 'black', label = f'{int(numtran)} transits', capsize = 2.5, alpha = 0.75)


    ### Chisq ###
    # plots the weighted mean and chisq value
    if chisq is not None:
        ax.axhline(mu, ls = 'dashed', c = 'seagreen', alpha = 0.65, label = 'flux weighted mean')
        fig.text(0.6, 0.2, f'$\chi^2 = {int(chisq)}$', fontsize = 12, bbox=dict(boxstyle='round',ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),))
    
    
    ### Labels and limits ###
    fig.suptitle('Chromatic Stellar Radii Effect')
    ax.set_xlabel('Wavelength (um)')
    ax.set_ylabel('Transit depth (ppm)')
    ax.legend(title = f'{name}',fontsize='small', fancybox=True, facecolor = 'lightgray')

    # setting better plot limits
    ax.set_ylim(np.min(smooth_csre_spec) - 1, fixed_depth + 1)
    fig.show()

    return fig, ax