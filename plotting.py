import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import importlib

def plotting(wave, smooth_csre_spec, errs, name):
    """
    Makes some nice plots.
    """
    
    # loading parameter file for the exoplanet system
    params = importlib.import_module(f'{name}_params')
    
    # loading the Rp and Rs from CSRE params
    Rp, Rs, numtran = params.csre_params['Rp'], params.csre_params['Rs'], params.csre_params['numtran'] 
    
    # fixed transit depth for comparison
    # includes no absorption effects from planet OR star!
    fixed_depth = (Rp / Rs)**2 * 1e6 # ppm
    
    fig, ax = plt.subplots(figsize = (9,4))
    
    x,y,yerr = errs
    ax.errorbar(x, y, yerr/np.sqrt(numtran), None, 'o', ms = 2, c = 'black', label = f'{int(numtran)} transits', capsize = 2.5, alpha = 0.75)
    ax.plot(wave, smooth_csre_spec, c = 'seagreen', label = 'CSRE', alpha = 0.75)
    ax.axhline(fixed_depth, c = 'black', ls = 'dotted', label = 'fixed transit depth')
    
    # Labels
    fig.suptitle('Chromatic Stellar Radii Effect')
    ax.set_xlabel('Wavelength (um)')
    ax.set_ylabel('Transit depth (ppm)')
    ax.legend(title = f'{name}',fontsize='small', fancybox=True, facecolor = 'lightgray')

    # setting plot limits for better visualizations
    ax.set_ylim(np.min(smooth_csre_spec) - 1, fixed_depth + 1)
    fig.show()

    return fig, ax