import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy.interpolate import interp1d

import importlib
import molmass

from main import chisquared

def plotting(name, wave, smooth_csre_spec, errs = None, numtran = 1, chisq = None, numsig = None, mu = None):
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
    
    ### Transit spectrum ###
    ax.plot(wave, smooth_csre_spec, c = 'black', label = 'CSRE', alpha = 1)
    ax.axhline(fixed_depth, c = 'black', ls = 'dotted', label = 'fixed transit depth')
    
    ### Errors ###
    if errs is not None:
        x,y,yerr = errs
        # plots multiple errorbars corresponding for different number of transits
        if isinstance(numtran, list):
            numtran.sort() # makes sure it plots smaller errorbars on top of larger ones
            alphas = np.logspace(np.log10(0.33),np.log10(1),len(numtran))
            for i in range(len(numtran)):
                ax.errorbar(x, y, yerr/np.sqrt(numtran[i]), None, 'o', ms = 2, c = 'slateblue', label = f'{int(numtran[i])} transits', capsize = 2.5, alpha = alphas[i]) 
        else:
            ax.errorbar(x, y, yerr/np.sqrt(numtran), None, 'o', ms = 2, c = 'crimson', label = f'{int(numtran)} transits', capsize = 2.5, alpha = 0.75)

    ### Chisq ###
    # includes chisq value on plot
    if chisq is not None:
        fig.text(0.6, 0.2, f'$\chi^2 = {int(chisq)}$', fontsize = 12, bbox=dict(boxstyle='round',ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),))
        
    ### flux weighted mean ###
    # plots mu
    if mu is not None:
        ax.axhline(mu, ls = 'dashed', c = 'black', alpha = 0.65, label = 'flux weighted mean')
    
    
    ### Labels and limits ###
    #fig.suptitle('Chromatic Stellar Radii Effect')
    ax.set_xlabel('Wavelength (um)')
    ax.set_ylabel('Transit depth (ppm)')
    ax.legend(title = f'{name}',fontsize='small', fancybox=True, facecolor = 'lightgray', loc = 'upper right')

    # setting better plot limits
    ax.set_xlim(np.min(wave), np.max(wave))
    ax.set_ylim(np.min(smooth_csre_spec) - 1, fixed_depth + 1)
    fig.show()

    return fig, ax


def sigma_lines(fig, ax, x, y, numsig):
    """
    Handles the plotting of the lines for 
    1sigma, 2sigma, 3sigma, etc. on the chisq plot.
    
    """
    
    for i in range(len(numsig)):
        # chisq = (ns)^2
        hline_value = numsig[i]**2
        
        # min and max for lines
        xfunc = interp1d(y,x)
        hline_max = xfunc(hline_value)

        ax.hlines(hline_value, 0, hline_max, ls = 'dotted', color = 'black')
        ax.vlines(hline_max, 0, hline_value, ls = 'dotted', color = 'black')

        ax.text(hline_max, hline_value, f'${numsig[i]} \sigma$', horizontalalignment = 'right', verticalalignment='bottom')
        
        print(f"{numsig[i]}-sigma level: {math.ceil(hline_max)} transits")
        
        # using for plot limits
        if i == 0:
            xmin = 0.95 * hline_max
            ymin = 0.95 * hline_value
        if i == len(numsig) - 1:
            
            xmax = 1.05 * hline_max
            ymax = 1.05 * hline_value
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    
    return fig, ax


def plot_depth_changes(name, wave, spec, errs, molecules, N = 5):
    """
    Calculates and plots the change in transit depth according to the analytic equation:
    
    deltadelta = 2*(NH/Rs)*delta           where delta = (Rp/Rs)^2
    
    This is done for each molecule included in the ExoTransmit fit, 
    and they are all plotted with respect to the flux-weighted average (mu)
    of the CSRE, such that they cover the range (mu - deltadelta, mu + deltadelta).
    
    ================== Params ==================
    
    name       : name of system ; must match what is on {name}_params.py
    wave       : array or list of wavelengths (microns)
    spec       : transit spectrum including the CSRE
    errs       : error array output by CSRE function
    molecules  : list of molecules included in your ExoTransmit fit
    N          : number of scale heights to use, should be of order unity (Default = 5)
    
    ================== Returns =================
    
    fig  : figure object
    ax   : axes object
    """
    
    params = importlib.import_module(f'{name}_params') # for loading param values
    
    # initializing dict of molecular masses
    masses_dict = {}
    for molecule in molecules:
        masses_dict[molecule] = molmass.Formula(molecule).mass
    
    # sorting by increasing mass (so smallest deltadeltas get plotted first)
    keys = list(masses_dict.keys())
    values = list(masses_dict.values())
    sorted_value_index = np.argsort(values)
    masses_dict = {keys[i]: values[i] for i in sorted_value_index}
    amu = sp.constants.physical_constants['atomic mass constant'][0] # [kg]

    # constants and param values
    g      = 10**params.pandexo_params['star_logg'] / 100 # [m/s^2]
    T      = params.pandexo_params['star_temp'] # [K]
    k      = sp.constants.k # [J / K]
    Rp, Rs = params.csre_params['Rp'], params.csre_params['Rs'] # [m]
    depth  = (Rp / Rs) ** 2 * 1e6 # [ppm]

    # calculating the change in transit depth due to each molecule
    deltadeltas = []
    for molecule in masses_dict:
        u = masses_dict[molecule] * amu
        H = k * T / (u * g)

        deltadelta = 2 * depth * (N*H) / Rs
        deltadeltas.append(deltadelta)

    # using a fancy color map for the deldels
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Wistia(np.linspace(1,0,int(len(molecules)))))

    mu = chisquared(name, errs)[-1]
    fig, ax = plotting(name, wave, spec, mu = mu)

    # plotting all the deltadeltas from each molecule
    for i in range(len(deltadeltas)):
        ax.fill_between(wave, mu - deltadeltas[i], mu + deltadeltas[i])
    # adding deltadelta label on side
    ax.text(5.35,mu,']', fontsize = 14, verticalalignment = 'bottom')
    ax.text(5.425,mu,r'$\Delta\delta$', fontsize = 14, verticalalignment = 'bottom')
    
    return fig, ax