import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import importlib

# my functions
from smoothing import adap_smooth
from run_pandexo import run_pandexo
from binning import bin_err

# path to transit spectra
spectra_path = '/Users/coffey/Downloads/kipping/Exo_Transmit/Spectra'

def CSRE(name):
    """
    Calculates the chromatic stellar radii effect (CSRE) for the
    exoplanetary system specified by name. Must have a corresponding
    parameter file of the form 'name.py' in the same directory!
    
    ================== Params ==================
    
    Rp          : radius of the transiting planet [m]
    Rs          : radius of the host star [m]
    Rref        : radius of the reference star used in ExoTransmit runs [m]
    planet_name : name of planet (string)
    star_name   : name of star (string)
    flat        : whether or not to include planet's absorption features (True assumes a flat transit spectrum)
    R           : resolution to smooth final spectrum to, we recommend 100 to match JWST NIRSpec's PRISM grating
    
    ================== Returns =================
    
    wave              : array or list of wavelengths (microns)
    smooth_csre_spec  : transit spectrum including the CSRE, smoothed to chosen resolution
    errs              : (x,y,yerr) of error bars on transit spectrum (ppm)
    """
    
    # loading parameter file for the exoplanet system
    params = importlib.import_module(f'{name}_params')
    
    # loading the CSRE params
    Rp, Rs, Rref, planet_name, star_name, flat, R, numtran = params.csre_params.values()
    
    star_spectrum = np.loadtxt(f'{spectra_path}/transmission_{star_name}.dat', skiprows = 2).T  # ExoTransmit 0.3 - 30 um
    
    # slicing wavelength range to JWST NIRSpec (0.6 - 5.3 um)
    wave_index = slice(np.where(star_spectrum[0] < 0.6e-6)[0][-1], np.where(star_spectrum[0] > 5.3e-6)[0][0])
    wave       = star_spectrum[0][wave_index] * 1e6  # microns
    star_spec  = star_spectrum[1][wave_index] / 100  # convert from %
    
    # Multiplying out the Sun's radius
    star_spec_noRref = np.sqrt(star_spec) * Rref
    
    ### Default mode (assume flat planet spectrum) ###
    if flat == True:
        # combining the spectra for our new (Rp/Rs)^2
        csre_spec  = (Rp / star_spec_noRref)**2 * 1e6  # ppm 

        # smoothing to JWST resolution (R ~ 100 for PRISM)
        smooth_csre_spec = adap_smooth(wave, csre_spec, R)
    
    
    ### include planet absorption ###
    if flat == False:
        planet_spectrum = np.loadtxt(f'{spectra_path}/transmission_{planet_name}.dat', skiprows = 2).T # ExoTransmit 0.3 - 30 um
        planet_spec     = planet_spectrum[1][wave_index] / 100  # convert from %
        
        # Multiplying out the Sun's radius
        planet_spec_noRref = np.sqrt(planet_spec) * Rref
        
        # combining the spectra for our new (Rp/Rs)^2
        csre_spec  = (planet_spec_noRref / star_spec_noRref)**2 * 1e6  # ppm 
    
        # smoothing to JWST resolution (R ~ 100 for PRISM is default)
        smooth_csre_spec = adap_smooth(wave, csre_spec, R)
        
    ### Error bars ###
    errs = errorbars(params, wave, smooth_csre_spec)    
    
    return wave, smooth_csre_spec, errs


def errorbars(params, wave, smooth_csre_spec, R = 25):
    """
    Places errorbars on the final transit spectrum.
    
    ================== Params ==================
    
    params            : parameters for PandExo run (taken from pandexo_params.py input file)
    wave              : array or list of wavelengths (microns)
    smooth_csre_spec  : transit spectrum including the CSRE, smoothed to chosen resolution
    R                 : resolution of error bars (defaulted at 25 to reduce plot clutter)
    
    ================== Returns =================
    
    x     : wavelength positions of error bars (microns)
    y     : transit depths at x wavelength positions (ppm)
    yerr  : array or list of errors on transit spectrum (ppm)
    """
    
    # Running PandExo
    err_wave, err = run_pandexo(params)
    
    # rebinning to a lower R
    new_err_wave, new_err, dellambs = bin_err(err_wave, err, R)
    
    # need to interp to find y value for error bars on the transit spectrum
    interp_func = interp1d(wave, smooth_csre_spec)

    x    = new_err_wave
    y    = interp_func(new_err_wave)
    yerr = new_err
    
    return x, y, yerr












