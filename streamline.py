import numpy as np
import matplotlib.pyplot as plt

# my functions
from smoothing import adap_smooth
from run_pandexo import run_pandexo
from binning import bin_errs, calc_bin_widths, white_light_curve


# path to transit spectra
spectra_path = '/Users/coffey/Downloads/kipping/Exo_Transmit/Spectra'

def CSRE(Rp, Rs, Rref, planet_name, star_name, flat = True, res = 100):
    """
    Calculates the chromatic stellar radii effect (CSRE).
    
    ================== Params ==================
    
    Rp          : radius of the transiting planet [m]
    Rs          : radius of the host star [m]
    Rref        : radius of the reference star used in ExoTransmit runs [m]
    planet_name : name of planet (string)
    star_name   : name of star (string)
    flat        : whether or not to include planet's absorption features (True assumes a flat transit spectrum)
    res         : resolution to smooth final spectrum to, we recommend 100 to match JWST NIRSpec's PRISM grating
    
    ================== Returns =================
    
    smooth_csre_spec  : transit spectrum including the CSRE, smoothed to chosen resolution.
    
    """
    
    star_spectrum = np.loadtxt(f'{spectra_path}/transmission_{star_name}.dat', skiprows = 2).T  # ExoTransmit 0.3 - 30 um
    
    # slicing wavelength range to JWST NIRSpec (0.6 - 5.3 um)
    wave_index = slice(np.where(star_spectrum[0] < 0.6e-6)[0][-1], np.where(star_spectrum[0] > 5.3e-6)[0][0])
    wave       = star_spectrum[0][wave_index] * 1e6  # microns
    star_spec  = star_spectrum[1][wave_index] / 100  # convert from %
    
    # Multiplying out the Sun's radius
    star_spec_noRref = np.sqrt(star_spec) * Rref
    
    ### Default mode ###
    if flat == True:
        # combining the spectra for our new (Rp/Rs)^2
        csre_spec  = (Rp / star_spec_noRref)**2 * 1e6  # ppm 

        # smoothing to JWST resolution (R ~ 100 for PRISM)
        smooth_csre_spec = adap_smooth(wave, csre_spec, R = res)
    
    
    ### include planet absorption ###
    if flat == False:
        planet_spectrum = np.loadtxt(f'{spectra_path}/transmission_{planet_name}.dat', skiprows = 2).T # ExoTransmit 0.3 - 30 um
        planet_spec     = planet_spectrum[1][wave_index] / 100  # convert from %
        
        # Multiplying out the Sun's radius
        planet_spec_noRref = np.sqrt(planet_spec) * Rref
        
        # combining the spectra for our new (Rp/Rs)^2
        csre_spec  = (planet_spec_noRref / star_spec_noRref)**2 * 1e6  # ppm 
    
        # smoothing to JWST resolution (R ~ 100 for PRISM is default)
        smooth_csre_spec = adap_smooth(wave, csre_spec, R = res)
    
    return wave, smooth_csre_spec

    
    
def errors(wave, smooth_csre_spec):
    """
    Places errorbars on the final transit spectrum.
    """
    
    
    err_wave, err = run_pandexo('trappist', 11.354, 'J', star_temp = 2566, star_metal = -1.4, star_logg = 5.276,
                       trans_dur = 0.6, star_radius = 0.119, planet_radius = 0.0996)
    
    R = 25
    new_errs, new_wave, dellambs = bin_errs(err_wave, err, R)
    
    # shifting the wavelength values for the error bars to 
    # closest values for the wavelength array corresponding to the spectrum
    # so that I can plot the error bars onto the spectrum itself

    newer_wave = []
    indices = []
    for i in range(len(new_wave)):
        index = (np.abs(wave - new_wave[i])).argmin()
        indices.append(index)
        newer_wave.append(wave[index])
        
    
    x = newer_wave
    y = np.take(smooth_csre_spec, indices)
    yerr = new_errs
    
    return x, y, yerr


def plotting(wave, smooth_csre_spec, Rp, Rs, err = False):
    """
    Makes some nice plots.
    """
    
    # fixed transit depth for comparison
    # includes no absorption effects from planet OR star!
    fixed_depth = (Rp / Rs)**2 * 1e6 # ppm
    
    fig, ax = plt.subplots(figsize = (9,4))
        
    ax.plot(wave, smooth_csre_spec, c = 'seagreen', label = 'CSRE')
    ax.axhline(y = fixed_depth, c = 'black', ls = 'dotted', label = 'fixed transit depth')
    
    if err == True:
        x,y,yerr = errors(wave, smooth_csre_spec)
        ax.errorbar(x, y, yerr, None, 'o', ms = 2, c = 'gray', label = '1 transit', capsize = 2.5, alpha = 0.5)
        ax.errorbar(x, y, yerr / np.sqrt(10), None, 'o', ms = 2, c = 'gray', label = '10 transits', capsize = 2.5, alpha = 0.75)
        ax.errorbar(x, y, yerr / np.sqrt(100), None, 'o', ms = 2, c = 'gray', label = '100 transits', capsize = 2.5, alpha = 1)

        
    # Labels
    fig.suptitle('Chromatic Stellar Radii Effect')
    ax.set_xlabel('Wavelength (um)')
    ax.set_ylabel('Transit depth (ppm)')
    ax.legend()

    fig.show()

    return fig, ax












