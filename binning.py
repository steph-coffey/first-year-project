import numpy as np


def bin_err(wave, err, R):
    """
    Bins error bars to a resolution R:
    
    R = lamb / deltalamb
    
    Thus, the bin size changes for each wavelength.
    
    ================== Params ==================
    
    wave  : array or list of wavelengths (microns)
    err   : array or list of errors on transit spectrum (ppm)
    R     : spectral resolution (JWST NIRSPEC: ~100 (prism), ~1000, or ~2700)
    
    ============================================
    
    """
    
    bin_centers = []
    new_err     = []
    bin_widths  = [] # full bin widths
    
    # starting from the last bin and working backwards
    # to line up the rightmost bin edge with the last 
    # point in the given wavelength array, use this formula
    # to get the last bin center
    last_bin_center  = (2*R)/(2*R + 1) * wave[-1]
    last_dellamb     = last_bin_center / R
    last_bin_indices = np.where(wave >= last_bin_center - last_dellamb/2)
    
    # standard error of the weighted mean
    last_weights      = 1 / (err[last_bin_indices] ** 2)
    last_standard_err = np.sqrt(1 / np.sum(last_weights))
    
    # storing to lists
    bin_centers.append(last_bin_center)
    new_err.append(last_standard_err)
    bin_widths.append(last_dellamb)
    
    # calculate next bin center
    bin_center = last_bin_center * (2*R - 1)/(2*R + 1)
    dellamb    = bin_center / R
    
    # repeat calculation above until left edge of first bin falls outside wavelength array
    while bin_center - dellamb/2 >= wave[0]:
        
        dellamb      = bin_center / R
        bin_indices  = np.where(np.logical_and(wave >= bin_center - dellamb/2, wave < bin_center + dellamb/2))
        weights      = 1 / (err[bin_indices] ** 2)
        standard_err = np.sqrt(1 / np.sum(weights))
        
        bin_centers.append(bin_center)
        new_err.append(standard_err)
        bin_widths.append(dellamb)
        
        # calculate next bin center and dellamb
        bin_center = bin_center * (2*R - 1)/(2*R + 1)
        dellamb    = bin_center / R
        
    # for first bin, reset bin center to first element of wavelength array
    # bin_width is now just the rest of the wavelength array
    # take difference between previous bin edge and first point in the 
    # wavelength array to be the dellamb/2
    dellamb     = 2 * ((bin_center + dellamb/2) - wave[0])
    bin_center  = wave[0]
    bin_indices = np.where(np.logical_and(wave >= bin_center, wave < bin_center + dellamb/2))
    
    weights      = 1 / (err[bin_indices] ** 2)
    standard_err = np.sqrt(1 / np.sum(weights))

    bin_centers.append(bin_center)
    new_err.append(standard_err)
    bin_widths.append(dellamb)
    
    # flip lists to be order or increasing wavelength
    bin_centers.reverse()
    new_err.reverse()
    bin_widths.reverse()
    
    return bin_centers, new_err, bin_widths


def white_light_curve(wave, err):
    """
    Calculates the white light curve for a given spectrum. 
    This reduces the spectrum to a single point. Useful as a
    comparison between the error binning in Pandexo and in my code.
    
    ================== Params ==================
    
    wave  : array or list of wavelengths (microns)
    err   : array or list of errors on transit spectrum (ppm)
    
    ============================================
    """
    
    # loops over entire error array
    weights      = 1 / (err ** 2)
    standard_err = np.sqrt(1 / np.sum(weights))
    
    # bin center at the center of wavelength array
    bin_center = np.mean([wave[0],wave[-1]])
    bin_width  = wave[-1] - wave[0]
    
    return bin_center, standard_err, bin_width