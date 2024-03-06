import numpy as np

def calc_bin_widths(wave, R):
    """
    Calculates deltalambda of a wavelength array, given the resolution:
    
    R = lamb / deltalamb
    
    This is used to place x errorbars on the points and ensure they don't overlap.
    
    ================== Params ==================
    
    wave   : array or list of wavelengths (microns)
    R      : spectral resolution (JWST NIRSPEC: ~100 (prism), ~1000, or ~2700)
    
    ============================================
    """
    
    bin_widths = []
    
    for lamb in wave:
        dellamb = lamb / R
        bin_widths.append(dellamb)
        
    return bin_widths

def bin_errs(wave, errs, R):
    """
    Bins error bars to a resolution R:
    
    R = lamb / deltalamb
    
    Thus, the bin size changes for each wavelength.
    
    ================== Params ==================
    
    wave   : array or list of wavelengths (microns)
    errs   : array or list of errors on transit spectrum (ppm)
    R      : spectral resolution (JWST NIRSPEC: ~100 (prism), ~1000, or ~2700)
    
    ============================================
    
    """
    
    new_errs    = []
    bin_centers = []
    bin_widths  = [] # full bin widths
    
    # treating first bin a bit different
    dellamb0     = wave[0] / R # bin center is 0th element
    bin_indices0 = np.where(wave <= wave[0] + dellamb0/2)
    
    # standard error of the weighted mean
    weights0      = 1 / (errs[bin_indices0] ** 2)
    standard_err0 = np.sqrt(1 / np.sum(weights0))
    
    # storing to lists
    new_errs.append(standard_err0)
    bin_centers.append(wave[0])
    bin_widths.append(dellamb0)
    
    # calculate next bin center
    bin_center = (2*R*wave[bin_indices0[0][-1] + 1]) / (2*R - 1)
    """print('bin min ', 'na')
    print('bin center ', wave[0])
    print('bin max ', wave[wave <= wave[0] + dellamb0/2][-1])
    print('next bin min ', wave[bin_indices0[0][-1] + 1])
    print("---------------------------------------")"""
    # repeat calculation above until last bin
    while bin_center <= wave[-1]:
        
        dellamb      = bin_center / R
        bin_indices  = np.where(np.logical_and(wave >= bin_center - dellamb/2, wave <= bin_center + dellamb/2))
        weights      = 1 / (errs[bin_indices] ** 2)
        standard_err = np.sqrt(1 / np.sum(weights))
        """print('bin min ', wave[wave >= bin_center - dellamb/2][0])
        print('bin center ', bin_center)
        print('bin max ', wave[wave <= bin_center + dellamb/2][-1])
        print('next bin min ', wave[bin_indices[0][-1] + 1])
        print("---------------------------------------")"""
        new_errs.append(standard_err)
        bin_centers.append(bin_center)
        bin_widths.append(dellamb)
        #print(bin_indices)
        """print('bin center: ', bin_center)
        print('bin width: ', dellamb)
        print('bin min: ', bin_center - dellamb/2)
        print('bin max: ', bin_center + dellamb/2)
        print('--------------------------------------')"""
        
        # calculate next bin center
        bin_center = (2*R*wave[bin_indices[0][-1] + 1]) / (2*R - 1)
        
    # for last bin, reset bin center to last element of wavelength array
    # bin_width is now just the rest of the wavelength array
    bin_center   = wave[-1]
    bin_indices  = np.where(wave >= bin_centers[-1] + bin_widths[-1]/2)
    dellamb      = 2 * (bin_center - wave[bin_indices[0][0]])
    weights      = 1 / (errs[bin_indices] ** 2)
    standard_err = np.sqrt(1 / np.sum(weights))

    new_errs.append(standard_err)
    bin_centers.append(bin_center)
    bin_widths.append(dellamb)

    
    return new_errs, bin_centers, bin_widths
