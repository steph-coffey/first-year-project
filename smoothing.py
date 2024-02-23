import numpy as np

def adap_smooth(wave, spec, R):
    """
    Smooths to a fixed spectral resolution according to:
    
    R = lamb / deltalamb
    
    Thus, the bin size changes for each wavelength.
    
    ================== Params ==================
    
    wave   : array or list of wavelengths
    spec   : array or list of transit depths
    R      : spectral resolution (JWST NIRSPEC: ~100 (prism), ~1000, or ~2700)
    
    ============================================
    
    """
    
    # Finding center of first bin
    # using equation lamb = 2R*np.min(lambs) / (2R-1)
    lamb0 = (2*R*wave[0]) / (2*R - 1)
    dellamb = lamb0 / R
    
    # binning wavelengths w/n +- dellamb/2 of lamb0
    bin_indices = np.where(np.logical_and(wave >= lamb0 - dellamb/2, wave < lamb0 + dellamb/2))
    
    # taking mean in bin
    mean0 = spec[bin_indices].mean()
    
    # creating list to store the smooth spectrum
    smooth_spec = [mean0]
    
    # repeating everything above for every bin
    for i in range(1,len(spec)):
        lamb = wave[i]
        dellamb = lamb / R
        bin_indices = np.where(np.logical_and(wave >= lamb - dellamb/2, wave < lamb + dellamb/2))
        mean = spec[bin_indices].mean()

        smooth_spec.append(mean)
    
    return smooth_spec
