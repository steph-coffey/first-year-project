import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import os

import pandexo.engine.justplotit as jpi
import pickle as pk

def run_pandexo(params):
    """
    Runs Pandexo and returns wavelength axis (microns) and error bars (ppm).
    For our purposes, we don't need Pandexo's transit spectrum model.
    Comments on specifics in my pandexo notebook and/or their documentation.
    
    ==================================== Params ====================================
    
    sys_name      : name of transiting planet system
    mag           : magnitude of system, at reference wavelength
    ref_wave      : J mag = 1.25, H = 1.6, K = 2.22, etc.
    star_temp     : in K
    star_metal    : as log Fe/H
    star_logg     : log surface gravity cgs
    trans_dur     : in hours
    star_radius   : in Solar radii
    planet_radius : in Jupiter radii
    
    ================================================================================

    """
    
    # loading PandExo params
    sys_name, mag, ref_wave, star_temp, star_metal, \
    star_logg, trans_dur, star_radius, planet_radius = params.pandexo_params.values()
    
    ref_wave_dict = {'J' : 1.25, 'H' : 1.6, 'K' : 2.22}

    ### Loading blank exo dictionary
    exo_dict = jdi.load_exo_dict()

    ### Observational properties
    exo_dict['observation']['sat_level']     = 80
    exo_dict['observation']['sat_unit']      = '%'
    exo_dict['observation']['noccultations'] = 1
    exo_dict['observation']['R']             = None
    exo_dict['observation']['baseline_unit'] = 'total'
    exo_dict['observation']['baseline']      = (5 * trans_dur) *60.0*60.0
    exo_dict['observation']['noise_floor']   = 14

    ### Star properties
    exo_dict['star']['type']     = 'phoenix'
    exo_dict['star']['mag']      = mag
    exo_dict['star']['ref_wave'] = ref_wave_dict[ref_wave]
    exo_dict['star']['temp']     = star_temp
    exo_dict['star']['metal']    = star_metal
    exo_dict['star']['logg']     = star_logg
    exo_dict['star']['radius']   = star_radius
    exo_dict['star']['r_unit']   = 'R_sun'

    ### Planet properties
    exo_dict['planet']['type']             = 'constant'
    exo_dict['planet']['transit_duration'] = trans_dur *60.0*60.0
    exo_dict['planet']['radius']           = planet_radius
    exo_dict['planet']['r_unit']           = 'R_jup'
    exo_dict['planet']['td_unit']          = 's'
    exo_dict['planet']['w_unit']           = 'um'
    exo_dict['planet']['f_unit']           = 'rp^2/r*^2'
    
    ### The actual run
    result_dict = jdi.run_pandexo(exo_dict, ['NIRSpec Prism'], output_file = f'{sys_name}_pandexo.p')

    ### Pulling out wavelength and error values (copied from source)
    wave = result_dict['FinalSpectrum']['wave']
    err = result_dict['FinalSpectrum']['error_w_floor']*1e6
    wave = wave[~np.isnan(err)]
    err = err[~np.isnan(err)]

    return wave, err



