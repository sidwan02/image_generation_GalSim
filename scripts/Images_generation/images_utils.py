import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import galsim
import scipy

from cosmos_params import *

import photutils
from photutils.centroids import centroid_com

sys.path.insert(0,'../tools_generation/')
from tools_generation import utils

rng = galsim.BaseDeviate(None)


############ PARAMETER MEASUREMENTS

def get_data(gal, gal_image, psf_image, param_or_real='param'):
    '''
    Return redshift, moments_sigma, ellipticities and magnitude of the galaxy gal

    Parameters:
    ----------
    gal: galaxy from which you want to extract parameters
    gal_image: image of the galaxy gal
    psf_image: psf on the image gal_image, necessary to extract ellipticites and moments_sigma
    param_or_real: the measurement is for a parametric image or a real image
    '''
    shear_est = 'KSB' #'REGAUSS' for e (default) or 'KSB' for g
    res = galsim.hsm.EstimateShear(gal_image, psf_image, shear_est=shear_est, strict=True)
    if param_or_real == 'param':
        mag = gal.calculateMagnitude(filters['r'].withZeropoint(28.13))
    else:
        mag = np.nan
    if res.error_message == "":
        if shear_est != 'KSB':
            return [gal.SED.redshift, res.moments_sigma, res.corrected_e1, res.corrected_e2, mag]
        else:
            return [gal.SED.redshift, res.moments_sigma, res.corrected_g1, res.corrected_g2, mag]
    else:
        return [gal.SED.redshift, np.nan, np.nan, np.nan, mag]


############ SHIFTING GALAXIES

def shift_gal(gal, method='uniform', shift_x0=0., shift_y0=0., max_dx=0.1, min_r = fwhm_lsst/2., max_r = 2.):
    """
    Return galaxy shifted according to the chosen shifting method
    
    Parameters:
    ----------
    gal: galaxy to shift (GalSim object)
    method: method to use for shifting
    shift_x0: shift of centered/brightest galaxy to shift others according to its coordinates
    shift_y0: shift of centered/brightest galaxy to shift others according to its coordinates
    max_dx: dx maximum when using uniform shift
    min_r: minimum radius of annulus
    max_r: maximum radius of annulus
    """
    if method == 'noshift':
        shift_x = 0.
        shift_y = 0.
    elif method == 'uniform':
        shift_x = np.random.uniform(-max_dx,+max_dx)
        shift_y = np.random.uniform(-max_dx,+max_dx)
    elif method == 'annulus':
        r = np.sqrt(np.random.uniform(min_r**2, max_r**2))
        theta = np.random.uniform(0., 2*np.pi)
        shift_x = r * np.cos(theta)
        shift_y = r * np.sin(theta)
    elif method == 'uniform+betaprime':
        r = np.clip(scipy.stats.betaprime.rvs(*beta_prime_parameters), 0., 0.6)
        theta = np.random.uniform(0., 2*np.pi)
        shift_x = r * np.cos(theta) + np.random.uniform(-max_dx,+max_dx)
        shift_y = r * np.sin(theta) + np.random.uniform(-max_dx,+max_dx)
    else:
        raise ValueError
    shift_x += shift_x0
    shift_y += shift_y0
    return gal.shift((shift_x,shift_y)), (shift_x,shift_y)



########### PEAK DETECTION

def peak_detection(denormed_img, band, shifts, img_size, npeaks, nb_blended_gal, training_or_test, dist_cut):
    '''
    Return coordinates of the centroid of the closest galaxy from the center

    Parameters:
    ----------
    denormed_img: denormalized image of the blended galaxies
    band: filter in which the detection is done
    shifts: initial shifts used for the image generation
    img_size: size of the image
    npeaks: maximum number of peaks to return for photutils find_peaks function
    nb_blended_gal: number of blended galaxies in the image
    training_or_test: choice of training or test sample being generated
    dist_cut: cut in distance to check if detected galaxy is not too close from its neighbours
    '''
    gal = denormed_img
    df_temp = photutils.find_peaks(gal, threshold=5*np.sqrt(sky_level_pixel[band]), npeaks=npeaks, centroid_func=centroid_com)
    if df_temp is not None:
        df_temp['x_peak'] = (df_temp['x_centroid']-((img_size/2.)-0.5))*pixel_scale[band]
        df_temp['y_peak'] = (df_temp['y_centroid']-((img_size/2.)-0.5))*pixel_scale[band]
        df_temp.sort('peak_value', reverse=True)
        # Distances of true centers to brightest peak
        qq = [np.sqrt(float((shifts[j,0]-df_temp['x_peak'][0])**2+ (shifts[j,1]-df_temp['y_peak'][0])**2)) for j in range(nb_blended_gal)]
        idx_closest = np.argmin(qq)
        if nb_blended_gal>1:
            # Distance from peak galaxy to others
            qq_prime = [np.sqrt(float((shifts[idx_closest,0]-shifts[j,0])**2+ (shifts[idx_closest,1]-shifts[j,1])**2)) if j!=idx_closest else np.inf for j in range(nb_blended_gal)]
            idx_closest_to_peak_galaxy = np.argmin(qq_prime)
            if training_or_test != 'test':
                if not np.all(np.array(qq_prime) > dist_cut):
                    print('TRAINING CUT: closest is not central and others are too close')
                    return False
        else:
            idx_closest_to_peak_galaxy = np.nan
        return idx_closest, idx_closest_to_peak_galaxy, df_temp[0]['x_centroid'], df_temp[0]['y_centroid'], df_temp[0]['x_peak'], df_temp[0]['y_peak'], len(df_temp)
    else:
        return False

########## DRAWING OF IMAGE WITH GALSIM

def draw_images(galaxies_psf, band, img_size, filter_name,sky_level_pixel, real_or_param = 'param'):
    '''
    Return single galaxy noiseless images as well as the blended noisy one

    Parameters:
    ----------
    galaxies_psf: PSF to add on the image
    band: filter number in which the image is drawn
    img_size: size of the drawn image
    filter_name: name of the filter
    sky_level_pixel: sky level pixel for noise realization
    real_or_param: the galaxy generation use real image or parametric model
    '''
    # Create image in r bandpass filter to do the peak detection
    blend_img = galsim.ImageF(img_size, img_size, scale=pixel_scale[band])

    images = []
    
    for j, gal in enumerate(galaxies_psf):
        temp_img = galsim.ImageF(img_size, img_size, scale=pixel_scale[band])
        # Parametric image
        if real_or_param == 'param':
            gal.drawImage(filters[filter_name], image=temp_img)
        # Real image
        elif real_or_param == "real":
            gal.drawImage(image=temp_img)
        images.append(temp_img)
        blend_img += temp_img
    # add noise
    poissonian_noise = galsim.PoissonNoise(rng, sky_level=sky_level_pixel)
    blend_img.addNoise(poissonian_noise)

    return images, blend_img
