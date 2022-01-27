# Import packages

import numpy as np
import sys
import os
import galsim

from cosmos_params import *

import photutils
from photutils.centroids import centroid_com

import utils
from images_utils import get_fit_data, get_data, shift_gal, peak_detection, draw_images

rng = galsim.BaseDeviate(None)

# IMAGES NUMPY ARRAYS GENERATION
# CASE OF PARAMETRIC IMAGES - SIMULATION


def image_generator_sim(cosmos_cat_dir,
                        training_or_test,
                        isolated_or_blended,
                        used_idx=None,
                        nmax_blend=4,
                        max_try=3,
                        mag_cut=28.,
                        method_first_shift='noshift',
                        method_others_shift='uniform',
                        max_dx=3.2,
                        max_r=2.,
                        do_peak_detection=True,
                        center_brightest=True,
                        max_stamp_size=64,
                        psf_lsst_fixed=True):
    """
    Return numpy arrays: noiseless and noisy image of single galaxy and of blended galaxies as well as the pandaframe including data about the image and the shifts in the test sample generation configuration

    Parameters:
    ----------
    cosmos_cat_dir: COSMOS catalog directory
    training_or_test: choice for generating a training or testing dataset
    isolated_or_blended: choice for generation of samples of isolated galaxy images or blended galaxies images
    used_idx: indexes to use in the catalog (to use different parts of the catalog for training/validation/test)
    nmax_blend: maximum number of galaxies in a blended galaxies image
    max_try: maximum number of try before leaving the function (to avoir infinite loop)
    mag_cut: cut in magnitude to select function below this magnitude
    method_first_shift: chosen method for shifting the centered galaxy
    do_peak_detection: boolean to do the peak detection
    """
    # Define PSF
    PSF_lsst, fwhm_lsst = psf_lsst(psf_lsst_fixed=psf_lsst_fixed)
    PSF = [PSF_euclid_nir]*3 + [PSF_euclid_vis] + [PSF_lsst]*6
    print("PSF[6]: ", PSF[6])
    # Import the COSMOS catalog
    cosmos_cat = galsim.COSMOSCatalog(
        'real_galaxy_catalog_25.2.fits', dir=cosmos_cat_dir)
    counter = 0
    np.random.seed()  # important for multiprocessing !

    assert training_or_test in ['training', 'validation', 'test']
    assert isolated_or_blended in ['blended', 'isolated']

    while counter < max_try:
        try:
            ud = galsim.UniformDeviate()

            if np.shape(nmax_blend) == ():
                nb_blended_gal = nmax_blend
            else:
                nb_blended_gal = np.random.randint(
                    nmax_blend[0], nmax_blend[1])
                nmax_blend = nmax_blend[1]
            data = {}
            galaxies = []
            mag = []
            mag_ir = []
            j = 0
            while j < nb_blended_gal:
                # Chose the part of the catalog used for generation
                if used_idx is not None:
                    idx = np.random.choice(used_idx)
                else:
                    idx = np.random.randint(cosmos_cat.nobject)
                # Generate galaxy
                gal = cosmos_cat.makeGalaxy(
                    idx, gal_type='parametric', chromatic=True, noise_pad_size=0)
                # Get data from fit (parametric model)
                data['e1_fit_'+str(j)], data['e2_fit_'+str(j)], data['weight_fit_' +
                                                                     str(j)] = get_fit_data(cosmos_cat_dir, idx)
                # Compute the magnitude of the galaxy
                _mag_temp = gal.calculateMagnitude(
                    filters['r'].withZeropoint(28.13))
                # Magnitude cut
                if _mag_temp < mag_cut:
                    gal = gal.rotate(ud() * 360. * galsim.degrees)
                    galaxies.append(gal)
                    mag.append(_mag_temp)
                    mag_ir.append(gal.calculateMagnitude(
                        filters['H'].withZeropoint(24.92-22.35*coeff_noise_h)))
                    j += 1
            for i in range(nmax_blend-nb_blended_gal):
                data['e1_fit_'+str(nb_blended_gal+i)], data['e2_fit_'+str(nb_blended_gal+i)
                                                            ], data['weight_fit_'+str(nb_blended_gal+i)] = [np.nan, np.nan, np.nan]

            # Compute ellipticities and magnitude for galaxies in r band before the shifting.
            psf_image = PSF[6].drawImage(
                nx=max_stamp_size, ny=max_stamp_size, scale=pixel_scale[6])
            images = []
            galaxies_psf = [galsim.Convolve(
                [gal*coeff_exp[6], PSF[6]]) for gal in galaxies]
            for j, gal in enumerate(galaxies_psf):
                temp_img = galsim.ImageF(
                    max_stamp_size, max_stamp_size, scale=pixel_scale[6])

                gal.drawImage(filters['r'], image=temp_img)
                images.append(temp_img)

            for z in range(nb_blended_gal):
                data['redshift_'+str(z)], data['moment_sigma_'+str(z)], data['e1_ksb_'+str(
                    z)], data['e2_ksb_'+str(z)], data['mag_'+str(z)] = get_data(galaxies[z], images[z], psf_image)
            if nb_blended_gal < nmax_blend:
                for z in range(nb_blended_gal, nmax_blend):
                    data['redshift_'+str(z)], data['moment_sigma_'+str(z)], data['e1_ksb_'+str(
                        z)], data['e2_ksb_'+str(z)], data['mag_'+str(z)] = 10., 10., 10., 10., 10.

            # Optionally, find the brightest and put it first in the list
            if center_brightest:
                _idx = np.argmin(mag)
                galaxies.insert(0, galaxies.pop(_idx))
                mag.insert(0, mag.pop(_idx))
                mag_ir.insert(0, mag_ir.pop(_idx))

            # Shifts galaxies
            shift = np.zeros((nmax_blend, 2))
            if center_brightest == False:
                # Shift the lowest magnitude galaxy
                galaxies[0], shift[0] = shift_gal(
                    galaxies[0], method=method_first_shift, max_dx=max_dx, max_r=max_r)
            # Shift all the other galaxies
            for j, gal in enumerate(galaxies[1:]):
                galaxies[j+1], shift[j+1] = shift_gal(
                    gal, method=method_others_shift, max_dx=max_dx, max_r=max_r)

            # Compute distances of the neighbour galaxies to the lowest magnitude galaxy
            if nb_blended_gal > 1:
                distances = [shift[j][0]**2+shift[j][1] **
                             2 for j in range(1, nb_blended_gal)]
                idx_closest_to_peak_galaxy = np.argmin(distances)+1
            else:
                idx_closest_to_peak_galaxy = 0

            if training_or_test == 'test':
                galaxy_noiseless = np.zeros(
                    (nmax_blend, 10, max_stamp_size, max_stamp_size))
            else:
                galaxy_noiseless = np.zeros(
                    (10, max_stamp_size, max_stamp_size))
            blend_noisy = np.zeros((10, max_stamp_size, max_stamp_size))

            # Realize peak detection in r-band filter if asked
            if do_peak_detection:
                band = 6
                galaxies_psf = [galsim.Convolve(
                    [gal*coeff_exp[band], PSF[band]]) for gal in galaxies]

                images, blend_img = draw_images(
                    galaxies_psf, band, max_stamp_size*2, 'r', sky_level_pixel[band])
                blend_noisy_temp = blend_img.array.data
                peak_detection_output = peak_detection(
                    blend_noisy_temp, band, shift, max_stamp_size*2, 4, nb_blended_gal, training_or_test, dist_cut=0.65/2.)
                if not peak_detection_output:
                    print('No peak detected')
                    raise RuntimeError
                else:
                    idx_closest_to_peak, idx_closest_to_peak_galaxy, center_pix_x, center_pix_y, center_arc_x, center_arc_y, n_peak = peak_detection_output

                # Modify galaxies and shift accordingly
                galaxies = [gal.shift(-center_arc_x, -center_arc_y)
                            for gal in galaxies]
                shift[:nb_blended_gal] -= np.array(
                    [center_arc_x, center_arc_y])

            # Now draw image in all filters
            for i, filter_name in enumerate(filter_names_all):
                galaxies_psf = [galsim.Convolve(
                    [gal*coeff_exp[i], PSF[i]]) for gal in galaxies]
                images, blend_img = draw_images(
                    galaxies_psf, i, max_stamp_size, filter_name, sky_level_pixel[i])
                if isolated_or_blended == 'isolated' or not do_peak_detection:
                    idx_closest_to_peak = 0
                    n_peak = 1

                if training_or_test == 'test':
                    galaxy_noiseless[0][i] = images[idx_closest_to_peak].array.data
                    if isolated_or_blended == 'blended':
                        for m in range(1, nb_blended_gal):
                            if m <= idx_closest_to_peak:
                                galaxy_noiseless[m][i] = images[m-1].array.data
                            elif m > idx_closest_to_peak:
                                galaxy_noiseless[m][i] = images[m].array.data
                else:
                    galaxy_noiseless[i] = images[idx_closest_to_peak].array.data
                blend_noisy[i] = blend_img.array.data
            break

        except RuntimeError as e:
            print(e)

    # For testing, return unormalized images and data
    data['fwhm_lsst'] = fwhm_lsst
    data['nb_blended_gal'] = nb_blended_gal
    data['mag'] = mag[0]
    data['mag_ir'] = mag_ir[0]
    if nb_blended_gal > 1:
        data['closest_mag'] = mag[idx_closest_to_peak_galaxy]
        data['closest_mag_ir'] = mag_ir[idx_closest_to_peak_galaxy]
        data['closest_x'] = shift[idx_closest_to_peak_galaxy][0]
        data['closest_y'] = shift[idx_closest_to_peak_galaxy][1]
    else:
        data['closest_mag'] = np.nan
        data['closest_mag_ir'] = np.nan
        data['closest_x'] = np.nan
        data['closest_y'] = np.nan
    data['idx_closest_to_peak'] = idx_closest_to_peak
    data['n_peak_detected'] = n_peak
    data['SNR'] = utils.SNR(galaxy_noiseless, sky_level_pixel, band=6)[1]
    data['SNR_peak'] = utils.SNR_peak(
        galaxy_noiseless, sky_level_pixel, band=6)[1]
    return galaxy_noiseless, blend_noisy, data, shift


# CASE OF REAL IMAGES
def image_generator_real(cosmos_cat_dir,
                         training_or_test,
                         isolated_or_blended,
                         used_idx=None,
                         nmax_blend=4,
                         max_try=3,
                         mag_cut=28.,
                         method_first_shift='noshift',
                         method_others_shift='uniform',
                         max_dx=3.2,
                         max_r=2.,
                         do_peak_detection=True,
                         center_brightest=True,
                         max_stamp_size=64,
                         psf_lsst_fixed=True):
    """
    Return numpy arrays: noiseless and noisy image of single galaxy and of blended galaxies as well as the pandaframe including data about the image and the shifts in the test sample generation configuration

    Parameters:
    ----------
    cosmos_cat_dir: COSMOS catalog directory
    training_or_test: choice for generating a training or testing dataset
    isolated_or_blended: choice for generation of samples of isolated galaxy images or blended galaxies images
    used_idx: indexes to use in the catalog (to use different parts of the catalog for training/validation/test)
    nmax_blend: maximum number of galaxies in a blended galaxies image
    max_try: maximum number of try before leaving the function (to avoir infinite loop)
    mag_cut: cut in magnitude to select function below this magnitude
    method_first_shift: chosen method for shifting the centered galaxy
    do_peak_detection: boolean to do the peak detection
    """
    # Define PSF
    PSF_lsst = psf_lsst(psf_lsst_fixed=False)
    PSF = [PSF_euclid_nir]*3 + [PSF_euclid_vis] + [PSF_lsst]*6
    # print("PSF[6]: ", PSF[6])
    # Import the COSMOS catalog
    cosmos_cat = galsim.COSMOSCatalog(
        'real_galaxy_catalog_25.2.fits', dir=cosmos_cat_dir)
    counter = 0
    np.random.seed()  # important for multiprocessing !

    assert training_or_test in ['training', 'validation', 'test']
    assert isolated_or_blended in ['blended', 'isolated']

    while counter < max_try:
        try:
            ud = galsim.UniformDeviate()
            real_gal_list = []

            if np.shape(nmax_blend) == ():
                nb_blended_gal = nmax_blend
            else:
                nb_blended_gal = np.random.randint(
                    nmax_blend[0], nmax_blend[1])
                nmax_blend = nmax_blend[1]
            data = {}
            galaxies = []
            mag = []
            mag_ir = []
            j = 0
            while j < nb_blended_gal:
                # Chose the part of the catalog used for generation
                if used_idx is not None:
                    idx = np.random.choice(used_idx)
                else:
                    idx = np.random.randint(cosmos_cat.nobject)
                # Generate galaxy
                gal = cosmos_cat.makeGalaxy(
                    idx, gal_type='parametric', chromatic=True, noise_pad_size=0)
                # Compute the magnitude of the galaxy
                _mag_temp = gal.calculateMagnitude(
                    filters['r'].withZeropoint(28.13))
                # Magnitude cut
                if _mag_temp < mag_cut:
                    gal = gal.rotate(ud() * 360. * galsim.degrees)
                    galaxies.append(gal)
                    mag.append(_mag_temp)
                    mag_ir.append(gal.calculateMagnitude(
                        filters['H'].withZeropoint(24.92-22.35*coeff_noise_h)))
                    j += 1

                # Take the real galaxy image only if parametric galaxy is actually created
                if len(galaxies) == (len(real_gal_list)+1):
                    bp_file = os.path.join(os.path.join(os.path.dirname(
                        os.path.realpath(__file__)), '../data/share_galsim/'), 'wfc_F814W.dat.gz')
                    bandpass = galsim.Bandpass(
                        bp_file, wave_type='ang').thin().withZeropoint(25.94)
                    real_gal = cosmos_cat.makeGalaxy(idx, gal_type='real',
                                                     noise_pad_size=max_stamp_size*pixel_scale_lsst)
                    real_gal_list.append(real_gal)

            # Compute ellipticities and magnitude for galaxies in r band before the shifting.

            # for some reason, unlike in image_generator_sim, PSF is produced as a list of tuples for real image_generator_real

            # # ADDED ===
            # psf_6 = PSF[6]
            # # =========

            for x in PSF:
                print(type(x))

            PSF = list(map(lambda x: x[0] if type(x) == tuple else x, PSF))
            print("PSF: ", PSF)

            psf_image = PSF[6].drawImage(
                nx=max_stamp_size, ny=max_stamp_size, scale=pixel_scale[6])
            images = []
            galaxies_psf = [galsim.Convolve(
                [real_gal*coeff_exp[6], PSF[6]]) for real_gal in real_gal_list]
            for j, gal in enumerate(galaxies_psf):
                temp_img = galsim.ImageF(
                    max_stamp_size, max_stamp_size, scale=pixel_scale[6])

                gal.drawImage(image=temp_img)  # filters['r'],
                images.append(temp_img)

            for z in range(nb_blended_gal):
                res, data['mag_'+str(z)] = get_data(real_gal_list[z],
                                                    images[z], psf_image, param_or_real='real'), mag[z]
                data['redshift_'+str(z)] = res[0]
                data['moment_sigma_'+str(z)] = res[1]
                data['e1_'+str(z)] = res[2]
                data['e2_'+str(z)] = res[3]
            if nb_blended_gal < nmax_blend:
                for z in range(nb_blended_gal, nmax_blend):
                    data['redshift_'+str(z)], data['moment_sigma_'+str(z)], data['e1_ksb_'+str(
                        z)], data['e2_ksb_'+str(z)], data['mag_'+str(z)] = 10., 10., 10., 10., 10.

            # Optionally, find the brightest and put it first in the list
            if center_brightest:
                _idx = np.argmin(mag)
                galaxies.insert(0, galaxies.pop(_idx))
                real_gal_list.insert(0, real_gal_list.pop(_idx))
                mag.insert(0, mag.pop(_idx))
                mag_ir.insert(0, mag_ir.pop(_idx))

            # Shifts galaxies
            shift = np.zeros((nmax_blend, 2))
            if center_brightest == False:
                # Shift the lowest magnitude galaxy
                real_gal_list[0], shift[0] = shift_gal(
                    real_gal_list[0], method=method_first_shift, max_dx=max_dx, max_r=max_r)
            # Shift all the other galaxies
            for j, gal in enumerate(real_gal_list[1:]):
                real_gal_list[j+1], shift[j+1] = shift_gal(
                    gal, method=method_others_shift, max_dx=max_dx, max_r=max_r)

            # Compute distances of the neighbour galaxies to the lowest magnitude galaxy
            if nb_blended_gal > 1:
                distances = [shift[j][0]**2+shift[j][1] **
                             2 for j in range(1, nb_blended_gal)]
                idx_closest_to_peak_galaxy = np.argmin(distances)+1
            else:
                idx_closest_to_peak_galaxy = 0

            if training_or_test == 'test':
                galaxy_noiseless = np.zeros(
                    (nmax_blend, 10, max_stamp_size, max_stamp_size))
                galaxy_noiseless_real = np.zeros(
                    (nmax_blend, 10, max_stamp_size, max_stamp_size))
            else:
                galaxy_noiseless = np.zeros(
                    (10, max_stamp_size, max_stamp_size))
                galaxy_noiseless_real = np.zeros(
                    (10, max_stamp_size, max_stamp_size))
            blend_noisy = np.zeros((10, max_stamp_size, max_stamp_size))
            blend_noisy_real = np.zeros((10, max_stamp_size, max_stamp_size))

            # Realize peak detection in r-band filter if asked
            if do_peak_detection:
                band = 6
                galaxies_psf = [galsim.Convolve(
                    [real_gal*coeff_exp[band], PSF[band]]) for real_gal in real_gal_list]

                images, blend_img = draw_images(
                    galaxies_psf, band, max_stamp_size*2, 'r', sky_level_pixel[band])
                blend_noisy_temp = blend_img.array.data
                peak_detection_output = peak_detection(
                    blend_noisy_temp, band, shift, max_stamp_size*2, 4, nb_blended_gal, training_or_test, dist_cut=0.65/2.)
                if not peak_detection_output:
                    print('No peak detected')
                    raise RuntimeError
                else:
                    idx_closest_to_peak, idx_closest_to_peak_galaxy, center_pix_x, center_pix_y, center_arc_x, center_arc_y, n_peak = peak_detection_output

                # Modify galaxies and shift accordingly
                galaxies = [gal.shift(-center_arc_x, -center_arc_y)
                            for gal in galaxies]
                shift[:nb_blended_gal] -= np.array(
                    [center_arc_x, center_arc_y])

            # Draw real images
            galaxies_real_psf = [galsim.Convolve(
                [real_gal*coeff_exp[6], PSF_lsst]) for real_gal in real_gal_list]
            images_real, _ = draw_images(
                galaxies_real_psf, 6, max_stamp_size, 'r', sky_level_pixel[6], real_or_param='real')

            # Now draw image in all bands
            for i, filter_name in enumerate(filter_names_all):
                galaxies_psf = [galsim.Convolve(
                    [gal*coeff_exp[i], PSF[i]]) for gal in galaxies]
                images, blend_img = draw_images(
                    galaxies_psf, i, max_stamp_size, filter_name, sky_level_pixel[i])
                if isolated_or_blended == 'isolated' or not do_peak_detection:
                    idx_closest_to_peak = 0
                    n_peak = 1

                if training_or_test == 'test':
                    if isolated_or_blended == 'blended':
                        for m in range(1, nb_blended_gal):
                            if m <= idx_closest_to_peak:
                                galaxy_noiseless[m][i] = images[m-1].array.data
                            elif m > idx_closest_to_peak:
                                galaxy_noiseless[m][i] = images[m].array.data
                else:
                    galaxy_noiseless[i] = images[idx_closest_to_peak].array.data
                blend_noisy[i] = blend_img.array.data

                # Rescale real images by flux
                images_real_array = np.zeros(
                    (len(images_real), max_stamp_size, max_stamp_size))
                for jj, image_real in enumerate(images_real):
                    img_temp = images[jj]
                    image_real -= np.min(image_real.array)
                    images_real_array[jj] = image_real.array * \
                        np.sum(img_temp.array)/np.sum(image_real.array)

                # real galaxies
                if training_or_test == 'test':
                    galaxy_noiseless_real[0][i] = images_real_array[idx_closest_to_peak].data
                    if isolated_or_blended == 'blended':
                        for m in range(1, nb_blended_gal):
                            if m <= idx_closest_to_peak:
                                galaxy_noiseless_real[m][i] = images_real_array[m-1].data
                            elif m > idx_closest_to_peak:
                                galaxy_noiseless_real[m][i] = images_real_array[m].data
                else:
                    galaxy_noiseless_real[i] = images_real_array[idx_closest_to_peak].data
                for image_real_array in images_real_array:
                    blend_noisy_real[i] += image_real_array

                # Add noise
                blend_noisy_real_temp = galsim.Image(
                    blend_noisy_real[i], dtype=np.float64)
                poissonian_noise = galsim.PoissonNoise(
                    rng, sky_level=sky_level_pixel[i])
                blend_noisy_real_temp.addNoise(poissonian_noise)
                blend_noisy_real[i] = blend_noisy_real_temp.array.data

            break

        except RuntimeError as e:
            print(e)

    data['fwhm_lsst'] = fwhm_lsst
    data['nb_blended_gal'] = nb_blended_gal
    data['mag'] = mag[0]
    data['mag_ir'] = mag_ir[0]
    if nb_blended_gal > 1:
        data['closest_mag'] = mag[idx_closest_to_peak_galaxy]
        data['closest_mag_ir'] = mag_ir[idx_closest_to_peak_galaxy]
        data['closest_x'] = shift[idx_closest_to_peak_galaxy][0]
        data['closest_y'] = shift[idx_closest_to_peak_galaxy][1]
    else:
        data['closest_mag'] = np.nan
        data['closest_mag_ir'] = np.nan
        data['closest_x'] = np.nan
        data['closest_y'] = np.nan
    data['idx_closest_to_peak'] = idx_closest_to_peak
    data['n_peak_detected'] = n_peak
    data['SNR'] = utils.SNR(galaxy_noiseless, sky_level_pixel, band=6)[1]
    data['SNR_peak'] = utils.SNR_peak(
        galaxy_noiseless, sky_level_pixel, band=6)[1]
    return galaxy_noiseless_real, blend_noisy_real, data, shift
