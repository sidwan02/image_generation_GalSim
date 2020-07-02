# Import packages

import numpy as np
import sys
import os
import galsim
from multiprocess import *
import pandas as pd
from tqdm import tqdm, trange

sys.path.insert(0,'../../scripts/tools_generation/')
import tools_generation
from tools_generation import utils

from images_generator import image_generator_sim, image_generator_real

# The script is used as, eg,
# >> python main_generation_cosmos.py centered/ simulation training isolated false 10 1000
# to produce 10 files in the training sample with 1000 images each of isolated galaxy centered on the image with no shift.
case = str(sys.argv[1]) # directory. Examples: centered/
gal_type = str(sys.argv[2]) #simulation or real/
training_or_test = str(sys.argv[3]) # this is a directory . Can be 'training', 'test' or 'validation'. 
isolated_or_blended = str(sys.argv[4]) #isolated blended
do_peak_detection = str(sys.argv[5]).lower() == 'true'
N_files = int(sys.argv[6]) # Nb of files to generate
N_per_file = int(sys.argv[7]) # Number of galaxies per file
assert training_or_test in ['training', 'validation', 'test']

# Fixed parameters:
max_try = 100 # maximum number of try before leaving the function (to avoir infinite loop)
mag_cut = 27.5 # cut in magnitude to select galaxies below this magnitude
max_stamp_size = 64 # Size of patch to generate
nmax_blend = (2,20) # Number of galaxies on an image if integer, or interval for sampling if tuple
center_brightest = True # Center the brightest galaxy (i.e. the galaxy with the lowest magnitude)
# If center_brightest = False : choose with method to use to shift the brightest
method_shift_brightest = 'noshift'
# And then you need to choose the method to shift the other galaxies as a function of the position of the brightest on the image
method_shift_others = 'uniform'
max_dx = 3.2 #in arcseconds, limit to use for uniform shifting: the center of the shifted galaxy will be shifted from the center or from the brightest galaxy from a random number between [-max_dx ; max_dx] arcsecond
max_r = 2. #in arcseconds, limit to use for annulus shifting: galaxy is shifted in an annulus around the center of the image or of the brightest galaxy which has for minimum radius fwhm_lsst/2 and for maximum radius max_r

# Method to shift centered galaxy
if isolated_or_blended == 'isolated':
    # where to save images and data
    save_dir = '/sps/lsst/users/barcelin/data/test/' + case + training_or_test
    # what to call those files
    root = 'galaxies_isolated_20191024_'
    # Maximum number of galaxies on the image. Here, "isolated" so only 1 galaxy.
    nmax_blend = 1
elif isolated_or_blended == 'blended':
    # where to save images and data
    save_dir = '/sps/lsst/users/barcelin/data/test/' + case + training_or_test
    # what to call those files
    root = 'galaxies_blended_20191024_'
    # Maximum number of galaxies on the image. Here, "isolated" so only 1 galaxy.
    nmax_blend = nmax_blend
else:
    raise NotImplementedError
# Path to the catalog
cosmos_cat_dir = '/sps/lsst/users/barcelin/COSMOS_25.2_training_sample'
# Loading the COSMOS catalog
cosmos_cat = galsim.COSMOSCatalog('real_galaxy_catalog_25.2.fits', dir=cosmos_cat_dir) 
# Select galaxies to keep for the test sample
if training_or_test == 'test':
    used_idx = np.arange(5000)
# Rest of the galaxies used for training and validation
else:
    used_idx = np.arange(5000,cosmos_cat.nobjects)

# keys for data objects
keys = []
for i in range (nmax_blend[1]):
    keys = keys + ['redshift_'+str(i), 'moment_sigma_'+str(i), 'e1_'+str(i), 'e2_'+str(i), 'mag_'+str(i)]
keys = keys + ['nb_blended_gal', 'SNR', 'SNR_peak', 'mag', 'mag_ir', 'closest_x', 'closest_y', 'closest_mag', 'closest_mag_ir',  'idx_closest_to_peak', 'n_peak_detected']

for icat in trange(N_files):
    # Run params
    root_i = root+str(icat)

    galaxies = []
    shifts = []

    #if training_or_test == 'test':
        # If test, create Pandas DataFrame to return properties of test galaxies
    # Here we save data for all datasets
    df = pd.DataFrame(index=np.arange(N_per_file), columns=keys)
    
    # Depending of type of galaxies you wand (simulation or real galaxies) use the correct generating function
    if gal_type == 'simulation':
        res = utils.apply_ntimes(image_generator_sim, N_per_file, (cosmos_cat_dir, training_or_test, isolated_or_blended, used_idx, nmax_blend, max_try, mag_cut, method_shift_brightest, method_shift_others, max_dx, max_r, do_peak_detection, center_brightest, max_stamp_size))
    elif gal_type == 'real':
        res = utils.apply_ntimes(image_generator_real, N_per_file, (cosmos_cat_dir, training_or_test, isolated_or_blended, used_idx, nmax_blend, max_try, mag_cut, method_shift_brightest, method_shift_others, max_dx, max_r, do_peak_detection, center_brightest, max_stamp_size))

    
    for i in trange(N_per_file):
        # Save data and shifts for all training, validation and test files
        gal_noiseless, blend_noisy, data, shift = res[i]
        assert set(data.keys()) == set(keys)
        df.loc[i] = [data[k] for k in keys]
        shifts.append(shift)
        galaxies.append((gal_noiseless, blend_noisy))

    # Save noisy blended images and denoised single central galaxy images
    np.save(os.path.join(save_dir, root_i+'_images.npy'), galaxies)
    # Save data and shifts
    df.to_csv(os.path.join(save_dir, root_i+'_data.csv'), index=False)
    np.save(os.path.join(save_dir, root_i+'_shifts.npy'), np.array(shifts))
    
    del galaxies, res, shifts, df