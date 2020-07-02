# Import packages
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import galsim
import multiprocessing
import time
from tqdm import tqdm, trange
import pathlib
from pathlib import Path

import plot

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


############### SNR COMPUTATION IN R-BAND FILTER

def SNR_peak(gal_noiseless, sky_background_pixel, band=6, snr_min=2):
    '''
    Return SNR computed with brightest peak value

    Parameters:
    ---------
    gal_noiseless: noiseless image of the isolated galaxy
    sky_background_pixel: sky background level per pixel
    band: filter number of r-band filter
    snr_min: minimum snr to keep the image
    '''
    # Make sure images have shape [nband, nx, ny] and sky_background_pixel has length nband
    assert len(sky_background_pixel) == gal_noiseless.shape[0]
    assert gal_noiseless.shape[1] == gal_noiseless.shape[2]
    
    snr = np.max(gal_noiseless[band])/sky_background_pixel[band]
    return (snr>snr_min), snr


def SNR(gal_noiseless, sky_background_pixel, band=6, snr_min=5):
    '''
    Return SNR

    Parameters:
    ---------
    gal_noiseless: noiseless image of the isolated galaxy
    sky_background_pixel: sky background level per pixel
    band: filter number of r-band filter
    snr_min: minimum snr to keep the image
    '''
    # Make sure images have shape [nband, nx, ny] and sky_background_pixel has length nband
    assert len(sky_background_pixel) == gal_noiseless.shape[0]
    assert gal_noiseless.shape[1] == gal_noiseless.shape[2]
    
    signal = gal_noiseless[band]
    variance = signal+sky_background_pixel[band] # for a Poisson process, variance=mean
    snr = np.sqrt(np.sum(signal**2/variance))
    return (snr>snr_min), snr




################## COMPUTE BLENDEDNESS 

def compute_blendedness_single(image1, image2):
    """
    Return blendedness computed with two images of single galaxy created with GalSim

    Parameters
    ----------
    image1, image2: GalSim images convolved with PSF
    """
    if isinstance(image1, galsim.image.Image):
        im1 = np.array(image1.array.data)
        im2 = np.array(image2.array.data)
    else:
        im1 = image1
        im2 = image2

    blnd = np.sum(im1*im2)/np.sqrt(np.sum(im1**2)*np.sum(im2**2))
    return blnd

def compute_blendedness_total(img_central, img_others):
    """
    Return blendedness computed with image of target galaxy and all neighbour galaxies

    Parameters
    ----------
    img_central, img_others: GalSim images convolved with PSF
    """
    if isinstance(img_central, galsim.image.Image):
        ic = np.array(img_central.array.data)
        io = np.array(img_others.array.data)
    else :
        ic = img_central
        io = img_others
    itot = ic + io
    blnd = 1. - np.sum(ic*ic)/np.sum(itot*ic)
    return blnd

def compute_blendedness_aperture(img_central, img_others, radius):
    """
    Return blendedness computed with image of target galaxy and all neighbour galaxies on a circle of radius "radius"

    Parameters
    ----------
    img_central, img_others: GalSim images convolved with PSF
    radius: radius of the circle used to compute blend rate
    """
    if isinstance(img_central, galsim.image.Image):
        ic = np.array(img_central.array.data)
        io = np.array(img_others.array.data)
    else :
        ic = img_central
        io = img_others
    h, w = ic.shape
    mask = plot.createCircularMask(h, w, center=None, radius=radius)
    flux_central = np.sum(ic*mask.astype(float))
    flux_others = np.sum(io*mask.astype(float))
    return flux_others / (flux_central+flux_others)



##############   MULTIPROCESSING    ############
def apply_ntimes(func, n, args, verbose=True, timeout=None):
    """
    Applies `n` times the function `func` on `args` (useful if, eg, `func` is partly random).
    Parameters
    ----------
    func : function
        func must be pickable, see https://docs.python.org/2/library/pickle.html#what-can-be-pickled-and-unpickled .
    n : int
    args : any
    timeout : int or float
        If given, the computation is cancelled if it hasn't returned a result before `timeout` seconds.
    Returns
    -------
    type
        Result of the computation of func(iter).
    """
    pool = multiprocessing.Pool()

    multiple_results = [pool.apply_async(func, args) for _ in range(n)]

    pool.close()
    
    return [res.get(timeout) for res in tqdm(multiple_results, desc='# castor.parallel.apply_ntimes', disable = True)]

