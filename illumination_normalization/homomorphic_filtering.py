#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
A simple homomorphic filtering of an image followed by an histogram
equalization
"""

from os import path
import numpy as np
import skimage
from skimage.exposure import rescale_intensity, equalize
from skimage.util.shape import view_as_windows
import skimage.io as io

from generate_log import generate_log

EPSILON = 1e-5


def homomorphic_filtering(img_fname,
                          filter_shape=(11,11),
                          sigma_x=1.5,
                          sigma_y=1.5):

    assert path.isfile(img_fname)
    filter_shape = tuple(filter_shape)
    assert len(filter_shape) == 2
    assert filter_shape[0] > 0
    assert filter_shape[1] > 0
    assert sigma_x > 0.
    assert sigma_y > 0.

    # -- readin the image into a numpy array
    img = io.imread(img_fname, as_grey=True)

    # -- some parameters
    h0, w0 = img.shape[:2]
    fh, fw = filter_shape

    # -- transforming to floating points and rescaling intensity
    img = skimage.img_as_float(img)
    img = rescale_intensity(img)
    img = np.where(img < EPSILON, EPSILON, img)

    # -- taking the log of the image
    img = np.log(img)

    # -- getting the Laplacian of Gaussian filter
    log_f = generate_log(filter_shape, sigma_x=sigma_x, sigma_y=sigma_y)

    # -- correlation
    nimg = view_as_windows(img, filter_shape)
    nimg = nimg.reshape(-1, fh * fw)
    nimg = np.dot(nimg, log_f.ravel())
    nimg = nimg.reshape((h0 - fh + 1, w0 - fw + 1))

    # -- taking the image exponential
    img = np.exp(nimg)

    # -- finally we perform an histogram equalization
    img = equalize(img)

    return img
