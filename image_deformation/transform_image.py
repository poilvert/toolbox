#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
Given a transformation function and an image, this program computes
the transformed image.
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline


def transform_2D_greyscale_image(img, inverse_transformation,
                                 fill_val=0.):
    """
    ``inverse_transformation`` is expected to take two 1D array as input
    and return two 1D array as output
    """

    assert img.ndim == 2
    h, w = img.shape
    img_dtype = img.dtype

    X, Y = np.mgrid[0.:h, 0.:w].astype(np.float32)
    x, y = X.ravel(), Y.ravel()

    # -- first we build an interpolation object from the intensity
    #    function of the input image
    int_fct = RectBivariateSpline(np.arange(h).astype(np.float32),
                                  np.arange(w).astype(np.float32),
                                  img.astype(np.float32))

    # -- then we find the set of "valid" points (by "valid" we mean
    #    the set of points which are images of points in the window
    #    [0, h - 1] * [0, w - 1] by the transformation)
    mask, xa, ya = _get_valid_points(x, y, inverse_transformation, h, w)

    # -- we finally compute the transformed image
    final_img = (fill_val * np.ones((h, w))).astype(np.float32)
    final_img[mask] = int_fct.ev(xa, ya)

    return final_img.astype(img_dtype)


def _get_valid_points(x, y, inverse_transformation, h, w):
    """
    This routine computes all the images of (x, y) by the transformation
    ``inverse_transformation``. Then it selects which are the points for
    which the computed images are within a given window [0, h - 1] * [0, w - 1]

    ``inverse_transformation`` is expected to take two 1D array as input
    and return two 1D array as output
    """

    out_x, out_y = inverse_transformation(x, y)

    x_min, x_max = np.float32(0.), np.float32(h - 1)
    y_min, y_max = np.float32(0.), np.float32(w - 1)

    x_mask = (out_x >= x_min) * (out_x <= x_max)
    y_mask = (out_y >= y_min) * (out_y <= y_max)

    mask_1d = x_mask * y_mask
    mask = mask_1d.reshape(int(h), int(w))

    return mask, out_x[mask_1d], out_y[mask_1d]
