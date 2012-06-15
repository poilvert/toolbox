#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
Given a transformation function, this program computes
the transformed images.
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline


def transform_images(imgs, transformation):

    assert imgs.ndim == 4
    ni, h, w, d = imgs.shape

    x = np.arange(h).astype('f')
    y = np.arange(w).astype('f')

    # -- first we compute the new coordinates


    if d == 1:

        for i in xrange(ni):

            # -- get image
            img = imgs[i, :, :, 0]

            # -- create interpolated function
            fct = RectBivariateSpline(x, y, img)

            # --
