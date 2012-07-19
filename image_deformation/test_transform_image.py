#!/usr/bin/env python

# Authors: Nicolas Pinto <nicolas.pinto@gmail.com>
#          Nicolas Poilvert <nicolas.poilvert@gmail.com>

# Licence: BSD

"""
Test suite for ``transform_image``
"""

import skimage
import skimage.io as io

from transform_library import inverse_translate
from transform_library import inverse_inplane_rotate
from transform_library import inverse_zoom
from transform_library import inverse_fovea_like
from transform_image import transform_2D_greyscale_image

DEFAULT_OBAMA_PICTURE = 'Obama.png'


def test_translate_Obama():

    Obama_greyscale = skimage.img_as_float(
                          io.imread(DEFAULT_OBAMA_PICTURE, as_grey=True)
                          )
    io.imshow(Obama_greyscale)
    io.show()

    translated_Obama = transform_2D_greyscale_image(Obama_greyscale,
                                                    inverse_translate)
    io.imshow(translated_Obama)
    io.show()

    return


def test_inplane_rotate_Obama():

    Obama_greyscale = skimage.img_as_float(
                          io.imread(DEFAULT_OBAMA_PICTURE, as_grey=True)
                          )
    io.imshow(Obama_greyscale)
    io.show()

    rotated_Obama = transform_2D_greyscale_image(Obama_greyscale,
                                                    inverse_inplane_rotate)
    io.imshow(rotated_Obama)
    io.show()

    return


def test_zoom_Obama():

    Obama_greyscale = skimage.img_as_float(
                          io.imread(DEFAULT_OBAMA_PICTURE, as_grey=True)
                          )
    io.imshow(Obama_greyscale)
    io.show()

    zoomed_Obama = transform_2D_greyscale_image(Obama_greyscale,
                                                    inverse_zoom)
    io.imshow(zoomed_Obama)
    io.show()

    return


def test_fovea_Obama():

    Obama_greyscale = skimage.img_as_float(
                          io.imread(DEFAULT_OBAMA_PICTURE, as_grey=True)
                          )
    io.imshow(Obama_greyscale)
    io.show()

    foveated_Obama = transform_2D_greyscale_image(Obama_greyscale,
                                                    inverse_fovea_like)
    io.imshow(foveated_Obama)
    io.show()

    return


if __name__ == "__main__":
    #test_translate_Obama()
    #test_inplane_rotate_Obama()
    #test_zoom_Obama()
    test_fovea_Obama()
