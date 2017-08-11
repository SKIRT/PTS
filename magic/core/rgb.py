#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.rgb Contains the RGBImage class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import imageio
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.vector import PixelShape

# -----------------------------------------------------------------

class RGBImage(object):

    """
    This class...
    """

    def __init__(self, red, green, blue, **kwargs):

        """
        The constructor ...
        :param red:
        :param green:
        :param blue:
        :param kwargs:
        """

        # Set the colour channels
        self.red = red
        self.green = green
        self.blue = blue

        # Other attributes
        self.path = kwargs.pop("path", None)

    # -----------------------------------------------------------------

    @classmethod
    def from_array(cls, array):

        """
        This function ...
        :param array:
        :return:
        """

        # Check
        if array.shape[-1] != 3: raise ValueError("Not an RGB image")

        # Get channels
        red = array[:, :, 0]
        green = array[:, :, 1]
        blue = array[:, :, 2]

        # Create RGB image
        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Read image
        data = imageio.imread(path)

        # Create from 3D array
        image = cls.from_array(data)

        # Set path
        image.path = path

        # Return the image
        return image

    # -----------------------------------------------------------------

    @classmethod
    def from_frame(cls, frame, interval="pts", scale="log", peak_alpha=1., colours="red"):

        """
        This function ...
        :param frame:
        :param interval:
        :param scale:
        :param peak_alpha:
        :param colours:
        :return:
        """

        from .rgba import frame_to_components
        red, green, blue, alpha = frame_to_components(frame, interval=interval, scale=scale, alpha=False, peak_alpha=peak_alpha, colours=colours, absolute_alpha=False)
        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return PixelShape(self.red.shape[0], self.red.shape[1])

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.shape.x

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.shape.y

    # -----------------------------------------------------------------

    def invert_colors(self):

        """
        This function ...
        :return:
        """

        # RED
        self.invert_red()

        # GREEN
        self.invert_green()

        # BLUE
        self.invert_blue()

    # -----------------------------------------------------------------

    def invert_red(self):

        """
        This function ...
        :return:
        """

        self.red = 255 - self.red

    # -----------------------------------------------------------------

    def invert_green(self):

        """
        This function ...
        :return:
        """

        self.green = 255 - self.green

    # -----------------------------------------------------------------

    def invert_blue(self):

        """
        This function ...
        :return:
        """

        self.blue = 255 - self.blue

    # -----------------------------------------------------------------

    def asarray(self):

        """
        Thisf unction ...
        :return:
        """

        return components_to_rgb(self.red, self.green, self.blue)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        import matplotlib.pyplot as plt
        plt.imshow(self.asarray())
        plt.show()

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Write
        imageio.imwrite(path, self.asarray())

        # Update path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        self.saveto(self.path)

# -----------------------------------------------------------------

def components_to_rgb(red, green, blue):

    """
    This function ...
    :param red:
    :param green:
    :param blue:
    :param alpha:
    :return:
    """

    # MAKE THE IMAGE ARRAY
    # Stack, create the image array
    arrays = [red, green, blue]
    image = np.stack(arrays, axis=-1)

    # Return
    return image

# -----------------------------------------------------------------
