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
from ...core.basics.colour import parse_colour

# -----------------------------------------------------------------

class DamagedImageFileError(Exception):

    """
    This class ...
    """

    def __init__(self, message, path=None):

        """
        Thisf unction ...
        :param message:
        :param path:
        """

        # Call the base class constructor with the parameters it needs
        super(DamagedImageFileError, self).__init__(message)

        # The image file path
        self.path = path

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
        try: data = imageio.imread(path)
        except ValueError as e:
            message = str(e)
            if "image file is truncated" in message:
                raise DamagedImageFileError("Image file is damaged", path=path)
            else: raise e

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
        red, green, blue, alpha = frame_to_components(frame, interval=interval, scale=scale, alpha=None, peak_alpha=peak_alpha, colours=colours)
        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @classmethod
    def from_mask(cls, mask, colour="black", background_color="white"):

        """
        This function ...
        :param mask:
        :param colour:
        :param background_color:
        :return:
        """

        red, green, blue = mask_to_components(mask, colour=colour, background_color=background_color)
        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @classmethod
    def from_alpha_mask(cls, mask, colour="black", background_color="white"):

        """
        This function ...
        :param mask:
        :param colour:
        :param background_color:
        :return:
        """

        red, green, blue = alpha_mask_to_components(mask, colour=colour, background_color=background_color)
        return cls(red, green, blue)

    # -----------------------------------------------------------------

    @property
    def shape(self):
        return PixelShape(self.red.shape[0], self.red.shape[1])

    # -----------------------------------------------------------------

    @property
    def xsize(self):
        return self.shape.x

    # -----------------------------------------------------------------

    @property
    def ysize(self):
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

def make_components(shape, colour):

    """
    This function ...
    :param shape:
    :param colour:
    :return:
    """

    # Create all ones
    ones = np.ones(shape, dtype=np.uint8)

    colour = parse_colour(colour)
    red = ones * colour.red
    green = ones * colour.green
    blue = ones * colour.blue

    # Return the components
    return red, green, blue

# -----------------------------------------------------------------

def mask_to_components(mask, colour="black", background_color="white"):

    """
    This function ...
    :param mask:
    :param colour:
    :param background_color:
    :return:
    """

    # FLIP UP-DOWN
    data = np.flipud(mask.data)

    if colour == "red":

        red = data.astype(np.uint8) * 255
        green = np.zeros_like(data, dtype=np.uint8)
        blue = np.zeros_like(data, dtype=np.uint8)

    elif colour == "green":

        red = np.zeros_like(data, dtype=np.uint8)
        green = data.astype(np.uint8) * 255
        blue = np.zeros_like(data, dtype=np.uint8)

    elif colour == "blue":

        red = np.zeros_like(data, dtype=np.uint8)
        green = np.zeros_like(data, dtype=np.uint8)
        blue = data.astype(np.uint8) * 255

    else:

        ones = np.ones_like(data, dtype=np.uint8)

        colour = parse_colour(colour)
        red = ones * colour.red
        green = ones * colour.green
        blue = ones * colour.blue

    #inverted = mask.inverse()
    inverted = np.logical_not(data)
    background_color = parse_colour(background_color)
    red[inverted] = background_color.red
    green[inverted] = background_color.green
    blue[inverted] = background_color.blue

    # Return the components
    return red, green, blue

# -----------------------------------------------------------------

def alpha_mask_to_components(mask, colour="black", background_color="white"):

    """
    This function ...
    :param mask:
    :param colour:
    :param background_color:
    :return:
    """

    # FLIP UP-DOWN
    data = np.flipud(mask.data)
    inverted = 255 - data

    background_color = parse_colour(background_color)

    if colour == "red":

        red = background_color.red * inverted + data
        green = background_color.green * inverted
        blue = background_color.blue * inverted

    elif colour == "green":

        red = background_color.red * inverted
        green = background_color.green * inverted + data
        blue = background_color.blue * inverted

    elif colour == "blue":

        red = background_color.red * inverted
        green = background_color.green * inverted
        blue = background_color.blue * inverted + data

    else:

        ones = np.ones_like(data, dtype=np.uint8)
        colour = parse_colour(colour)

        red = background_color.red * inverted + colour.red * data
        green = background_color.green * inverted + colour.green * data
        blue = background_color.blue * inverted + colour.blue * data

    # Return the components
    return red, green, blue

# -----------------------------------------------------------------
