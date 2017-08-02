#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.rgba Contains the RGBAImage class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import imageio

# Import the relevant PTS classes and modules
from ...core.basics.colour import parse_colour

# -----------------------------------------------------------------

class RGBAImage(object):

    """
    This class...
    """

    def __init__(self, data, **kwargs):

        """
        The constructor ...
        :param data:
        :param kwargs:
        """

        # Attributes
        self.data = data
        self.path = kwargs.pop("path", None)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        data = imageio.imread(path)
        return cls(data, path=path)

    # -----------------------------------------------------------------

    @classmethod
    def from_frame(cls, frame, interval="pts", scale="log", alpha=True, peak_alpha=1., colours="red", absolute_alpha=False):

        """
        This function ...
        :param frame:
        :return:
        """

        data = frame_to_rgba(frame, interval=interval, scale=scale, alpha=alpha, peak_alpha=peak_alpha, colours=colours, absolute_alpha=absolute_alpha)
        return cls(data)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return (self.data.shape[0], self.data.shape[1])

    # -----------------------------------------------------------------

    @property
    def red(self):

        """
        This function ...
        :return:
        """

        return self.data[:, :, 0]

    # -----------------------------------------------------------------

    @property
    def green(self):

        """
        This function ...
        :return:
        """

        return self.data[:, :, 1]

    # -----------------------------------------------------------------

    @property
    def blue(self):

        """
        This function ...
        :return:
        """

        return self.data[:, :, 2]

    # -----------------------------------------------------------------

    @property
    def alpha(self):

        """
        This function ...
        :return:
        """

        return self.data[:, :, 3]

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Write
        imageio.imwrite(path, self.data)

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

def frame_to_rgba(frame, interval="pts", scale="log", alpha=True, peak_alpha=1., colours="red", absolute_alpha=False):

    """
    This function ...
    :param frame:
    :param interval:
    :param scale:
    :param alpha:
    :param peak_alpha:
    :param colours:
    :param absolute_alpha:
    :return:
    """

    # Import standard modules
    import numpy as np
    from matplotlib.cm import get_cmap

    # Import astronomical modules
    from astropy.visualization import SqrtStretch, LogStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from astropy.visualization import MinMaxInterval, ZScaleInterval

    # -----------------------------------------------------------------

    # Get data and replace nans and infs
    data = np.copy(frame.data)
    # data[np.isnan(data)] = 0.0
    # data[np.isinf(data)] = 0.0
    data[np.isinf(data)] = float("nan")

    # FLIP UP-DOWN
    data = np.flipud(data)

    # INTERVAL
    if interval == "zscale":
        vmin, vmax = ZScaleInterval().get_limits(data)
    elif interval == "pts":
        # Determine the maximum value in the box and the mimimum value for plotting
        # print("here")
        vmin = max(np.nanmin(data), 0.)
        vmax = 0.5 * (np.nanmax(data) + vmin)
    elif interval == "minmax":
        vmin, vmax = MinMaxInterval().get_limits(data)
    else:
        from ...core.tools import parsing
        try:
            vmin, xmax = parsing.real_tuple(interval)
        except ValueError:
            raise ValueError("Cannot interpret the interval")

    # Normalization
    if scale == "log":
        norm = ImageNormalize(stretch=LogStretch(), vmin=vmin, vmax=vmax)
    elif scale == "sqrt":
        norm = ImageNormalize(stretch=SqrtStretch(), vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Invalid option for 'scale'")

    # Normalize
    normalized = norm(data)

    # Determine transparency
    if absolute_alpha:
        transparency = np.ones_like(normalized)
        transparency[np.isnan(data)] = 0.0
        transparency[data == 0.] = 0.0
    elif alpha:
        transparency = peak_alpha * normalized / np.nanmax(normalized)
    else:
        transparency = np.ones_like(normalized)

    # CREATE THE CHANNEL ARRAYS

    # Red image
    if colours == "red":

        # NxMx4
        red = normalized * 255
        blue = np.zeros_like(red)
        green = np.zeros_like(red)
        alpha = transparency * 255

    # Blue image
    elif colours == "blue":

        red = np.zeros_like(normalized)
        blue = normalized * 255
        green = np.zeros_like(red)
        alpha = transparency * 255

    # Green image
    elif colours == "green":

        red = np.zeros_like(normalized)
        blue = np.zeros_like(normalized)
        green = normalized * 255
        alpha = transparency * 255

    # More intricate colour or colour map
    else:

        # Try to parse the colour
        try:

            colour = parse_colour(colours)
            red = normalized * colour.red
            green = normalized * colour.green
            blue = normalized * colour.blue
            alpha = transparency * 255

        # Assume colour map
        except ValueError:

            # Get the colour map
            cmap = get_cmap(colours)
            rgba = cmap(normalized)
            red = rgba[:, :, 0] * 255
            green = rgba[:, :, 1] * 255
            blue = rgba[:, :, 2] * 255
            alpha = transparency * 255

    # MAKE THE IMAGE ARRAY
    # Stack, create the image array
    arrays = [red, green, blue, alpha]
    image = np.stack(arrays, axis=-1)

    # Return
    return image

# -----------------------------------------------------------------
