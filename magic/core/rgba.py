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
import numpy as np

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.basics.colour import parse_colour
from ..dist_ellipse import distance_ellipse
from ..basics.vector import PixelShape
from ..tools import plotting

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

        # Read image
        data = imageio.imread(path)

        # CHECK WHETHER RGBA
        #print(data)
        #print(data.shape)
        if data.shape[-1] != 4: data = rgb_to_rgba(data)

        #print(data[:,:,3])

        # Create RGBA image
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

        #return (self.data.shape[0], self.data.shape[1])
        return PixelShape(self.data.shape[0], self.data.shape[1])

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

    def invert_colors(self):

        """
        This function ...
        :return:
        """

        self.data[:, :, 0:3] = 255 - self.data[:, :, 0:3]

    # -----------------------------------------------------------------

    def soften_edges(self, region, factor_range):

        """
        This function ...
        :param region:
        :param range:
        :param factor_range:
        :return:
        """

        #if not

        #astro = data.astronaut()

        #center = Position(self.xsize, self.ysize)
        center = region.center

        angle = - region.angle + Angle(-90., "deg")
        #angle = region.angle

        # Determine the ratio of semimajor and semiminor
        ratio = region.semiminor / region.semimajor
        radius = distance_ellipse(self.shape, center, ratio, angle) / region.semiminor

        #mask = region.to_mask(self.xsize, self.ysize)
        #plotting.plot_box(radius)
        #plotting.plot_mask(radius > 1)
        #plotting.plot_mask(mask)
        #plotting.plot_mask(radius < 1)

        #rows, cols = np.mgrid[:l_row, :l_col]
        #radius = np.sqrt((rows - l_row / 2) ** 2 + (cols - l_col / 2) ** 2)

        #r_min, r_max = factor_range.min * radius.max(), factor_range.max * radius.max()
        #r_min, r_max = factor_range.min, factor_range.max
        #print(r_min, r_max)

        outside_max = radius > factor_range.max
        inside_min = radius < factor_range.min

        #plotting.plot_mask(outside_max, title="outside max")
        #plotting.plot_mask(inside_min, title="inside min")

        # METHOD THAT WAS SUPPOSED TO WORK:
        # l_row, l_col, nb_channel = self.data.shape
        # alpha_channel = np.zeros((l_row, l_col), dtype=np.uint8)
        #
        # alpha_channel[inside_min] = 1
        # alpha_channel[outside_max] = 0
        #
        # #plotting.plot_box(alpha_channel)
        # #plotting.plot_mask(alpha_channel)
        #
        # gradient_zone = np.logical_and(radius >= factor_range.min, radius <= factor_range.max)
        #
        # #plotting.plot_mask(gradient_zone)
        #
        test = (factor_range.max - radius) / factor_range.span
        # #plotting.plot_box(test)
        # plotting.plot_mask(test < 1)
        #
        # #alpha_channel[gradient_zone] = (factor_range.max - radius[gradient_zone]) / factor_range.span
        # alpha_channel[gradient_zone] = (factor_range.max - radius[gradient_zone])
        # #print(np.min(alpha_channel), np.max(alpha_channel))
        #
        # #plotting.plot_box(alpha_channel)
        #
        # alpha_channel *= 255


        #alpha_channel = np.zeros((l_row, l_col), dtype=np.uint8)

        alpha_channel = test
        alpha_channel[inside_min] = 1
        alpha_channel[outside_max] = 0
        #alpha_channel *= 255
        #alpha_channel = alpha_channel.astype(np.uint8)

        #plotting.plot_box(alpha_channel)
        #print(alpha_channel)

        #feathered = np.empty((l_row, l_col, nb_channel + 1), dtype=np.uint8)

        #print(self.alpha)

        #self.show()
        #mask = self.alpha == 255
        #plotting.plot_mask(mask)
        #plotting.plot_mask(self.alpha == 255, title="alpha")

        # MULTIPLY WITH THE ORIGINAL ALPHA CHANNEL (AVOID MAKING ALREADY TRANSPARENT SECTIONS SUDDENLY NOT TRANSPARENT ANYMORE)
        alpha_channel = self.alpha * alpha_channel

        #feathered[..., :3] = astro[:]
        #feathered[..., -1] = alpha_channel[:]
        #self.data[..., -1] = alpha_channel[:]
        self.data[..., -1] = alpha_channel

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        import matplotlib.pyplot as plt

        import matplotlib.image as mpimg
        #img = mpimg.imread('file-name.png')

        plt.imshow(self.data)
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
    if scale == "log": norm = ImageNormalize(stretch=LogStretch(), vmin=vmin, vmax=vmax)
    elif scale == "sqrt": norm = ImageNormalize(stretch=SqrtStretch(), vmin=vmin, vmax=vmax)
    else: raise ValueError("Invalid option for 'scale'")

    # Normalize
    normalized = norm(data)

    # Determine transparency
    if absolute_alpha:
        transparency = np.ones_like(normalized, dtype=np.uint8)
        transparency[np.isnan(data)] = 0
        transparency[data == 0] = 0
    elif alpha: transparency = peak_alpha * normalized / np.nanmax(normalized)
    else: transparency = np.ones_like(normalized, dtype=np.uint8)

    # CREATE THE CHANNEL ARRAYS

    # Red image
    if colours == "red":

        # NxMx4
        red = normalized * 255
        blue = np.zeros_like(red, dtype=np.uint8)
        green = np.zeros_like(red, dtype=np.uint8)
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

    # Make the image
    return components_to_rgba(red, green, blue, alpha)

# -----------------------------------------------------------------

def rgb_to_rgba(rgb):

    """
    This function ...
    :param rgb:
    :return:
    """

    red = rgb[:, :, 0]
    green = rgb[:, :, 1]
    blue = rgb[:, :, 2]

    # Create RGBa
    return components_to_rgba(red, green, blue)

# -----------------------------------------------------------------

def components_to_rgba(red, green, blue, alpha=None):

    """
    This function ...
    :param red:
    :param green:
    :param blue:
    :param alpha:
    :return:
    """

    # Create alpha channel if not defined
    if alpha is None: alpha = 255 * np.ones_like(red, dtype=np.uint8)

    # MAKE THE IMAGE ARRAY
    # Stack, create the image array
    arrays = [red, green, blue, alpha]
    image = np.stack(arrays, axis=-1)

    # Return
    return image

# -----------------------------------------------------------------
