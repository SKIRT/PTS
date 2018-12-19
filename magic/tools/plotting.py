#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.plotting Contains convenient plotting functions.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import warnings
from abc import ABCMeta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.widgets import Slider
from collections import OrderedDict
from scipy.stats.kde import gaussian_kde
from matplotlib import markers
from matplotlib.path import Path

# Import astronomical modules
from astropy.visualization.stretch import SqrtStretch, LogStretch, LinearStretch, HistEqStretch, AsinhStretch
from astropy.visualization.stretch import SinhStretch, PowerStretch, PowerDistStretch, InvertedPowerDistStretch
from astropy.visualization.stretch import SquaredStretch, InvertedLogStretch, InvertedHistEqStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import MinMaxInterval, ZScaleInterval
from photutils import CircularAperture

# Import the relevant PTS classes and modules
from ...core.tools import types, sequences
from ...core.basics.log import log
from ...core.basics import containers
from . import scatter_density # NOQA

# -----------------------------------------------------------------

colours = ['Pink','LightPink','HotPink','DeepPink','PaleVioletRed','MediumVioletRed','Red','LightSalmon','Salmon',
           'DarkSalmon','LightCoral','IndianRed','Crimson','FireBrick','DarkRed','Red','Orange','OrangeRed','Tomato',
           'Coral','DarkOrange','Orange','Gold','Yellow','Yellow','LightYellow','LemonChiffon','LightGoldenrodYellow',
           'PapayaWhip','Moccasin','PeachPuff','PaleGoldenrod','Khaki','DarkKhaki','Brown','Cornsilk','BlanchedAlmond',
           'Bisque','NavajoWhite','Wheat','BurlyWood','Tan','RosyBrown','SandyBrown','Goldenrod','DarkGoldenrod','Peru',
           'Chocolate','SaddleBrown','Sienna','Brown','Maroon','DarkOliveGreen','Olive','OliveDrab','YellowGreen','LimeGreen',
           'Lime','LawnGreen','Chartreuse','GreenYellow','SpringGreen','MediumSpringGreen','LightGreen','PaleGreen','DarkSeaGreen',
           'MediumSeaGreen','SeaGreen','ForestGreen','Green','DarkGreen','Cyan','MediumAquamari','Aqua','Cyan','LightCyan',
           'PaleTurquoise','Aquamarine','Turquoise','MediumTurquois','DarkTurquoise','LightSeaGreen','CadetBlue','DarkCyan',
           'Teal','LightSteelBlue','PowderBlue','LightBlue','SkyBlue','LightSkyBlue','DeepSkyBlue','DodgerBlue',
           'CornflowerBlue','SteelBlue','RoyalBlue','Blue','MediumBlue','DarkBlue','Navy','MidnightBlue','Lavender',
           'Thistle','Plum','Violet','Orchid','Fuchsia','Magenta','MediumOrchid','MediumPurple','BlueViolet','DarkViolet',
           'DarkOrchid','DarkMagenta','Purple','Indigo','DarkSlateBlue','SlateBlue','MediumSlateBlue','White','Snow',
           'Honeydew','MintCream','Azure','AliceBlue','GhostWhite','WhiteSmoke','Seashell','Beige','OldLace','FloralWhite',
           'Ivory','AntiqueWhite','Linen','LavenderBlush','MistyRose','Gainsboro','LightGray','Silver','DarkGray','Gray',
           'DimGray','LightSlateGray','SlateGray','DarkSlateGray','Black']

pretty_colours = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']

line_styles = ['-', '--', '-.', ':']

# -----------------------------------------------------------------

# DEFINE COLOR MAPS FROM DS9

from matplotlib.cm import register_cmap, cmap_d

ds9a = {'red': lambda v : np.interp(v, [0, 0.25, 0.5, 1],
                                        [0, 0, 1, 1]),
         'green': lambda v : np.interp(v, [0, 0.25, 0.5, 0.77, 1],
                                          [0, 1, 0, 0, 1]),
         'blue': lambda v : np.interp(v, [0, 0.125, 0.5, 0.64, 0.77, 1],
                                         [0, 0, 1, 0.5, 0, 0])}

ds9b = {'red': lambda v : 4 * v - 1,
        'green': lambda v : 4 * v - 2,
        'blue': lambda v : np.select([v < 0.25, v < 0.5, v < 0.75, v <= 1],
                                      [4 * v, -4 * v + 2, 0, 4 * v - 3])}

# Note that this definition slightly differs from ds9cool, but make more sense to me...
ds9cool = {'red': lambda v : 2 * v - 1,
           'green': lambda v : 2 * v - 0.5,
           'blue': lambda v : 2 * v}

ds9i8 = {'red': lambda v : np.where(v < 0.5, 0, 1),
        'green': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                       v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 1, 0, 1, 0, 1, 0, 1]),
        'blue': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                      v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 0, 1, 1, 0, 0, 1, 1])}

ds9aips0 = {'red': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.475, 0, 0.373, 0, 0, 1, 1, 1]),
            'green': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0, 0, 0.655, 0.596, 0.965, 1, 0.694, 0]),
            'blue': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.608, 0.785, 0.925, 0, 0, 0, 0, 0])}

ds9rainbow = {'red': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [1, 0, 0, 1, 1]),
              'green': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [0, 0, 1, 1, 0]),
              'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [1, 1, 0, 0])}

# This definition seems a bit strange...
ds9he = {'red': lambda v : np.interp(v, [0, 0.015, 0.25, 0.5, 1],
                                        [0, 0.5, 0.5, 0.75, 1]),
         'green': lambda v : np.interp(v, [0, 0.065, 0.125, 0.25, 0.5, 1],
                                          [0, 0, 0.5, 0.75, 0.81, 1]),
         'blue': lambda v : np.interp(v, [0, 0.015, 0.03, 0.065, 0.25, 1],
                                         [0, 0.125, 0.375, 0.625, 0.25, 1])}

ds9heat = {'red': lambda v : np.interp(v, [0, 0.34, 1], [0, 1, 1]),
           'green': lambda v : np.interp(v, [0, 1], [0, 1]),
           'blue': lambda v : np.interp(v, [0, 0.65, 0.98, 1], [0, 0, 1, 1])}

# Set aliases, where colormap exists in matplotlib
cmap_d['ds9bb'] = cmap_d['afmhot']
cmap_d['ds9grey'] = cmap_d['gray']

# Register all other colormaps
register_cmap('ds9b', data=ds9b)
register_cmap('ds9cool', data=ds9cool)
register_cmap('ds9a', data=ds9a)
register_cmap('ds9i8', data=ds9i8)
register_cmap('ds9aips0', data=ds9aips0)
register_cmap('ds9rainbow', data=ds9rainbow)
register_cmap('ds9he', data=ds9he)
register_cmap('ds9heat', data=ds9heat)

# -----------------------------------------------------------------

def get_stretch(name, data=None, parameter=None):

    """
    This function ...
    :param name:
    :param data:
    :param parameter:
    :return:
    """

    if name == "auto":

        if data is None: raise ValueError("Data must be passed when scale = 'auto'")

        #from ...core.basics.distribution import Distribution
        #from ...core.plot.distribution import plot_distribution

        #distr = Distribution.from_data("linear", data)
        #plot_distribution(distr, title="linear")

        mean = np.nanmean(data)
        median = np.nanmedian(data)
        stddev = np.nanstd(data)
        if mean == 0 or median == 0: reldiff = 1
        else: reldiff = abs(mean - median) / abs(mean)

        logdata = np.log(data)
        #logdistr = Distribution.from_data("log", logdata)
        #plot_distribution(logdistr, title="log")

        isvalid = np.isfinite(logdata)
        logmean = np.nanmean(logdata[isvalid])
        logmedian = np.nanmedian(logdata[isvalid])
        logstddev = np.nanstd(logdata[isvalid])
        if logmean == 0 or logmedian == 0: logreldiff = 1
        else: logreldiff = abs(logmean - logmedian) / abs(logmean)

        #print(mean, median, stddev, reldiff)
        #print(logmean, logmedian, logstddev, logreldiff)

        # Log or linear?
        if logreldiff < reldiff: return LogStretch()
        else: return LinearStretch()

    elif name == "log":
        if parameter is not None: stretch = LogStretch(a=parameter)
        else: stretch = LogStretch() # default
    elif name == "sqrt": stretch = SqrtStretch()
    elif name == "linear": stretch = LinearStretch()
    elif name == "histeq":
        if data is None: raise ValueError("Data has to be passed for histogram equalization scale")
        stretch = HistEqStretch(data)
    elif name == "asinh":
        if parameter is not None: stretch = AsinhStretch(parameter)
        else: stretch = AsinhStretch() # default
    elif name == "sinh":
        if parameter is not None: stretch = SinhStretch(parameter)
        else: stretch = SinhStretch() # default
    elif name == "power":
        if parameter is None: raise ValueError("Parameter has to be passed for power stretch")
        stretch = PowerStretch(a=parameter)
    elif name == "powerdist":
        if parameter is not None: stretch = PowerDistStretch(a=parameter)
        else: stretch = PowerDistStretch() # default
    elif name == "invpowerdist":
        if parameter is not None: stretch = InvertedPowerDistStretch(a=parameter)
        else: stretch = InvertedPowerDistStretch() # default
    elif name == "squared": stretch = SquaredStretch()
    elif name == "invlog":
        if parameter is None: raise ValueError("Parameter has to be passed for inverted log stretch")
        stretch = InvertedLogStretch(parameter)
    elif name == "invhisteq":
        if data is None: raise ValueError("Data has to be passed for inverted hist equalization scale")
        stretch = InvertedHistEqStretch(data)
    else: raise ValueError("Invalid scale '" + name + "'")

    # Return
    return stretch

# -----------------------------------------------------------------

def get_normalization(scale_name, vmin, vmax, data=None, scale_parameter=None):

    """
    Thins function ...
    :param scale_name:
    :param vmin:
    :param vmax:
    :param data:
    :param scale_parameter:
    :return:
    """

    # Get stretch
    stretch = get_stretch(scale_name, data=data, parameter=scale_parameter)

    # Create normalization
    norm = ImageNormalize(stretch=stretch, vmin=vmin, vmax=vmax)

    # Return the normalization
    return norm

# -----------------------------------------------------------------

def color_dict(gradient):

    """
    Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on
    """

    from ...core.basics.colour import rgb_to_hex

    return {"hex":[rgb_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}

# -----------------------------------------------------------------

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=307):

    """
    Returns a gradient list of (n) colors between
        two hex colors. start_hex and finish_hex
        should be the full six-digit color string,
        inlcuding the number sign ("#FFFFFF")
    """

    from ...core.basics.colour import hex_to_rgb

    # Starting and ending colors in RGB form
    s = hex_to_rgb(start_hex)
    f = hex_to_rgb(finish_hex)

    # Initilize a list of the output colors with the starting color
    RGB_list = [s]

    # Calculate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):

        # Interpolate RGB vector for color at the current value of t
        curr_vector = [int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) for j in range(3)]

        # Add it to our list of output colors
        RGB_list.append(curr_vector)

    # Return
    return RGB_list

# -----------------------------------------------------------------

# Have colormaps separated into categories:
# http://matplotlib.org/examples/color/colormaps_reference.html
cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
         ('Sequential (2)', [
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper']),
         ('Diverging', [
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
         ('Qualitative', [
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3']),
         ('Miscellaneous', [
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]

# -----------------------------------------------------------------

def plot_color_gradients(cmap_category, cmap_list, nrows):

    """
    This function ...
    :param cmap_category:
    :param cmap_list:
    :param nrows:
    :return:
    """

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    fig, axes = plt.subplots(nrows=nrows)
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    axes[0].set_title(cmap_category + ' colormaps', fontsize=14)

    # Loop over
    for ax, name in zip(axes, cmap_list):

        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes: ax.set_axis_off()

# -----------------------------------------------------------------

def show_gradients():

    """
    This function ...
    :return
    """

    nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps)

    for cmap_category, cmap_list in cmaps:
        plot_color_gradients(cmap_category, cmap_list, nrows)
    plt.show()

# -----------------------------------------------------------------

def plot_coordinates_on_image(image, x_coordinates, y_coordinates, path=None, format=None):

    """
    This function ...
    :param dimensions:
    :param x_coordinates:
    :param y_coordinates:
    :return:
    """

    plt.figure(figsize=(15, 15))

    plt.imshow(image, origin="lower", interpolation="nearest")
    plt.plot(x_coordinates, y_coordinates, 'ro')

    if path is None: plt.show()
    else: plt.savefig(path, format=format)

    plt.close()

# -----------------------------------------------------------------

def plot_mask_on_image(image, mask, title=None, path=None, format=None):

    """
    This function ...
    :param image:
    :param mask:
    :param title:
    :param path:
    :param format:
    :return:
    """

    pass

# -----------------------------------------------------------------

def plot_mask(mask, **kwargs):

    """
    This function ...
    :param mask:
    :param kwargs:
    :return:
    """

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Create prepared mask
    mask, kwargs = create_prepared_mask(mask, **kwargs)

    # Get properties
    title = kwargs.pop("title", None)
    path = kwargs.pop("path", None)
    format = kwargs.pop("format", None)
    show_axes = kwargs.pop("show_axes", True)
    transparent = kwargs.pop("transparent", False)
    axes = kwargs.pop("axes", None)
    cmap = kwargs.pop("colormap", "Greys")
    xsize = kwargs.pop("xsize", 7)
    ysize = kwargs.pop("ysize", 7)

    # Get the data
    if isinstance(mask, np.ndarray): maskdata = mask
    elif hasattr(mask, "data"): maskdata = mask.data
    else: maskdata = np.asarray(mask)

    # Check if not completely masked or completely unmasked
    all_masked = np.all(maskdata)
    all_unmasked = np.all(np.logical_not(maskdata))
    if all_masked:
        warnings.warn("The mask is completely filled")
        cmap = cmap + "_r"
    if all_unmasked: warnings.warn("The mask is completely empty")

    nxpixels = maskdata.shape[1]
    nypixels = maskdata.shape[0]

    # Create figure if necessary, get the axes
    only_axes = False
    if axes is None:
        plt.figure(figsize=(xsize, ysize))
        plt.xlim(0, nxpixels - 1)
        plt.ylim(0, nypixels - 1)
        axes = plt.gca()
    else: only_axes = True

    # Plot the mask
    aspect = "equal"
    image = axes.imshow(maskdata, origin="lower", interpolation="nearest", cmap=cmap, aspect=aspect)
    output.image = image

    # Add region
    region = kwargs.pop("region", None)
    regions = kwargs.pop("regions", None)
    if region is not None:

        from ..region.composite import CompositeRegion
        if isinstance(region, CompositeRegion):
            for patch in region.to_mpl_patches(): axes.add_patch(patch)
        else:
            # Add patch
            patch = region.to_mpl_patch()
            axes.add_patch(patch)

    # Add multiple regions
    if regions is not None:
        for patch in regions.to_mpl_patches(): axes.add_patch(patch)

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Show axis?
        if not show_axes: plt.axis('off')

        # Add title
        if title is not None: plt.title(title)
        else: plt.title("Black means True")

        # Show or save as file
        if path is None: plt.show()
        else: plt.savefig(path, format=format, transparent=transparent)

        # Close the figure
        plt.close()

    # Return
    return output

# -----------------------------------------------------------------

def create_prepared_frame(frame, **kwargs):

    """
    This function ...
    :param frame:
    :param kwargs:
    :return:
    """

    # Crop?
    crop_to = kwargs.pop("crop_to", None)
    cropping_factor = kwargs.pop("cropping_factor", 1.)

    # Crop the frame
    if crop_to is not None: frame = frame.cropped_to(crop_to, factor=cropping_factor, out_of_bounds="expand")
    else: frame = frame.copy()

    # Truncate?
    truncate_outside = kwargs.pop("truncate_outside", None)
    if truncate_outside is not None:
        from ..region.region import SkyRegion, PixelRegion
        if isinstance(truncate_outside, PixelRegion): truncate_outside = truncate_outside.to_mask(frame.xsize, frame.ysize)
        elif isinstance(truncate_outside, SkyRegion): truncate_outside = truncate_outside.to_pixel(frame.wcs).to_mask(frame.xsize, frame.ysize)
        mask = truncate_outside.inverse()
        # Mask
        frame[mask] = float("nan")

    # Mask negatives
    if kwargs.pop("mask_negatives", False):
        frame[frame < 0] = float("nan")

    # Return the new frame
    return frame, kwargs

# -----------------------------------------------------------------

def create_prepared_mask(mask, **kwargs):

    """
    This function ...
    :param mask:
    :param kwargs:
    :return:
    """

    from ..core.mask import Mask
    from ..basics.mask import Mask as oldMask

    # Crop?
    crop_to = kwargs.pop("crop_to", None)
    cropping_factor = kwargs.pop("cropping_factor", 1.)

    # Crop the frame
    if crop_to is not None:
        if isinstance(mask, oldMask):
            mask = Mask(mask)
            mask.crop_to(crop_to, factor=cropping_factor, out_of_bounds="expand") # don't make copy again
        else: mask = mask.cropped_to(crop_to, factor=cropping_factor, out_of_bounds="expand")
    else:
        if isinstance(mask, oldMask): mask = Mask(mask)
        else: mask = mask.copy()

    # Truncate?
    truncate_outside = kwargs.pop("truncate_outside", None)
    if truncate_outside is not None:
        from ..region.region import SkyRegion, PixelRegion
        if isinstance(truncate_outside, PixelRegion): truncate_outside = truncate_outside.to_mask(mask.xsize, mask.ysize)
        elif isinstance(truncate_outside, SkyRegion): truncate_outside = truncate_outside.to_pixel(mask.wcs).to_mask(mask.xsize, mask.ysize)
        truncate_mask = truncate_outside.inverse()
        # Mask the mask
        mask[truncate_mask] = False

    # Return the new mask
    return mask, kwargs

# -----------------------------------------------------------------

def plot_datacube(datacube, **kwargs):

    """
    This function ...
    :param datacube:
    :param kwargs:
    :return:
    """

    # Get the frames
    cube = datacube.asarray(axis=0)

    # Set options
    wavelengths = datacube.wavelengths(unit="micron", add_unit=True)
    kwargs["labels"] = [str(wavelength) for wavelength in wavelengths]
    kwargs["slider_name"] = "Wavelength index"

    # Plot
    plot_cube(cube, **kwargs)

# -----------------------------------------------------------------

def plot_image_masks(image, **kwargs):

    """
    This function ...
    :param image:
    :param kwargs:
    :return:
    """

    kwargs["labels"] = image.mask_names
    plot_masks(image.masks, **kwargs)

# -----------------------------------------------------------------

def plot_masks(masks, **kwargs):

    """
    This function ...
    :param masks:
    :param kwargs:
    :return:
    """

    # Get the number of planes
    nmasks = len(masks)
    xsize = masks[0].shape[1]
    ysize = masks[0].shape[0]

    # Get settings
    labels = kwargs.pop("labels", None)
    slider_name = kwargs.pop("slider_name", "Index")
    if labels is not None and len(labels) != nmasks: raise ValueError("Number of labels must be equal to number of masks")
    title = kwargs.pop("title", None)
    show_axes = kwargs.pop("show_axes", True)

    # Create figure
    figsize = (5, 8)
    figure = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1])
    axes = plt.subplot(gs[0])
    sliderax = plt.subplot(gs[1])

    # Plot first mask
    output = plot_mask(masks[0], axes=axes)
    image = output.image

    # Show axis?
    if not show_axes: axes.axis('off')

    # Add title
    if title is not None: axes.set_title(title)

    # Add text
    if labels is not None: t = axes.text(0.5 * xsize, 0.01 * ysize, labels[0], verticalalignment='bottom', horizontalalignment='center', color="white")
    else: t = None

    # Create slider
    slider = DiscreteSlider(sliderax, slider_name, 0, nmasks - 1, increment=1, valinit=0)

    # Create update function
    def update(index):
        image.set_data(masks[index])
        if labels is not None: t.set_text(labels[index])
        figure.canvas.draw_idle()

    # Set update function
    slider.on_changed(update)

    # Tight layout
    plt.tight_layout()

    # Show
    plt.show()

    # Close the figure
    plt.close()

# -----------------------------------------------------------------

class DiscreteSlider(Slider):

    """A matplotlib slider widget with discrete steps."""

    def __init__(self, *args, **kwargs):

        """Identical to Slider.__init__, except for the "increment" kwarg.
        "increment" specifies the step size that the slider will be discritized
        to.
        """

        kwargs["valfmt"] = "%d"

        self.inc = kwargs.pop('increment', 0.5)
        Slider.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    def set_val(self, val):

        """
        This function ...
        :param val:
        :return:
        """

        discrete_val = int(val / self.inc) * self.inc
        # We can't just call Slider.set_val(self, discrete_val), because this
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: return
        for cid, func in self.observers.iteritems(): func(discrete_val)

# -----------------------------------------------------------------

def plot_cube(cube, **kwargs):

    """
    This function ...
    :param kwargs:
    :return:
    """

    # Get the number of planes
    nplanes = cube.shape[0]
    xsize = cube[0].shape[1]
    ysize = cube[1].shape[0]

    # Get settings
    labels = kwargs.pop("labels", None)
    slider_name = kwargs.pop("slider_name", "Index")
    if labels is not None and len(labels) != nplanes: raise ValueError("Number of labels must be equal to number of planes in the cube")
    share_normalization = kwargs.pop("share_normalization", True)
    title = kwargs.pop("title", None)
    colorbar = kwargs.pop("colorbar", True)
    show_axes = kwargs.pop("show_axes", True)
    path = kwargs.pop("path", None)
    transparent = kwargs.pop("transparent", False)

    # Create figure
    figsize = (5, 8)
    figure = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1])
    axes = plt.subplot(gs[0])
    sliderax = plt.subplot(gs[1])

    # Determine intervals
    if not share_normalization:
        intervals = []
        for index in range(nplanes):
            data = cube[index]
            vmin, vmax = get_vmin_vmax(data, **kwargs)
            intervals.append((vmin, vmax))
    else: intervals = None

    # Plot first plane
    output = plot_box(cube[0], axes=axes, **kwargs)
    if share_normalization:
        norm = output.norm
        vmin = output.vmin
        vmax = output.vmax
    else: vmin = vmax = None

    # Show axis?
    if not show_axes: axes.axis('off')

    # Add title
    if title is not None: axes.set_title(title)

    # Add colorbar?
    if colorbar: cb = plt.colorbar(output.image, ax=axes)
    else: cb = None

    # Add text
    if labels is not None: t = axes.text(0.5 * xsize, 0.01 * ysize, labels[0], verticalalignment='bottom', horizontalalignment='center', color="white")
    else: t = None

    # Create slider
    slider = DiscreteSlider(sliderax, slider_name, 0, nplanes-1, increment=1, valinit=0)

    # Create update function
    def update(index):

        if share_normalization: output.image.set_data(cube[index])
        else:
            output.image.set_data(cube[index])
            output.image.set_clim(vmin=intervals[index][0], vmax=intervals[index][1])
            if colorbar:
                cb.set_clim(vmin=intervals[index][0], vmax=intervals[index][1])
                cb.draw_all()

        if labels is not None: t.set_text(labels[index])
        figure.canvas.draw_idle()

    # Set update function
    slider.on_changed(update)

    # Tight layout
    plt.tight_layout()

    # Show or save as file
    if path is None: plt.show()
    else: plt.savefig(path, format=format, transparent=transparent)

    # Close the figure
    plt.close()

# -----------------------------------------------------------------

def plot_frame(frame, **kwargs):

    """
    This function ...
    :param frame:
    :param kwargs:
    :return:
    """

    # Create prepared frame
    frame, kwargs = create_prepared_frame(frame, **kwargs)

    # Plot
    return plot_box(frame, **kwargs)

# -----------------------------------------------------------------

def get_vmin_vmax(data, interval="pts", around_zero=False, symmetric=False, normalize_in=None, soft_min=False,
                  soft_max=False, soft_min_scaling=1., soft_max_scaling=1., symmetric_method="mean",
                  check_around_zero=True, wcs=None, zmin=None, zmax=None, soft_zmin=False, soft_zmax=False,
                  percentile_low=0.5, percentile_high=99.9, logscale=False, minmax_scaling=0.5, min_relative_max_log=1e-4):

    """
    This function ...
    :param data:
    :param interval:
    :param around_zero:
    :param symmetric:
    :param normalize_in:
    :param soft_min:
    :param soft_max:
    :param soft_min_scaling:
    :param soft_max_scaling:
    :param symmetric_method:
    :param check_around_zero:
    :param wcs:
    :param zmin:
    :param zmax:
    :param soft_zmin:
    :param soft_zmax:
    :param percentile_low:
    :param percentile_high:
    :param logscale:
    :param minmax_scaling:
    :param min_relative_max_log:
    :return:
    """

    # Check parameters
    if symmetric and not around_zero: raise ValueError("Cannot enable 'symmetric' but not 'around_zero'")

    # DETERMINE NORMALIZE MIN AND MAX: ONLY FOR PTS INTERVAL METHOD FOR NOW
    from ..region.region import SkyRegion, PixelRegion

    # Normalize_in is passed
    if normalize_in is not None:
        if isinstance(normalize_in, SkyRegion):
            if wcs is None: raise ValueError("Cannot give sky region when the passed data doesn't have a coordinate system")
            normalize_in = normalize_in.to_pixel(wcs)
        if isinstance(normalize_in, PixelRegion): normalize_in = normalize_in.to_mask(data.shape[1], data.shape[0])

        # Get the mask data
        if isinstance(normalize_in, np.ndarray): mask_data = normalize_in
        else: mask_data = normalize_in.data

        pixels = data[mask_data]
        normalize_min = np.nanmin(pixels)
        normalize_max = np.nanmax(pixels)

    # Normalize over the entire image
    else:

        pixels = data.flatten()
        normalize_min = np.nanmin(data)
        normalize_max = np.nanmax(data)

    nnegatives = np.sum(pixels < 0)
    npositives = np.sum(pixels > 0)
    has_negatives = nnegatives > 0
    has_positives = npositives > 0

    # Check around_zero flag
    if check_around_zero and around_zero and not (has_negatives and has_positives): raise ValueError("Not around zero")

    # ZSCALE
    if interval == "zscale":

        vmin, vmax = ZScaleInterval().get_limits(pixels)
        below_zero = vmin < 0
        above_zero = vmax > 0

        # Around zero
        if around_zero and symmetric:
            if check_around_zero and not (below_zero and above_zero): raise ValueError("Not around zero")

            # Set max
            if symmetric_method == "mean": vmax = 0.5 * sum([abs(vmin), abs(vmax)])
            elif symmetric_method == "max": vmax = max([abs(vmin), abs(vmax)])
            elif symmetric_method == "min": vmax = min([abs(vmin), abs(vmax)])
            else: raise ValueError("Invalid symmetric method")

            # Set min
            vmin = - vmax

    # PERCENTILE
    elif interval == "percentile":

        vmin = np.nanpercentile(pixels, percentile_low)
        vmax = np.nanpercentile(pixels, percentile_high)

    # PTS INTERVAL
    elif interval == "pts":

        # Around zero
        if around_zero:

            vmin = minmax_scaling * normalize_min
            vmax = minmax_scaling * normalize_max

            # Symmetric?
            if symmetric:

                # Set max
                if symmetric_method == "mean": vmax = 0.5 * sum([abs(vmin), abs(vmax)]) # mean
                elif symmetric_method == "max": vmax = max([abs(vmin), abs(vmax)])
                elif symmetric_method == "min": vmax = min([abs(vmin), abs(vmax)])
                else: raise ValueError("Invalid symmetric method")

                # Set min
                vmin = - vmax

        # More positive values than negatives
        elif npositives > nnegatives:

            # Determine the maximum value in the box and the mimimum value for plotting
            if logscale:

                vmin = max(normalize_min, min_relative_max_log * normalize_max)
                # DO ALLOW ZERO FOR LOG? -> matplotlib seems to be able to handle it -> keep the same for now
                #vmin = 0.
                vmax = minmax_scaling * (normalize_max + vmin)

            else:

                vmin = max(normalize_min, 0.)
                vmax = minmax_scaling * (normalize_max + vmin)

        # More negative values than positives
        else:

            if logscale: raise ValueError("More negative values than positive values: are you sure you want to be plotting this data on a logarithmic scale?")

            vmax = min(normalize_max, 0.)
            vmin = minmax_scaling * (normalize_min + vmax)

    # MINIMAX
    elif interval == "minmax":

        vmin, vmax = MinMaxInterval().get_limits(pixels)
        below_zero = vmin < 0
        above_zero = vmax > 0

        # Around zero?
        if around_zero and symmetric:
            if check_around_zero and not (below_zero and above_zero): raise ValueError("Not around zero")

            # Set max
            if symmetric_method == "mean": vmax = 0.5 * sum([abs(vmin), abs(vmax)]) # mean
            elif symmetric_method == "max": vmax = max([abs(vmin), abs(vmax)])
            elif symmetric_method == "min": vmax = min([abs(vmin), abs(vmax)])
            else: raise ValueError("Invalid symmetric method")

            # Set min
            vmin = - vmax

    # List or tuple of 2 values (min and max)
    elif isinstance(interval, list) or isinstance(interval, tuple):

        vmin, vmax = interval

        if soft_min: vmin = max(vmin / soft_min_scaling, normalize_min)
        if soft_max: vmax = min(vmax * soft_max_scaling, normalize_max)

    # String -> parse
    elif isinstance(interval, basestring):

        #print(interval)
        from ...core.tools import parsing
        try: vmin, vmax = parsing.real_pair(interval)
        except ValueError: raise ValueError("Cannot interpret the interval")

        if soft_min: vmin = max(vmin/soft_min_scaling, normalize_min)
        if soft_max: vmax = min(vmax*soft_max_scaling, normalize_max)

    # Other
    else: raise ValueError("Invalid option for 'interval'")  # INVALID

    # Adjust to user-defined vmin and vmax
    if zmin is not None:
        if soft_zmin: vmin = max(vmin, zmin)
        else: vmin = zmin
    if zmax is not None:
        if soft_zmax: vmax = min(vmax, zmax)
        else: vmax = zmax

    # Return the interval
    return vmin, vmax

# -----------------------------------------------------------------

def plot_box(box, title=None, path=None, format=None, scale="log", interval="pts", cmap="viridis", colorbar=False,
             around_zero=False, symmetric=False, normalize_in=None, scale_parameter=None, show_axes=True,
             transparent=False, soft_min=False, soft_max=False, soft_min_scaling=1., soft_max_scaling=1.,
             region=None, regions=None, axes=None, xsize=7, ysize=7, interpolation="nearest", alpha=1, aspect="equal",
             symmetric_method="mean", check_around_zero=True, background_color=None, plot=None):

    """
    This function ...
    :param box:
    :param title:
    :param path:
    :param format:
    :param scale:
    :param interval:
    :param cmap:
    :param colorbar:
    :param around_zero:
    :param symmetric:
    :param normalize_in:
    :param scale_parameter:
    :param show_axes:
    :param transparent:
    :param soft_min:
    :param soft_max:
    :param soft_min_scaling:
    :param soft_max_scaling:
    :param region:
    :param regions:
    :param axes:
    :param xsize:
    :param ysize:
    :param interpolation:
    :param alpha:
    :param aspect:
    :param symmetric_method:
    :param check_around_zero:
    :param background_color:
    :param plot:
    :return:
    """

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Get the data
    if isinstance(box, np.ndarray): data = box
    else: data = box.data

    # Get coordinate system
    if hasattr(box, "wcs"): wcs = box.wcs
    else: wcs = None

    # Get dimension
    nxpix = data.shape[1]
    nypix = data.shape[0]

    # Get interval
    vmin, vmax = get_vmin_vmax(data, interval=interval, around_zero=around_zero, symmetric=symmetric,
                               normalize_in=normalize_in, soft_min=soft_min, soft_max=soft_max,
                               soft_min_scaling=soft_min_scaling, soft_max_scaling=soft_max_scaling,
                               symmetric_method=symmetric_method, check_around_zero=check_around_zero, wcs=wcs)
    #print(vmin, vmax, scale)
    output.vmin = vmin
    output.vmax = vmax

    # Get the normalization
    norm = get_normalization(scale, vmin, vmax, data=data, scale_parameter=scale_parameter)
    output.norm = norm

    # Create figure if necessary, get the axes
    if plot is not None: axes = plot.axes
    only_axes = False
    if axes is None:
        fig = plt.figure(figsize=(xsize,ysize))
        plt.xlim(0, nxpix - 1)
        plt.ylim(0, nypix - 1)
        axes = plt.gca()
        output.figure = fig
    else: only_axes = True
    output.axes = axes

    # Show the data
    extent = None
    image = axes.imshow(data, origin="lower", interpolation=interpolation, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap,
                alpha=alpha, aspect=aspect, extent=extent)
    output.image = image

    # Set background color
    if background_color is not None: axes.set_facecolor(background_color)

    # Add region
    if region is not None:

        from ..region.composite import CompositeRegion
        if isinstance(region, CompositeRegion):
            for patch in region.to_mpl_patches(): axes.add_patch(patch)
        else:
            # Add patch
            patch = region.to_mpl_patch()
            axes.add_patch(patch)

    # Add multiple regions
    if regions is not None:
        for patch in regions.to_mpl_patches(): axes.add_patch(patch)

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Show axis?
        if not show_axes: plt.axis('off')

        # Add title
        if title is not None: plt.title(title)

        # Add colorbar?
        if colorbar: plt.colorbar(image)

        # Show or save as file
        if path is None: plt.show()
        else: plt.savefig(path, format=format, transparent=transparent)

        # Close the figure
        plt.close(output.figure)

    # Return
    return output

# -----------------------------------------------------------------

def plot_residual_map(frame, is_percentual=False, title=None, plot=None, path=None):

    """
    Thisf unction ...
    :param frame:
    :param is_percentual:
    :param title:
    :param plot:
    :param path:
    :return:
    """

    # Plot
    cmap = 'RdBu'

    # Set interval
    if is_percentual: interval = (-100,100,)
    else: interval = (-1,1,)

    # Plot
    plot_map(frame, interval=interval, scale="linear", colorbar=True, cmap=cmap, path=path, title=title, plot=plot)

# -----------------------------------------------------------------

def plot_map(frame, interval="pts", scale="linear", colorbar=True, cmap="inferno", contours=False, ncontours=5,
             contours_color="white", path=None, background_color=None, title=None, plot=None, aspect="equal"):

    """
    This function ...
    :param frame:
    :param interval:
    :param scale:
    :param colorbar:
    :param cmap:
    :param contours:
    :param ncontours:
    :param contours_color:
    :param path:
    :param background_color:
    :param title:
    :param plot:
    :param aspect:
    :return:
    """

    # With contours
    if contours:

        # Plot with contours
        return plot_frame_contours(frame, interval=interval, scale=scale, colorbar=colorbar,
                                    data_cmap=cmap, plot_data=True, nlevels=ncontours,
                                    single_colour=contours_color, path=path, background_color=background_color, title=title,
                                    plot=plot, aspect=aspect)

    # No contours
    else: return plot_frame(frame, interval=interval, scale=scale, colorbar=colorbar, cmap=cmap, path=path,
                     background_color=background_color, title=title, plot=plot, aspect=aspect)

# -----------------------------------------------------------------

def plot_map_offset(frame, center, radius, offset_step, interval="pts", scale="linear", colorbar=True, cmap="inferno", contours=False, ncontours=5,
                    contours_color="white", path=None, background_color=None, title=None, plot=None, axes_label="Distance from center", aspect="equal"):

    """
    This function ...
    :param frame:
    :param center: CENTER IN PIXEL COORDINATES
    :param radius: RADIUS IN LENGTH
    :param offset_step:
    :param interval:
    :param scale:
    :param colorbar:
    :param cmap:
    :param contours:
    :param ncontours:
    :param contours_color:
    :param path:
    :param background_color:
    :param title:
    :param plot:
    :param axes_label:
    :param aspect:
    :return:
    """

    from ...core.tools import numbers
    from ...core.basics.range import RealRange

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Get physical pixelscale
    pixelscale = frame.average_physical_pixelscale

    # GET OFFSET STEP AND UNITS
    offset = offset_step.value
    offset_unit = offset_step.unit
    offset_unit_string = str(offset_unit)

    # Get physical radius, round up to
    radius = numbers.round_up_to_base(radius.to(offset_unit).value, base=offset) * offset_unit # round up to 5 kpc
    radius_pixels = numbers.round_up_to_int(radius.value / pixelscale.to(offset_unit).value)

    # Make cutout
    cutout = frame.cutout_around(center, radius_pixels, as_frame=True)
    # print(cutout.shape) # must preferentially be odd so that center is still in center

    # Get new pixel center
    center_pix = cutout.pixel_center
    output.x_min = cutout.get_meta("x_min")
    output.x_max = cutout.get_meta("x_max")
    output.y_min = cutout.get_meta("y_min")
    output.y_max = cutout.get_meta("y_max")
    output.center_pix = center_pix

    # Plot
    output.update(plot_map(cutout, interval=interval, scale=scale, colorbar=colorbar, cmap=cmap, contours=contours, ncontours=ncontours,
                    contours_color=contours_color, path=path, background_color=background_color, title=title, plot=plot, aspect=aspect))

    # Set axes
    plot.set_xlabel(axes_label + " [" + offset_unit_string + "]")
    plot.set_ylabel(axes_label + " [" + offset_unit_string + "]")

    # Determine ticks
    pixelscale_in_offset_unit = pixelscale.to(offset_unit).value

    # Get offsets
    # offsets = [-15, -10, -5, 0, 5, 10, 15]  # in kpc
    radius_range = RealRange(-radius.value, radius.value, inclusive=True)
    offsets = radius_range.linear_step(offset, symmetric=True, center=0.)
    output.offsets = offsets

    # labels = [str(offset) for offset in kpc_offsets]
    xticks = [center_pix.x + float(offset) / pixelscale_in_offset_unit for offset in offsets]
    yticks = [center_pix.y + float(offset) / pixelscale_in_offset_unit for offset in offsets]

    # Set ticks
    plot._plot.set_xticks(xticks)
    plot._plot.set_xticklabels(offsets)
    plot._plot.set_yticks(yticks)
    plot._plot.set_yticklabels(offsets)

    # Return the output
    return output

# -----------------------------------------------------------------

def plot_map_centered(frame, radius, offset, interval=None, cmap="viridis", plot=None, plot_radii=None, colorbar=False, background_color=None):

    """
    This function ...
    :param frame:
    :param radius:
    :param offset:
    :param interval:
    :param cmap:
    :param plot:
    :param plot_radii:
    :param colorbar:
    :param background_color:
    :return:
    """

    from matplotlib.patches import Circle

    # Assume center pixel is the center of the map
    center_pix = frame.pixel_center

    # Plot
    output = plot_map_offset(frame, center_pix, radius, offset, interval=interval, colorbar=colorbar, cmap=cmap, plot=plot, background_color=background_color)

    # Plot radii?
    if plot_radii is not None:

        # Loop over the radii
        for rad in plot_radii:

            # Create patch
            inner_radius = (rad / frame.average_physical_pixelscale).to("").value # in pixels
            inner_circle = Circle(output.center_pix.cartesian, inner_radius, fill=False, edgecolor="white")

            # Add circles
            output.axes.add_patch(inner_circle)

    # Return
    return output

# -----------------------------------------------------------------

def plot_frame_contours(frame, **kwargs):

    """
    This function ...
    :param frame:
    :param kwargs:
    :return:
    """

    # Create prepared frame
    frame, kwargs = create_prepared_frame(frame, **kwargs)

    # Plot
    return plot_contours(frame, **kwargs)

# -----------------------------------------------------------------

def plot_contours(box, nlevels=20, path=None, x_label="x", y_label="y", line_width=1, font_size=16, title=None,
                  format=None, cmap="jet", single_colour=None, labels=False, show_axes=True, transparent=False,
                  plot_data=False, axes=None, xsize=7, ysize=7, colorbar=False, scale="linear",

                  # For plotting data
                  interval="pts", data_cmap="viridis", around_zero=False, symmetric=False, soft_min=False, soft_max=False,
                  soft_min_scaling=1., soft_max_scaling=1., interpolation="nearest", alpha=1,
                  aspect="equal", symmetric_method="mean",
                  check_around_zero=True, background_color=None, plot=None
                  ):

    """
    This function ...
    :param box:
    :param nlevels:
    :param path:
    :param x_label:
    :param y_label:
    :param line_width:
    :param font_size:
    :param title:
    :param format:
    :param cmap:
    :param single_colour:
    :param labels:
    :param show_axes:
    :param transparent:
    :param plot_data:
    :param axes:
    :param xsize:
    :param ysize:
    :param colorbar:
    :param scale:
    :param interval:
    :param data_cmap:
    :param around_zero:
    :param symmetric:
    :param soft_min:
    :param soft_max:
    :param soft_min_scaling:
    :param soft_max_scaling:
    :param interpolation:
    :param alpha:
    :param aspect:
    :param symmetric_method:
    :param check_around_zero:
    :param background_color:
    :param plot:
    :return:
    """

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Get the data
    if isinstance(box, np.ndarray): data = box
    else: data = box.data

    # Get shape
    nxpix = data.shape[1]
    nypix = data.shape[0]

    # Create figure if necessary, get the axes
    if plot is not None: axes = plot.axes
    only_axes = False
    if axes is None:
        plt.figure(figsize=(xsize, ysize))
        plt.xlim(0, nxpix - 1)
        plt.ylim(0, nypix - 1)
        axes = plt.gca()
    else: only_axes = True
    output.axes = axes

    # Get interval
    vmin, vmax = get_vmin_vmax(data, interval=interval, around_zero=around_zero, symmetric=symmetric,
                               soft_min=soft_min, soft_max=soft_max, soft_min_scaling=soft_min_scaling,
                               soft_max_scaling=soft_max_scaling, symmetric_method=symmetric_method,
                               check_around_zero=check_around_zero)
    output.vmin = vmin
    output.vmax = vmax

    # Get the normalization
    norm = get_normalization(scale, vmin, vmax, data=data)
    output.norm = norm

    # Also plot the image data underneath
    if plot_data:

        # Set cmap
        if data_cmap is None: data_cmap = cmap

        # Show the data
        extent = None
        image = axes.imshow(data, origin="lower", interpolation=interpolation, vmin=vmin, vmax=vmax, norm=norm,
                            cmap=data_cmap, alpha=alpha, aspect=aspect, extent=extent)
        output.image = image

    # Set to None
    else: vmin = vmax = norm = image = None

    # Set background color
    if background_color is not None: axes.set_facecolor(background_color)

    # Define X and Y labels
    x = np.arange(nxpix)
    y = np.arange(nypix)

    # Square figure
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    axes.set_aspect('equal')

    # Contour and labels
    if single_colour is not None:
        colors = single_colour
        cmap = None
    else: colors = None
    cs = axes.contour(x, y, data, nlevels, colors=colors, cmap=cmap, linewidths=line_width, norm=norm)
    if labels: plt.clabel(cs, fontsize=font_size)

    # Axes labels
    axes.set_xlabel(x_label, fontsize=font_size)
    axes.set_ylabel(y_label, fontsize=font_size)

    # Make tick lines thicker
    for l in axes.get_xticklines(): l.set_markeredgewidth(line_width)
    for l in axes.get_yticklines(): l.set_markeredgewidth(line_width)

    # Make axis label larger
    for tick in axes.xaxis.get_major_ticks(): tick.label.set_fontsize(font_size)
    for tick in axes.yaxis.get_major_ticks(): tick.label.set_fontsize(font_size)

    # Make figure box thicker
    for s in axes.spines.values(): s.set_linewidth(line_width)

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Show axis?
        if not show_axes: plt.axis('off')

        # Add title
        if title is not None: plt.title(title)

        # Add colorbar?
        if colorbar and image is not None: plt.colorbar(image)

        # Show or save as file
        if path is None: plt.show()
        else: plt.savefig(path, format=format, transparent=transparent)

        # Close the figure
        plt.close()

    # Return
    return output

# -----------------------------------------------------------------

def plot_filled_frame_contours(frame, **kwargs):

    """
    This function ...
    :param frame:
    :param kwargs:
    :return:
    """

    # Create prepared frame
    frame, kwargs = create_prepared_frame(frame, **kwargs)

    # Plot
    plot_filled_contours(frame, **kwargs)

# -----------------------------------------------------------------

def plot_filled_contours(box, nlevels=20, title=None, path=None, format=None, cmap="jet", x_label="x", y_label="y",
                         show_axes=True, transparent=False):

    """
    This function ...
    :param box:
    :param nlevels:
    :param title:
    :param path:
    :param format:
    :param cmap:
    :param x_label:
    :param y_label:
    :param show_axes:
    :param transparent:
    :return:
    """

    #from matplotlib.colors import LogNorm

    # Get the data
    if isinstance(box, np.ndarray): data = box
    else: data = box.data
    z = data

    # Define x and y labels
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])

    # Setup figure and colorbar
    fig = plt.figure(figsize=(6, 6.5))

    # Define custom grid for image
    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(1, 1),
                     label_mode='L',
                     cbar_pad=0.05,
                     cbar_mode='single',
                     cbar_location='top')

    # Set axis labels
    grid[0].set_xlabel(x_label)
    grid[0].set_ylabel(y_label)

    # Display image contours and colorbar
    #im = grid[0].contourf(x, y, z, levels=levels, antialiased=False, cmap=cmap)
    im = grid[0].contourf(x, y, z, nlevels, antialiased=False, cmap=cmap)
    cbar = grid.cbar_axes[0].colorbar(im)

    # Set tick labels to colorbar (even numbers)
    #cbar.ax.set_xticks(levels[::2])

    # Show axes?
    if not show_axes: plt.axis("off")

    # Add title
    if title is not None: plt.title(title)

    # Show or save
    if path is None: plt.show()
    else: plt.savefig(path, format=format, transparent=transparent)

    # Close
    plt.close()

# -----------------------------------------------------------------

def plot_radial_profile(box, center, angle, ratio, nbins=20, path=None, title=None, measure="mean", max_radius=None, transparent=False):

    """
    This function ...
    :param box:
    :param center:
    :param angle:
    :param ratio:
    :param nbins:
    :param path:
    :param title:
    :param measure: sum, min, max, mean, median
    :param max_radius:
    :param transparent:
    :return:
    """

    # Get the data
    if isinstance(box, np.ndarray): data = box
    else: data = box.data

    from ..dist_ellipse import distance_ellipse
    from ...core.basics.range import RealRange
    from ..basics.coordinate import PixelCoordinate, SkyCoordinate
    from ..core.mask import intersection

    # Convert center to pixel coordinates
    if isinstance(center, PixelCoordinate): center_pix = center
    elif isinstance(center, SkyCoordinate):
        if not hasattr(box, "wcs"): raise ValueError("Passed data must have a coordinate system if center is given as sky coordinate")
        center_pix = center.to_pixel(box.wcs)
    else: raise ValueError("center must be pixel or sky coordinate")

    # Create distance ellipse frame
    distance_frame = distance_ellipse(data.shape, center_pix, ratio, angle)

    radius_list = []
    value_list = []
    #nmasked_list = []

    # Determine minimum distance
    min_distance = np.min(distance_frame)

    # Determine maximum distance
    if max_radius is not None:
        from astropy.units import Quantity
        if isinstance(max_radius, Quantity):
            if max_radius.unit.physical_type != "angle": raise ValueError("Unknown max radius unit")
            if not (hasattr(box, "pixelscale") and box.pixelscale is not None): raise ValueError("Passed data must have its pixelscale defined")
            max_radius = (max_radius / box.average_pixelscale).to("").value
        max_distance = max_radius
    else: max_distance = np.max(distance_frame)

    # Determine the step size
    step = (max_distance - min_distance) / float(nbins)

    # Debugging
    log.debug("Step size for plotting the radial profile is " + str(step) + " pixels")

    # Set the first range
    radius_range = RealRange(min_distance, min_distance + step)

    # Loop, shifting ranges of radius
    while True:

        # Check the range
        if radius_range.min > max_distance: break

        # Get the average radius
        radius_center = radius_range.center

        above_min_mask = distance_frame >= radius_range.min
        below_max_mask = distance_frame < radius_range.max

        # Create mask of pixels that are in the range
        range_mask = intersection(above_min_mask, below_max_mask)

        # Calculate the mean value in the range
        values = data[range_mask]
        if measure == "mean": value = np.nanmean(values)
        elif measure == "median": value = np.nanmedian(values)
        elif measure == "sum": value = np.nansum(values)
        elif measure == "min": value = np.nanmin(values)
        elif measure == "max": value = np.nanmax(values)
        else: raise ValueError("Invalid value for 'measure': " + measure)

        # Make a mask of all the pixels below the center radius
        #below_mask = distance_frame < radius_center

        # Add point
        radius_list.append(radius_center)
        value_list.append(value)

        # Shift the range
        radius_range += step

    # Create the plot
    plot_xy(radius_list, value_list, title=title, transparent=transparent, path=path)

# -----------------------------------------------------------------
# GET XY DATA
# -----------------------------------------------------------------

import operator
ops = { "+": operator.add, "-": operator.sub } # etc.

# -----------------------------------------------------------------

def has_add(name):
    return "+" in name

# -----------------------------------------------------------------

def has_sub(name):
    return "-" in name

# -----------------------------------------------------------------

def has_div(name):
    return "/" in name

# -----------------------------------------------------------------

def has_mul(name):
    return "*" in name

# -----------------------------------------------------------------

def has_operators(name):
    return has_add(name) or has_sub(name) or has_div(name) or has_mul(name)

# -----------------------------------------------------------------

def get_aux(curve, colname):

    """
    This function ...
    :param curve:
    :param colname:
    :return:
    """

    from ...core.units.unit import divide_units

    # Has division symbol
    if has_div(colname):

        colname_a, colname_b = colname.split("/")
        colname_a = colname_a.strip()
        colname_b = colname_b.strip()

        # Get units
        unit_a = curve.get_unit(colname_a)
        unit_b = curve.get_unit(colname_b)

        # Get data
        aux_a = curve.get_column_array(colname_a, unit=unit_a)
        aux_b = curve.get_column_array(colname_b, unit=unit_b)

        # Get unit and calculate data
        unit = divide_units(unit_a, unit_b)
        aux = aux_a / aux_b

        # Return
        return aux, unit

    # No symbol: pure column name
    else:

        aux = curve.get_array(colname)
        aux_unit = curve.get_unit(colname)
        return aux, aux_unit

# -----------------------------------------------------------------

def get_x(curve, colname=None, return_label=False, label=None):

    """
    This function ...
    :param curve:
    :param colname:
    :param return_label:
    :param label:
    :return:
    """

    from ...core.units.unit import divide_units

    # Check whether operators
    if colname is not None and has_div(colname):

        colname_a, colname_b = colname.split("/")
        colname_a = colname_a.strip()
        colname_b = colname_b.strip()

        # Get units
        unit_a = curve.get_unit(colname_a)
        unit_b = curve.get_unit(colname_b)

        # Get data
        x_a = curve.get_column_array(colname_a, unit=unit_a)
        x_b = curve.get_column_array(colname_b, unit=unit_b)

        # Get unit and calculate data
        unit = divide_units(unit_a, unit_b)
        x = x_a / x_b

        # Set label
        if label is not None: x_label = label + " [" + str(unit) + "]"
        else: x_label = colname + " [" + str(unit) + "]"

        # Return
        if return_label: return x, x_label
        else: return x

    # Just one column
    else: return get_x_column(curve, colname=colname, return_label=return_label, label=label)

# -----------------------------------------------------------------

def get_x_column(curve, colname=None, return_label=False, label=None):

    """
    This function ...
    :param curve:
    :param colname:
    :param return_label:
    :param label:
    :return:
    """

    # Get the units
    if colname is not None: x_unit = curve.get_unit(colname)
    else: x_unit = curve.x_unit

    # Get the data
    if colname is not None: x = curve.get_column_array(colname, unit=x_unit)
    else: x = curve.get_x(asarray=True, unit=x_unit)

    # Get the labels
    if label is not None: x_label = label
    elif colname is not None: x_label = colname
    else: x_label = curve.x_name
    if x_unit is not None: x_label += " [" + str(x_unit) + "]"

    # Return
    if return_label: return x, x_label.encode('utf-8')
    else: return x

# -----------------------------------------------------------------

def get_y(curve, colname=None, return_label=False, label=None):

    """
    This function ...
    :param curve:
    :param colname:
    :param return_label:
    :param label:
    :return:
    """

    from ...core.units.unit import divide_units

    # Check whether operators
    if colname is not None and has_div(colname):

        colname_a, colname_b = colname.split("/")
        colname_a = colname_a.strip()
        colname_b = colname_b.strip()

        # Get units
        unit_a = curve.get_unit(colname_a)
        unit_b = curve.get_unit(colname_b)

        # Get data
        y_a = curve.get_column_array(colname_a, unit=unit_a)
        y_b = curve.get_column_array(colname_b, unit=unit_b)

        # Get unit and data
        unit = divide_units(unit_a, unit_b)
        y = y_a / y_b

        # Set label
        if label is not None: y_label = label + " [" + str(unit) + "]"
        else: y_label = colname + " [" + str(unit) + "]"

        # Return
        if return_label: return y, y_label
        else: return y

    # Just one column
    else: return get_y_column(curve, colname=colname, return_label=return_label, label=label)

# -----------------------------------------------------------------

def get_y_column(curve, colname=None, return_label=False, label=None):

    """
    This function ...
    :param curve:
    :param colname:
    :param return_label:
    :param label:
    :return:
    """

    # Get the units
    if colname is not None: y_unit = curve.get_unit(colname)
    else: y_unit = curve.y_unit

    # Get the data
    if colname is not None: y = curve.get_column_array(colname, unit=y_unit)
    else: y = curve.get_y(asarray=True, unit=y_unit)

    # Get the labels
    if label is not None: y_label = label
    elif colname is not None: y_label = colname
    else: y_label = curve.y_name
    if y_unit is not None: y_label += " [" + str(y_unit) + "]"

    # Return
    if return_label: return y, y_label.encode('utf-8')
    else: return y

# -----------------------------------------------------------------

def get_xy(curve, return_labels=False, x_colname=None, y_colname=None, x_label=None, y_label=None):

    """
    This function ...
    :param curve:
    :param return_labels:
    :param x_colname:
    :param y_colname:
    :param x_label:
    :param y_label:
    :return:
    """

    # Get x
    x, x_label = get_x(curve, colname=x_colname, return_label=True, label=x_label)

    # Get y
    y, y_label = get_y(curve, colname=y_colname, return_label=True, label=y_label)

    # Return
    if return_labels: return x, y, x_label, y_label
    else: return x, y

# -----------------------------------------------------------------

def get_multiple_xy(curves, return_labels=False, return_units=False):

    """
    Thisf unction ...
    :param curves:
    :param return_labels:
    :param return_units:
    :return:
    """

    # Initialize dictionaries for x and y data
    x = OrderedDict()
    y = OrderedDict()

    # Get curves and names
    ncurves = len(curves)
    if types.is_dictionary(curves):

        names = curves.keys()
        curves = curves.values()

    elif types.is_sequence(curves):

        # Determine curve names
        y_labels = [curve.y_name for curve in curves]
        if sequences.all_different(y_labels): names = y_labels
        else: names = ["curve" + str(index) for index in range(ncurves)]

    # Invalid
    else: raise ValueError("Invalid argument")

    # Get the labels
    x_labels = [curve.x_name for curve in curves]
    y_labels = [curve.y_name for curve in curves]
    x_label = sequences.get_single(x_labels, method="common")
    y_label = sequences.get_single(y_labels, method="common")

    # Get the units
    x_units = [curve.x_unit for curve in curves]
    y_units = [curve.y_unit for curve in curves]
    x_unit = sequences.get_single(x_units, method="first_not_none")
    y_unit = sequences.get_single(y_units, method="first_not_none")

    # Loop over the curves
    for name, curve in zip(names, curves):

        # Get the data
        x_curve = curve.get_x(asarray=True, unit=x_unit)
        y_curve = curve.get_y(asarray=True, unit=y_unit)

        # Set to dictionary
        x[name] = x_curve
        y[name] = y_curve

    # Add units to labels
    if x_unit is not None: x_label += " [" + str(x_unit) + "]"
    if y_unit is not None: y_label += " [" + str(y_unit) + "]"

    # Return
    if return_labels:
        if return_units: return x, y, x_label, y_label, x_unit, y_unit
        else: return x, y, x_label, y_label
    else:
        if return_units: return x, y, x_unit, y_unit
        else: return x, y

# -----------------------------------------------------------------
# PLOTTING CURVES
# -----------------------------------------------------------------

def plot_curve(curve, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
               xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False,
               ynonzero=False, x_label=None, y_label=None, plot=None, vlines=None, hlines=None, color=None,
               x_color=None, y_color=None, vlinestyle="solid", hlinestyle="solid", vlinecolor=None, hlinecolor=None):

    """
    This function ...
    :param curve:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param x_label:
    :param y_label:
    :param plot:
    :param vlines:
    :param hlines:
    :param color:
    :param x_color:
    :param y_color:
    :param vlinestyle:
    :param hlinestyle:
    :param vlinecolor:
    :param hlinecolor:
    :return:
    """

    # Get x, y and labels
    x, y, _x_label, _y_label = get_xy(curve, return_labels=True)
    if x_label is None: x_label = _x_label
    if y_label is None: y_label = _y_label

    # Convert x limits if necessary
    if xlimits is not None:
        xlimits = list(xlimits) # for when tuple
        if hasattr(xlimits[0], "unit"): xlimits[0] = xlimits[0].to(curve.x_unit).value
        if hasattr(xlimits[1], "unit"): xlimits[1] = xlimits[1].to(curve.x_unit).value

    # Convert y limits if necessary
    if ylimits is not None:
        ylimits = list(ylimits) # for when tuple
        if hasattr(ylimits[0], "unit"): ylimits[0] = ylimits[0].to(curve.y_unit).value
        if hasattr(ylimits[1], "unit"): ylimits[1] = ylimits[1].to(curve.y_unit).value

    # Convert vertical lines if necessary
    if vlines is not None:
        vlines = [value.to(curve.x_unit).value if hasattr(value, "unit") else value for value in vlines]

    # Convert horizontal lines if necessary
    if hlines is not None:
        hlines = [value.to(curve.y_unit).value if hasattr(value, "unit") else value for value in hlines]

    # Plot
    plot_xy(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog, connect=True,
            xlimits=xlimits, ylimits=ylimits, xpositive=xpositive, ypositive=ypositive, xnonnegative=xnonnegative,
            ynonnegative=ynonnegative, xnonzero=xnonzero, ynonzero=ynonzero, plot=plot, vlines=vlines, hlines=hlines,
            color=color, x_color=x_color, y_color=y_color, vlinestyle=vlinestyle, hlinestyle=hlinestyle,
            vlinecolor=vlinecolor, hlinecolor=hlinecolor)

# -----------------------------------------------------------------

def plot_curves(curves, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
                xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False,
                ynonzero=False, x_label=None, y_label=None, show=None):

    """
    This function ...
    :param curves:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param x_label:
    :param y_label:
    :param show:
    :return:
    """

    # Get data
    x, y, _x_label, _y_label, x_unit, y_unit = get_multiple_xy(curves, return_labels=True, return_units=True)
    if x_label is None: x_label = _x_label
    if y_label is None: y_label = _y_label

    # Convert x limits if necessary
    if xlimits is not None:
        xlimits = list(xlimits) # for when tuple
        if hasattr(xlimits[0], "unit"): xlimits[0] = xlimits[0].to(x_unit).value
        if hasattr(xlimits[1], "unit"): xlimits[1] = xlimits[1].to(x_unit).value

    # Convert y limits if necessary
    if ylimits is not None:
        ylimits = list(ylimits)
        if hasattr(ylimits[0], "unit"): ylimits[0] = ylimits[0].to(y_unit).value
        if hasattr(ylimits[1], "unit"): ylimits[1] = ylimits[1].to(y_unit).value

    # Plot
    return plot_xy(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog, connect=True,
            xlimits=xlimits, ylimits=ylimits, xpositive=xpositive, ypositive=ypositive, xnonnegative=xnonnegative,
            ynonnegative=ynonnegative, xnonzero=xnonzero, ynonzero=ynonzero, show=show)

# -----------------------------------------------------------------
# PLOTTING SCATTER
# -----------------------------------------------------------------

def plot_scatter(scatter, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None, density=False,
                 xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False,
                 ynonzero=False, size=None, show=None):

    """
    This function ...
    :param scatter:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param density:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param size:
    :param show:
    :return:
    """

    # Get x, y and labels
    x, y, x_label, y_label = get_xy(scatter, return_labels=True)

    # Plot
    plot_xy(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog,
            connect=False, density=density, xlimits=xlimits, ylimits=ylimits, xpositive=xpositive, ypositive=ypositive,
            xnonnegative=xnonnegative, ynonnegative=ynonnegative, xnonzero=xnonzero, ynonzero=ynonzero, size=size, show=show)

# -----------------------------------------------------------------

def plot_scatter_density(scatter, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None, show=None,
                         colormaps=False, axes=None, plot=None, color=None, dpi=None, aux_colname=None, aux=None,
                         aux_name=None, aux_unit=None, aux_log=False, aux_limits=None, valid_points=None, x_colname=None,
                         y_colname=None, legend_location=None, density_log=False, aux_density=False, cmap=None,
                         x_label=None, y_label=None, label=None, remove=None, log_enhance=50.):

    """
    This function ...
    :param scatter:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param show:
    :param colormaps:
    :param axes:
    :param plot:
    :param color:
    :param dpi:
    :param aux_colname:
    :param aux:
    :param aux_name:
    :param aux_unit:
    :param aux_log:
    :param aux_limits:
    :param valid_points:
    :param x_colname:
    :param y_colname:
    :param legend_location:
    :param density_log:
    :param aux_density:
    :param cmap:
    :param x_label:
    :param y_label:
    :param label:
    :param remove:
    :param log_enhance:
    :return:
    """

    # Get x, y and labels
    #print(x_colname, y_colname)
    x, y, x_label, y_label = get_xy(scatter, return_labels=True, x_colname=x_colname, y_colname=y_colname, x_label=x_label, y_label=y_label)
    #print(x_label, y_label)

    # Get aux values
    if aux_colname is not None:
        if aux is not None: raise ValueError("Cannot pass both auxilary column name and auxilary values")
        aux, aux_unit = get_aux(scatter, aux_colname)
        if aux_name is None: aux_name = aux_colname

    # Filter?
    if valid_points is not None:
        x = x[valid_points]
        y = y[valid_points]
        if aux is not None: aux = aux[valid_points]

    # Remove?
    if remove is not None:

        remove_x = remove[0]
        remove_y = remove[1]

        # Initialize remove mask
        remove = np.ones_like(x, dtype=bool)

        if remove_x is not None:

            if remove_x.lower is not None:
                remove *= x > remove_x.lower

            if remove_x.upper is not None:
                remove *= x < remove_x.upper

        if remove_y is not None:

            if remove_y.lower is not None:
                remove *= y > remove_y.lower

            if remove_y.upper is not None:
                remove *= y < remove_y.upper

        # Add condition
        # mask *= np.logical_not(remove)
        mask = np.logical_not(remove)
        x = x[mask]
        y = y[mask]
        if aux is not None: aux = aux[mask]

    if aux is not None:
        aux_min = np.nanmin(aux)
        aux_max = np.nanmax(aux)
        print("aux min:", aux_min)
        print("aux max:", aux_max)

    # Plot
    return plot_xy_scatter_density(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog,
                                     xlimits=xlimits, ylimits=ylimits, show=show, colormaps=colormaps, axes=axes, plot=plot,
                                     color=color, dpi=dpi, aux=aux, aux_name=aux_name, aux_unit=aux_unit, aux_log=aux_log,
                                     aux_limits=aux_limits, legend_location=legend_location, density_log=density_log, aux_density=aux_density,
                                     cmap=cmap, label=label, log_enhance=log_enhance)

# -----------------------------------------------------------------

def plot_scatters(scatters, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None, density=False,
                  xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False,
                  ynonzero=False, size=None, show=None):

    """
    This function ...
    :param scatters:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param density:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param size:
    :param show:
    :return:
    """

    # Get data
    x, y, x_label, y_label = get_multiple_xy(scatters, return_labels=True)

    # Plot
    return plot_xy(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog,
            connect=False, density=density, xlimits=xlimits, ylimits=ylimits, xpositive=xpositive, ypositive=ypositive,
            xnonnegative=xnonnegative, ynonnegative=ynonnegative, xnonzero=xnonzero, ynonzero=ynonzero, size=size, show=show)

# -----------------------------------------------------------------

def plot_scatters_density(scatters, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
                            show=None, colormaps=False, colors=None, axes=None, plot=None, dpi=None,
                            legend_location=None, density_log=False):

    """
    This function ...
    :param scatters:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param show:
    :param colormaps:
    :param colors:
    :param axes:
    :param plot:
    :param dpi:
    :param legend_location:
    :param density_log:
    :return:
    """

    # Get data
    x, y, x_label, y_label = get_multiple_xy(scatters, return_labels=True)

    # Plot
    return plot_xy_scatter_density(x, y, title=title, path=path, x_label=x_label, y_label=y_label, xlog=xlog, ylog=ylog,
                                   xlimits=xlimits, ylimits=ylimits, show=show, colormaps=colormaps, colors=colors,
                                   axes=axes, plot=plot, dpi=dpi, legend_location=legend_location, density_log=density_log)

# -----------------------------------------------------------------
# PLOTTING DENSITY
# -----------------------------------------------------------------

def plot_density(points, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
                 nbins=200, contours=False, rug=False, xpositive=False, ypositive=False,
                xnonnegative=False, ynonnegative=False, xnonzero=False, ynonzero=False):

    """
    This function ...
    :param points:
    :param title:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param nbins:
    :param contours:
    :param rug:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :return:
    """

    # Get x, y and labels
    x, y, x_label, y_label = get_xy(points, return_labels=True)

    # Plot
    plot_xy_density(x, y, title=title, nbins=nbins, contours=contours, path=path, rug=rug,
                    xlog=xlog, ylog=ylog, xlimits=xlimits, ylimits=ylimits, xpositive=xpositive, ypositive=ypositive,
                    xnonnegative=xnonnegative, ynonnegative=ynonnegative, xnonzero=xnonzero, ynonzero=ynonzero)

# -----------------------------------------------------------------

def plot_densities(points, title=None, path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
                   nbins=200, contours=False, rug=False, method=None, show=None):

    """
    This function ...
    :param points: 
    :param title: 
    :param path: 
    :param xlog:
    :param ylog:
    :param xlimits: 
    :param ylimits: 
    :param nbins:
    :param contours:
    :param rug:
    :param method:
    :param show:
    :return: 
    """

    # Get data
    x, y, x_label, y_label = get_multiple_xy(points, return_labels=True)

    # Plot
    return plot_xy_density(x, y, title=title, nbins=nbins, contours=contours, path=path, rug=rug,
                    xlog=xlog, ylog=ylog, xlimits=xlimits, ylimits=ylimits, method=method, show=show)

# -----------------------------------------------------------------

def plot_xy(x, y, title=None, path=None, transparent=False, x_label=None, y_label=None, xlog=False,
            ylog=False, vlines=None, hlines=None, legend=True, xlimits=None, ylimits=None, connect=True,
            density=False, xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False,
            ynonzero=False, size=None, axes=None, plot=None, show=None, color=None, x_color=None, y_color=None,
            vlinestyle="solid", hlinestyle="solid", vlinecolor=None, hlinecolor=None):

    """
    Low-level function, only scalar values (no units)
    x and y as dictionaries or lists
    :param x:
    :param y:
    :param title:
    :param path:
    :param transparent:
    :param x_label:
    :param y_label:
    :param xlog:
    :param ylog:
    :param vlines:
    :param hlines:
    :param legend:
    :param xlimits:
    :param ylimits:
    :param connect:
    :param density:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param size:
    :param axes:
    :param plot:
    :param show:
    :param color:
    :param x_color:
    :param y_color:
    :param vlinestyle:
    :param hlinestyle:
    :param vlinecolor:
    :param hlinecolor:
    :return:
    """

    from pts.core.basics.map import Map
    from pts.core.basics.plot import sequential_colormaps

    output = Map()

    # Check
    if connect and density: raise ValueError("Cannot enable 'connect' and 'density' at the same time")

    # Create figure if necessary, get the axes
    if plot is not None: axes = plot.axes
    only_axes = False
    if axes is None:
        from ...core.basics.plot import MPLFigure
        # Create plot
        figure = MPLFigure()
        figure.transparent = transparent
        output.figure = figure
        plot = figure.create_one_plot()
        output.plot = plot
        axes = plot.axes
    else: only_axes = True
    output.axes = axes

    original_xlimits = xlimits
    original_ylimits = ylimits
    output.original_xlimits = original_xlimits
    output.original_ylimits = original_ylimits

    # Add the data
    has_multiple_labels = False
    if types.is_dictionary(x):

        # Checks
        if not types.is_dictionary(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.same_contents(x.keys(), y.keys()): raise ValueError("The curve names must agree in x and y data")

        #colormaps = iter(uniform_colormaps)
        #colormaps = iter(["Greys"] + distinguishable_colormaps)
        colormaps = iter(sequential_colormaps) # go from light to dark
        colors = iter(pretty_colours)

        has_multiple_labels = len(x) > 1

        # Loop over the curves
        for name in x:

            # Clean xy data
            # XLIMITS AND Y LIMITS ARE NOW IN LOG SCALE IF NECESSARY
            _x, _y, xlimits, ylimits = clean_xy_data(x[name], y[name], original_xlimits, original_ylimits, xlog=xlog, ylog=ylog,
                                                   xpositive=xpositive, ypositive=ypositive,
                                                   xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                                                   xnonzero=xnonzero, ynonzero=ynonzero, adjust_limits=False)

            # Plot
            if density:
                cmap = colormaps.next()
                _plot_xy(_x, _y, label=name, connect=connect, plot=plot, density=True, cmap=cmap, size=size)
            else:
                color = colors.next()
                _plot_xy(_x, _y, label=name, connect=connect, plot=plot, color=color, size=size)

    # Sequence
    elif types.is_sequence_or_array(x):

        # Checks
        if not types.is_sequence_or_array(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.equal_sizes(x, y): raise ValueError("The number of x and y points must agree")

        # CLEAN
        x, y, xlimits, ylimits = clean_xy_data(x, y, original_xlimits, original_ylimits, xlog=xlog, ylog=ylog,
                                             xpositive=xpositive, ypositive=ypositive,
                                             xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                                             xnonzero=xnonzero, ynonzero=ynonzero, adjust_limits=False)

        # Plot
        _plot_xy(x, y, connect=connect, density=density, plot=plot, size=size, color=color)

    # Invalid
    else: raise ValueError("Invalid type for x data: '" + str(type(x)) + "'")

    # Add vertical lines
    if vlines is not None:
        for k, vline in enumerate(vlines):
            #if xlog: vline = np.log10(vline)
            if vlinecolor is not None:
                if types.is_sequence(vlinecolor): plot.axvline(x=vline, color=vlinecolor[k], linestyle=vlinestyle)
                else: plot.axvline(x=vline, color=vlinecolor, linestyle=vlinestyle)
            else: plot.axvline(x=vline, linestyle=vlinestyle)

    # Add horizontal lines
    if hlines is not None:
        for l, hline in enumerate(hlines):
            #if ylog: hline = np.log10(hline)
            if hlinecolor is not None:
                if types.is_sequence(hlinecolor): plot.axhline(y=hline, color=hlinecolor[l], linestyle=hlinestyle)
                else: plot.axhline(y=hline, color=hlinecolor, linestyle=hlinestyle)
            else: plot.axhline(y=hline, linestyle=hlinestyle)

    # Set scales
    if xlog: axes.set_xscale("log")
    if ylog: axes.set_yscale("log")

    # Set axes limits
    if xlimits is not None: axes.set_xlim(xlimits[0], xlimits[1])
    if ylimits is not None: axes.set_ylim(ylimits[0], ylimits[1])

    # Set labels
    if x_label is not None:
        if x_color is not None: plot.set_xlabel(x_label, color=x_color)
        else: plot.set_xlabel(x_label)
    if y_label is not None:
        if y_color is not None: plot.set_ylabel(y_label, color=y_color)
        else: plot.set_ylabel(y_label)

    # Set tick colors?
    if x_color is not None: axes.tick_params("x", colors=x_color)
    if y_color is not None: axes.tick_params("y", colors=y_color)

    # Create legend
    if legend and has_multiple_labels: plot.legend()

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Add title
        if title is not None: figure.set_title(title) #plt.title(title)

        # Show or save
        if show is None and path is None: show = True
        if show: figure.show()
        if path is not None:
            output.path = path
            figure.saveto(path)

    # Return the output
    return output

# -----------------------------------------------------------------

def _plot_xy(x, y, label=None, connect=True, density=False, plot=None, cmap=None, color=None, size=None):

    """
    This function ...
    :param x:
    :param y:
    :param label:
    :param connect:
    :param density:
    :param plot:
    :param cmap:
    :param color:
    :param size:
    :return:
    """

    # Set default size
    if size is None: size = 30

    # Connect with lines
    if connect:

        if cmap is not None: raise ValueError("Cannot specify a colormap for line plot")

        if plot is not None: plot.plot(x, y, label=label, color=color)
        else: plt.plot(x, y, label=label, color=color)

    # Scatter points with density
    elif density: _plot_xy_density(x, y, )

    # Scatter points
    else:

        # Cannot specify cmap
        if cmap is not None: raise ValueError("Cannot specify a colormap when density is not plotted")

        # Plot scatter
        if plot is not None: plot.scatter(x, y, label=label, c=color, s=size)
        else: plt.scatter(x, y, label=label, c=color, s=size)

# -----------------------------------------------------------------

def log_scale(x, a=100):
    return np.log10(1+(a*x)) / np.log10(1+a)

# -----------------------------------------------------------------

def simple_log(x):
    return np.log10(x) / np.log10(x.max())

# -----------------------------------------------------------------

def alpha_from_densities(densities, tresfrac=0.05, scaleval=10., logscale=True):

    """
    This function ...
    :param densities:
    :param tresfrac:
    :param scaleval:
    :param logscale:
    :return:
    """

    densities = densities + (np.max(densities) * tresfrac)  # lower bound for alpha

    if logscale:

        if scaleval > 0:

            counts = densities / np.max(densities)
            alpha = log_scale(counts, scaleval)

        else: alpha = simple_log(densities)

    else: alpha = densities / np.max(densities)

    # Return alpha
    return alpha

# -----------------------------------------------------------------

def _plot_xy_density(x, y, axes, method="kde", min_count=1, rug=False, nbins=100, cmap="viridis", xlog=False, ylog=False,
                     xlimits=None, ylimits=None):

    """
    This function ...
    :param x:
    :param y:
    :param axes:
    :param method:
    :param min_count:
    :param rug:
    :param nbins:
    :param cmap:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :return:
    """

    # Contours
    if method == "contours":

        import seaborn as sns

        # Plot
        sns.kdeplot(x, y, ax=axes)

        if rug:
            sns.rugplot(x, color="g", ax=axes)
            sns.rugplot(y, vertical=True, ax=axes)

    # Density field
    elif method == "kde":

        # Warning
        warnings.warn("Calculating density of points: this can take a while ...")

        k = gaussian_kde([x, y])
        xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))

        # Plot density field
        axes.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cmap)

    # Invalid
    else: raise ValueError("Invalid method: '" + method + "'")

# -----------------------------------------------------------------

def plot_xy_density(x, y, title=None, nbins=200, contours=False, path=None, rug=False, transparent=False,
                    xlog=False, ylog=False, xlimits=None, ylimits=None, xpositive=False, ypositive=False,
                    xnonnegative=False, ynonnegative=False, xnonzero=False, ynonzero=False, vlines=None, hlines=None,
                    x_label=None, y_label=None, legend=True, axes=None, plot=None, method=None, show=None):

    """
    This function ...
    :param x:
    :param y:
    :param title:
    :param nbins:
    :param contours:
    :param path:
    :param rug:
    :param transparent:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param vlines:
    :param hlines:
    :param x_label:
    :param y_label:
    :param legend:
    :param axes:
    :param plot:
    :param method:
    :param show:
    :return:
    """

    from ...core.basics.map import Map
    from ...core.basics.plot import sequential_colormaps
    output = Map()

    # Define method
    if method is None:
        if contours or rug: method = "seaborn"
        else: method = "kde" # default method

    # Create figure if necessary, get the axes
    if plot is not None:
        axes = plot.axes
        output.plot = plot
    only_axes = False
    if axes is None:
        from ...core.basics.plot import MPLFigure
        # Create plot
        figure = MPLFigure()
        figure.transparent = transparent
        plot = figure.create_one_plot()
        axes = plot.axes
        output.figure = figure
        output.plot = plot
    else: only_axes = True
    output.axes = axes

    original_xlimits = xlimits
    original_ylimits = ylimits
    output.original_xlimits = original_xlimits
    output.original_ylimits = original_ylimits

    # Multiple data
    if types.is_dictionary(x):

        # Checks
        if not types.is_dictionary(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.same_contents(x.keys(), y.keys()): raise ValueError("The curve names must agree in x and y data")

        colormaps = iter(sequential_colormaps) # go from light to dark
        #colors = iter(pretty_colours)

        # Loop over the curves
        for name in x:

            # Clean xy data
            # XLIMITS AND Y LIMITS ARE NOW IN LOG SCALE IF NECESSARY
            _x, _y, xlimits, ylimits = clean_xy_data(x[name], y[name], original_xlimits, original_ylimits, xlog=xlog,
                                                     ylog=ylog, xpositive=xpositive, ypositive=ypositive,
                                                     xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                                                     xnonzero=xnonzero, ynonzero=ynonzero, adjust_limits=False)

            # Get colormap
            cmap = colormaps.next()

            # Plot density of points
            _plot_xy_density(_x, _y, axes, method=method, rug=rug, nbins=nbins, cmap=cmap, xlog=xlog, ylog=ylog, xlimits=xlimits, ylimits=ylimits)
            #figure.figure.canvas.draw()
            #figure.figure.canvas.draw_idle()

    # Single data set
    elif types.is_sequence_or_array(x):

        # Checks
        if not types.is_sequence_or_array(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.equal_sizes(x, y): raise ValueError("The number of x and y points must agree")

        # Clean xy data
        x, y, xlimits, ylimits = clean_xy_data(x, y, xlimits, ylimits, xlog=xlog, ylog=ylog,
                                               xpositive=xpositive, ypositive=ypositive,
                                               xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                                               xnonzero=xnonzero, ynonzero=ynonzero)

        # Plot density of points
        _plot_xy_density(x, y, axes, method=method, rug=rug, nbins=nbins, xlog=xlog, ylog=ylog, xlimits=xlimits, ylimits=ylimits)

    # Invalid
    else: raise ValueError("Invalid type for x data: '" + str(type(x)) + "'")

    # Add vertical lines
    if vlines is not None:
        for vline in vlines:
            #if xlog: vline = np.log10(vline)
            plot.axvline(x=vline)

    # Add horizontal lines
    if hlines is not None:
        for hline in hlines:
            #if ylog: hline = np.log10(hline)
            plot.axhline(y=hline)

    # Set scales
    if xlog: axes.set_xscale("log")
    if ylog: axes.set_yscale("log")

    # Set axes limits
    if xlimits is not None: axes.set_xlim(xlimits[0], xlimits[1])
    if ylimits is not None: axes.set_ylim(ylimits[0], ylimits[1])

    # Set axes labels
    if x_label is not None: plot.set_xlabel(x_label)
    if y_label is not None: plot.set_ylabel(y_label)

    # Create legend
    if legend: plot.legend()

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Add title
        if title is not None: figure.set_title(title) #plt.title(title)

        # Show or save
        if show is None and path is None: show = True
        if show: figure.show()
        if path is not None: figure.saveto(path)

    # Return the plot?
    #return figure, plot
    return output

# -----------------------------------------------------------------

def vmax_function(array):

    """
    This function ...
    :param array:
    :return:
    """

    # print(1, vmin, vmax)
    from astropy.stats import sigma_clip
    # print(array[array!=0])
    masked_array = sigma_clip(array[array != 0], sigma=4.0)
    values = masked_array.compressed()

    # print(values)
    #vmin = 0.
    # vmin = np.nanmin(values)
    vmax = np.nanmax(values)
    # print(2, vmin, vmax)

    # print(vmin, vmax)
    # vmin = -2
    # vmax = vmax / 5.

    # Return
    return vmax

# -----------------------------------------------------------------

def plot_xy_scatter_density(x, y, title=None, path=None, x_label=None, y_label=None, xlog=False, ylog=False,
                            xlimits=None, ylimits=None, show=None, colormaps=False, colors=None, axes=None, plot=None, dpi=None, color=None,
                            cmap=None, aux=None, aux_name=None, aux_unit=None, aux_log=False, aux_limits=None, density_log=False,
                            legend_location=None, aux_density=True, label=None, log_enhance=50.):

    """
    This function is a scatter density plotting function, using Astrofrog's matplotlib scatter density package
    :param x:
    :param y:
    :param title:
    :param path:
    :param x_label:
    :param y_label:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param show:
    :param colormaps:
    :param colors:
    :param axes:
    :param plot:
    :param dpi:
    :param color:
    :param cmap:
    :param aux:
    :param aux_name:
    :param aux_unit:
    :param aux_log:
    :param aux_limits:
    :param density_log:
    :param legend_location:
    :param aux_density:
    :param label:
    :param log_enhance:
    :return:
    """

    if legend_location is None: legend_location = "upper left"

    # Output
    from ...core.basics.map import Map
    output = Map()

    # Set data
    output.x = x
    output.y = y

    # Set dpi
    if dpi is None: dpi = 72

    from matplotlib.lines import Line2D
    from matplotlib.legend import Legend
    from matplotlib.cm import get_cmap

    from ...core.basics.plot import sequential_colormaps

    # Create figure if necessary, get the axes
    if plot is not None:
        output.plot = plot
        axes = plot.axes
    only_axes = False
    if axes is None:
        from ...core.basics.plot import MPLFigure
        # Create plot
        figure = MPLFigure()
        #figure.transparent = transparent
        plot = figure.create_one_plot(projection="scatter_density")
        axes = plot.axes
        output.figure = figure
        output.plot = plot
    else: only_axes = True
    output.axes = axes
    output.dpi = dpi
    output.xlog = xlog
    output.ylog = ylog

    # Set scales
    if xlog: axes.set_xscale("log")
    if ylog: axes.set_yscale("log")

    original_xlimits = xlimits
    original_ylimits = ylimits

    # Add the data
    if types.is_dictionary(x):

        # Checks
        if not types.is_dictionary(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.same_contents(x.keys(), y.keys()): raise ValueError("The curve names must agree in x and y data")
        if label is not None: raise ValueError("Cannot specify label with multiple datasets")

        if colormaps:
            colors = None
            if colormaps is True: colormaps = iter(sequential_colormaps)  # go from light to dark
            else: pass
            output.colormaps = dict()
        else:
            if colors is None: colors = iter(pretty_colours)
            else:
                if types.is_dictionary(colors): colors = iter(containers.sequence_from_dict(colors, x.keys(), default_iterator=iter(pretty_colours)))
                elif types.is_sequence(colors): colors = iter(colors)
                else: raise ValueError("Invalid type for 'colors': must be dict or sequence with length of the number of datasets")
            colormaps = None
            output.colors = dict()

        for_legend = OrderedDict()

        # Initialize dict
        output.scatters = dict()

        # Loop over the curves
        for name in x:

            # Clean xy data
            # XLIMITS AND Y LIMITS ARE NOW IN LOG SCALE IF NECESSARY
            _x, _y, xlimits, ylimits = clean_xy_data(x[name], y[name], original_xlimits, original_ylimits, xlog=xlog, ylog=ylog, adjust_limits=False)

            # Get colors
            if colormaps is not None:

                # Get next colormap
                cmap_name = colormaps.next()
                cmap = get_cmap(cmap_name)
                base_color = cmap(0)
                average_color = cmap(0.5)
                full_color = cmap(1)

                # Plot scatter points
                scatter = axes.scatter_density(_x, _y, cmap=cmap, dpi=dpi, label=name, vmin=0., vmax=vmax_function, c_alpha_logdensity=density_log) # size=size

                # Set color map name
                output.colormaps[name] = cmap_name

            # Single colors
            elif colors is not None:

                # Get tnext color
                color = colors.next()
                base_color = average_color = full_color = color

                # Plot scatter points
                scatter = axes.scatter_density(_x, _y, color=color, dpi=dpi, label=name, vmin=0., vmax=vmax_function, c_alpha_logdensity=density_log) # size=size

                # Set color name
                output.colors[name] = color

            # Error
            else: raise RuntimeError("Something went wrong")

            # Create fake patch
            patch = Line2D([], [], marker='o', markersize=4, label=name, linewidth=0, markeredgewidth=1, color=average_color)

            # Add legend patch
            for_legend[name] = patch

            # Add scatter
            output.scatters[name] = scatter

        # Create and add legend for the different datasets
        legend = Legend(axes, for_legend.values(), for_legend.keys(), loc=legend_location)
        axes.add_artist(legend)

    # Single data set
    elif types.is_sequence_or_array(x):

        # Checks
        if not types.is_sequence_or_array(y): raise ValueError("The type of x and y data must be equal")
        if not sequences.equal_sizes(x, y): raise ValueError("The number of x and y points must agree")

        # Clean xy data
        #print(x.shape, y.shape)
        #if aux is not None: print(aux.shape)
        x, y, xlimits, ylimits, valid = clean_xy_data(x, y, xlimits, ylimits, xlog=xlog, ylog=ylog, return_valid=True, get_limits=True)

        # Set valid mask
        output.valid = valid

        # Plot scatter points
        if aux is not None:

            # Plot with colormap
            if cmap is None: cmap = "inferno"
            aux = aux[valid]

            # Set interval
            if aux_limits is not None:
                min_aux = aux_limits[0]
                max_aux = aux_limits[1]
            else:
                min_aux = np.nanmin(aux)
                max_aux = np.nanmax(aux)

            # Set norm
            if aux_log:
                from astropy.visualization import LogStretch
                from astropy.visualization.mpl_normalize import ImageNormalize
                norm = ImageNormalize(vmin=min_aux, vmax=max_aux, stretch=LogStretch())
            else: norm = None

            # Set aux limits
            output.aux_limits = (min_aux, max_aux,)

            # Plot
            #print("label", label)
            scatter = axes.scatter_density(x, y, dpi=dpi, vmin=min_aux, vmax=max_aux, c=aux, cmap=cmap, norm=norm, label=label,
                                           c_with_alpha=aux_density, c_alpha_logdensity=density_log, c_alpha_log_enhance=log_enhance)
            scatter_proxy = create_marker_patch(label, "o", get_cmap(cmap)(0), edgewidth=0)

            # Set colormap and other
            output.cmap = cmap
            if aux_name is not None: output.aux_name = aux_name
            if aux_unit is not None: output.aux_unit = aux_unit

            # Add colorbar
            if not only_axes:
                if aux_name is None: aux_name = "Auxilary axis"
                if aux_unit is not None: aux_name = aux_name + " [" + str(aux_unit) + "]"
                figure.figure.colorbar(scatter, label=aux_name)

        # No auxilary axis: plot in color
        else:

            if color is None: color = "blue"
            #print("label", label)
            scatter = axes.scatter_density(x, y, color=color, dpi=dpi, vmin=0., vmax=vmax_function, label=label,
                                           c_alpha_logdensity=density_log)
            scatter_proxy = create_marker_patch(label, "o", color, edgewidth=0)

            # Set color
            output.color = color

        # Set scatter
        output.scatter = scatter
        output.scatter_proxy = scatter_proxy

    # Invalid
    else: raise ValueError("Invalid type for x data: '" + str(type(x)) + "'")

    # Set limits
    if xlimits is not None:
        output.xlimits = xlimits
        axes.set_xlim(xlimits[0], xlimits[1])
    if ylimits is not None:
        output.ylimits = ylimits
        axes.set_ylim(ylimits[0], ylimits[1])

    # Set axes labels
    if x_label is not None:
        output.x_label = x_label
        axes.set_xlabel(x_label)
    if y_label is not None:
        output.y_label = y_label
        axes.set_ylabel(y_label)

    # Axes were not provided: we are supposed to create the whole figure thingy and close it
    if not only_axes:

        # Add title
        if title is not None:
            output.title = title
            figure.set_title(title)  # plt.title(title)

        # Show or save
        if show is None and path is None: show = True
        if show: figure.show()
        if path is not None:
            output.path = path
            figure.saveto(path)

    # Return
    return output

# -----------------------------------------------------------------

def plot_joint(points, title=None, kind="scatter", path=None, xlog=False, ylog=False, xlimits=None, ylimits=None,
               xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False, xnonzero=False, ynonzero=False):
    
    """
    This function ...
    :param points:
    :param title:
    :param kind:
    :param path:
    :param xlog:
    :param ylog:
    :param xlimits:
    :param ylimits:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :return: 
    """

    # Get x, y and labels
    x, y, x_label, y_label = get_xy(points, return_labels=True)

    # Plot
    plot_xy_joint(x, y, kind=kind, title=title, path=path, xlog=xlog, ylog=ylog, xlimits=xlimits, ylimits=ylimits,
                  xpositive=xpositive, ypositive=ypositive, xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                  xnonzero=xnonzero, ynonzero=ynonzero)

# -----------------------------------------------------------------

def plot_xy_joint(x, y, kind="scatter", title=None, path=None, transparent=False, format=None, xlimits=None, ylimits=None,
                  xlog=False, ylog=False, xpositive=False, ypositive=False, xnonnegative=False, ynonnegative=False,
                  xnonzero=False, ynonzero=False):

    """
    This function ...
    :param x:
    :param y:
    :param kind:
    :param title:
    :param path:
    :param transparent:
    :param format:
    :param xlimits:
    :param ylimits:
    :param xlog:
    :param ylog:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :return:
    """

    import seaborn as sns

    # Clean xy dataif not fs.is_directory(dirpath): raise RuntimeError("The directory '" + dirpath + "' does not exist")
    x, y, xlimits, ylimits = clean_xy_data(x, y, xlimits, ylimits, xlog=xlog, ylog=ylog,
                                           xpositive=xpositive, ypositive=ypositive,
                                           xnonnegative=xnonnegative, ynonnegative=ynonnegative,
                                           xnonzero=xnonzero, ynonzero=ynonzero)

    # Create figure
    #plt.figure()

    # Plot
    sns.jointplot(x=x, y=y, kind=kind) #color="k")

    # Add title
    if title is not None: plt.title(title)

    # Show or save
    if path is None: plt.show()
    else: plt.savefig(path, format=format, transparent=transparent)

    # Close
    plt.close()

# -----------------------------------------------------------------

def plot_peak_model(box, x_peak, y_peak, model, title=None, vmin=None, vmax=None):

    """
    This function ...
    :param box:
    :param x_peak:
    :param y_peak:
    :param model:
    :param title:
    :param vmin:
    :param vmax:
    :return:
    """

    # Determine the maximum value in the box and the minimum value for plotting
    if vmin is None: vmin = max(np.nanmin(box), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(box) + vmin)

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    x_peak_pixel = int(round(x_peak))
    y_peak_pixel = int(round(y_peak))

    # Calculate the pixel value at the peak for the data, model and residual
    peak_data_value = box[y_peak_pixel,x_peak_pixel]
    peak_model_value = model(x_peak, y_peak)
    peak_residual_value = peak_data_value - peak_model_value

    # Plot the data with the best-fit model
    plt.figure(figsize=(10,3))
    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax, cmap="viridis")
    plt.plot(x_peak, y_peak, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data " + str(peak_data_value))
    plt.subplot(1,3,2)
    plt.imshow(model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax, cmap="viridis")
    plt.title("Model " + str(peak_model_value))
    plt.subplot(1,3,3)
    plt.imshow(box - model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax, cmap="viridis")
    plt.title("Residual " + str(peak_residual_value))

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    # Show the plot
    plt.show()

# -----------------------------------------------------------------

def plot_star(box, peak, model, title=None, vmin=None, vmax=None):

    """
    This function ...
    :param box:
    :param peak:
    :param model:
    :param title:
    :param vmin:
    :param vmax:
    :return:
    """

    # Normalization
    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    if vmin is None: vmin = max(np.nanmin(box), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(box) + vmin)

    # Evaluate the model and subtract it from the cutout
    evaluated = box.evaluate_model(model)
    subtracted = box - evaluated

    # Create a figure
    plt.figure(figsize=(10,3))

    # Plot the box
    plt.subplot(1,4,1)
    plt.imshow(box, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.plot(peak.x, peak.y, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.xsize-1)
    plt.ylim(0, box.ysize-1)
    plt.title("Cutout")

    # Plot the model
    plt.subplot(1,4,2)
    plt.imshow(evaluated, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
    plt.xlim(0, box.xsize-1)
    plt.ylim(0, box.ysize-1)
    plt.title("Model")

    # Plot the subtracted box on the same scale as the original box and model
    plt.subplot(1,4,3)
    plt.imshow(subtracted, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
    plt.xlim(0, box.xsize-1)
    plt.ylim(0, box.ysize-1)
    plt.title("Residual")

    # Plot the subtracted box on a narrower color scale
    plt.subplot(1,4,4)
    sp = plt.imshow(subtracted, origin='lower', interpolation="nearest", cmap="viridis")
    plt.xlim(0, box.xsize-1)
    plt.ylim(0, box.ysize-1)
    plt.title("Residual")
    plt.colorbar(sp, format="%.2f")

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    # Show the plot
    plt.show()

# -----------------------------------------------------------------

def plot_peaks(box, x_peaks, y_peaks, radius=None, title=None, vmin=None, vmax=None):

    """
    This function plots the data with peaks marked ...
    :param box:
    :param x_peaks:
    :param y_peaks:
    :param radius:
    :param title:
    :param vmin:
    :param vmax:
    :return:
    """

    # Determine the maximum value in the box and the minium value for plotting
    if vmin is None: vmin = max(np.nanmin(box), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(box) + vmin)

    # Set the normalization
    norm = ImageNormalize(stretch=SqrtStretch())

    # Make the plot
    plt.figure(figsize=(8,2.5))
    plt.imshow(box, origin='lower', norm=norm, interpolation='nearest', vmin=vmin, vmax=vmax, cmap="viridis")

    if radius is None: plt.plot(x_peaks, y_peaks, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    else:

        positions = (x_peaks, y_peaks)
        apertures = CircularAperture(positions, r=radius)
        apertures.plot(color='green', lw=1.5, alpha=0.5)

    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)

    if title is not None: plt.title(title)

    plt.show()

# -----------------------------------------------------------------

def plot_peak(box, x_peak, y_peak, radius=None, title=None):

    """
    This function plots the data with peaks marked ...
    :param box:
    :param x_peak:
    :param y_peak:
    :param radius:
    :param title:
    :return:
    """

    plot_peaks(box, [x_peak], [y_peak], radius=radius, title=title)

# -----------------------------------------------------------------

def plot_peaks_models(box, x_peaks, y_peaks, models, vmin=None, vmax=None):

    """
    This function plots the data with peaks marked and models subtracted
    :param box:
    :param x_peaks:
    :param y_peaks:
    :param models:
    :param vmin:
    :param vmax:
    :return:
    """

    # Determine the maximum value in the box and the minium value for plotting
    if vmin is None: vmin = max(np.nanmin(box), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(box) + vmin)

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    # Calculate the sum of all models
    total_model = models[0]
    for i in range(1, len(models)): total_model += models[i]

    # Make the plot
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax, cmap="viridis")
    plt.plot(x_peaks, y_peaks, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")
    plt.subplot(1,3,2)
    plt.imshow(total_model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax, cmap="viridis")
    plt.title("Model")
    plt.subplot(1,3,3)
    plt.imshow(box - total_model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax, cmap="viridis")
    plt.title("Residual")
    plt.show()

# -----------------------------------------------------------------

def plot_star_model(background, background_clipped, est_background, star, est_background_star, fitted_star, vmin=None, vmax=None):

    """
    This function ...
    :param background:
    :param background_clipped:
    :param est_background:
    :param star:
    :param est_background_star:
    :param fitted_star:
    :param vmin:
    :param vmax:
    :return:
    """

    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    if vmin is None: vmin = max(np.nanmin(background), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(background) + vmin)

    # Plot the data with the best-fit model
    plt.figure(figsize=(20,3))
    plt.subplot(1,7,1)
    plt.imshow(background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, background.shape[1]-1)
    plt.ylim(0, background.shape[0]-1)
    plt.title("Background")

    plt.subplot(1,7,2)
    plt.imshow(background_clipped, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, background_clipped.shape[1]-1)
    plt.ylim(0, background_clipped.shape[0]-1)
    plt.title("Sigma-clipped background")

    plt.subplot(1,7,3)
    plt.imshow(est_background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, est_background.shape[1]-1)
    plt.ylim(0, est_background.shape[0]-1)
    plt.title("Estimated background")

    plt.subplot(1,7,4)
    plt.imshow(star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star")

    plt.subplot(1,7,5)
    plt.imshow(star.data - est_background_star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star without background")

    plt.subplot(1,7,6)
    plt.imshow(fitted_star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, fitted_star.shape[1]-1)
    plt.ylim(0, fitted_star.shape[0]-1)
    plt.title("Fitted star")

    plt.subplot(1,7,7)
    plt.imshow(star.data - fitted_star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Residual")

    plt.show()

# -----------------------------------------------------------------

def plot_removal(cutout, mask, background, removed, title=None, vmin=None, vmax=None):

    """
    This function ...
    :param cutout:
    :param mask:
    :param background:
    :param removed:
    :param title:
    :param vmin:
    :param vmax:
    :return:
    """

    # Get raw data of mask as a numpy array
    if hasattr(mask, "data"): maskdata = mask.data
    else: maskdata = mask

    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    if vmin is None: vmin = max(np.nanmin(cutout), 0.)
    if vmax is None: vmax = 0.5 * (np.nanmax(cutout) + vmin)

    # Plot the data with the best-fit model
    plt.figure(figsize=(20,3))
    plt.subplot(1,4,1)
    plt.imshow(cutout, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, cutout.xsize-0.5)
    plt.ylim(-0.5, cutout.ysize-0.5)
    plt.title("Cutout")

    plt.subplot(1,4,2)
    plt.imshow(np.ma.masked_array(cutout, mask=maskdata), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, cutout.xsize-0.5)
    plt.ylim(-0.5, cutout.ysize-0.5)
    plt.title("Background mask")

    plt.subplot(1,4,3)
    plt.imshow(background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, background.xsize-0.5)
    plt.ylim(-0.5, background.ysize-0.5)
    plt.title("Estimated background")

    plt.subplot(1,4,4)
    plt.imshow(removed, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, background.xsize-0.5)
    plt.ylim(-0.5, background.ysize-0.5)
    plt.title("Cutout with star removed")

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    # Show the plot
    plt.show()

# -----------------------------------------------------------------

def plot_source(cutout, mask, background, peaks=None, title=None, show=True, scale="log", frame=None, vmin=None, vmax=None):

    """
    This function ...
    :param cutout:
    :param mask:
    :param background:
    :param peaks:
    :param title:
    :param show:
    :param scale:
    :param frame:
    :param vmin:
    :param vmax:
    :return:
    """

    if scale == "sqrt": norm = ImageNormalize(stretch=SqrtStretch())
    elif scale == "log": norm = ImageNormalize(stretch=LogStretch())
    else: raise ValueError("Invalid scale option")

    if frame is not None:

        gs1 = gridspec.GridSpec(3, 3)
        gs1.update(left=0.05, right=0.48, wspace=0.05)
        ax1 = plt.subplot(gs1[:-1, :])
        ax2 = plt.subplot(gs1[-1, :-1])
        ax3 = plt.subplot(gs1[-1, -1])

    # Determine the maximum value in the box and the minimum value for plotting
    #vmax = np.nanmax(cutout)
    #vmin = np.nanmin(cutout) if vmax <= 0 else 0.0

    vmin = np.nanmin(cutout)
    vmax = 0.5 * (np.nanmax(cutout) + vmin)

    #number = 6 if source_mask is not None else 5

    number = 5

    # Get raw data of mask as a numpy array
    if hasattr(mask, "data"): maskdata = mask.data
    else: maskdata = mask

    # Plot the data with the best-fit model
    plt.figure(figsize=(20,3))
    plt.subplot(1,number,1)
    plt.imshow(cutout, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, cutout.xsize-0.5)
    plt.ylim(-0.5, cutout.ysize-0.5)
    plt.title("Cutout")

    plt.subplot(1,number,2)
    plt.imshow(np.ma.masked_array(cutout, mask=maskdata), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, cutout.xsize-0.5)
    plt.ylim(-0.5, cutout.ysize-0.5)
    plt.title("Masked source")

    plt.subplot(1,number,3)
    plt.imshow(background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, background.xsize-0.5)
    plt.ylim(-0.5, background.ysize-0.5)
    plt.title("Estimated background")

    #plt.subplot(1,number,4)
    #plt.imshow(np.ma.masked_array(cutout, mask=mask.inverse()), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax)
    #plt.xlim(-0.5, cutout.xsize-0.5)
    #plt.ylim(-0.5, cutout.ysize-0.5)
    #plt.title("Masked background")

    plt.subplot(1,number,4)
    plt.imshow(cutout-background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    if peaks is not None: plt.plot(peaks[0], peaks[1], ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(-0.5, cutout.xsize-0.5)
    plt.ylim(-0.5, cutout.ysize-0.5)
    plt.title("Background subtracted")

    #plt.subplot(1,number,6)
    #plt.imshow(np.ma.masked_array(cutout-background, mask=mask.inverse()), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax)
    #plt.xlim(-0.5, cutout.xsize-0.5)
    #plt.ylim(-0.5, cutout.ysize-0.5)
    #plt.title("Background subtracted source")

    replaced = cutout.copy()
    replaced[maskdata] = background[maskdata]

    plt.subplot(1,number,5)
    plt.imshow(replaced, origin="lower", interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(-0.5, cutout.xsize - 0.5)
    plt.ylim(-0.5, cutout.ysize - 0.5)
    plt.title("Removed source (replaced by background)")

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    # Show the plot
    if show: plt.show()

# -----------------------------------------------------------------

def plot_background_subtraction(background, background_clipped, est_background, star, est_background_star):

    """
    This function ...
    :param background:
    :param background_clipped:
    :param est_background:
    :param star:
    :param est_background_star:
    :return:
    """

    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.nanmax(background)
    vmin = np.nanmin(background) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(20,3))
    plt.subplot(1,5,1)
    plt.imshow(background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, background.shape[1]-1)
    plt.ylim(0, background.shape[0]-1)
    plt.title("Background")

    plt.subplot(1,5,2)
    plt.imshow(background_clipped, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, background_clipped.shape[1]-1)
    plt.ylim(0, background_clipped.shape[0]-1)
    plt.title("Sigma-clipped background")

    plt.subplot(1,5,3)
    plt.imshow(est_background, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, est_background.shape[1]-1)
    plt.ylim(0, est_background.shape[0]-1)
    plt.title("Estimated background")

    plt.subplot(1,5,4)
    plt.imshow(star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star")

    plt.subplot(1,5,5)
    plt.imshow(star.data - est_background_star, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star without background")

    plt.show()

# -----------------------------------------------------------------

def plot_background_center(cutout, mask, peaks=None, title=None, show=True, scale="sqrt"):

    """
    This function ...
    :param cutout:
    :param mask:
    :param peaks:
    :param title:
    :param show:
    :param scale:
    :return:
    """

    if scale == "sqrt": norm = ImageNormalize(stretch=SqrtStretch())
    elif scale == "log": norm = ImageNormalize(stretch=LogStretch())
    else: raise ValueError("Invalid scale option")

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.nanmax(cutout)
    vmin = np.nanmin(cutout) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(10,4))
    plt.subplot(1,3,1)
    plt.imshow(cutout, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0.5, cutout.xsize-0.5)
    plt.ylim(0.5, cutout.ysize-0.5)
    plt.title("Cutout")

    # Get raw data of mask as a numpy array
    if hasattr(mask, "data"): maskdata = mask.data
    else: maskdata = mask

    plt.subplot(1,3,2)
    plt.imshow(np.ma.masked_array(cutout, mask=maskdata), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0.5, cutout.xsize-0.5)
    plt.ylim(0.5, cutout.ysize-0.5)
    plt.title("Masked source")

    plt.subplot(1,3,3)
    plt.imshow(np.ma.masked_array(cutout, mask=np.logical_not(maskdata)), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    if peaks is not None: plt.plot(peaks[0], peaks[1], ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0.5, cutout.xsize-0.5)
    plt.ylim(0.5, cutout.ysize-0.5)
    plt.title("Masked background")

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    # Show the plot
    if show: plt.show()

# -----------------------------------------------------------------

def plot_difference(box_a, box_b, share_colorscale=False, title=None):

    """
    This function ...
    :param box_a:
    :param box_b:
    :param share_colorscale:
    :param title:
    :return:
    """

    #norm = ImageNormalize(stretch=SqrtStretch())
    norm = ImageNormalize(stretch=LogStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    #vmax = np.nanmax(box_a)
    #vmin = np.nanmin(box_a) if vmax <= 0 else 0.0

    vmin = np.nanmin(box_a)
    vmax = 0.5 * (np.nanmax(box_b) + vmin)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    #plt.imshow(box_a, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.imshow(box_a, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, box_a.shape[1]-1)
    plt.ylim(0, box_a.shape[0]-1)
    plt.title("Data a")
    plt.subplot(1,3,2)
    #plt.imshow(box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.imshow(box_b, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
    plt.xlim(0, box_a.shape[1]-1)
    plt.ylim(0, box_a.shape[0]-1)
    plt.title("Data b")
    plt.subplot(1,3,3)

    if share_colorscale:

        plt.imshow(box_a - box_b, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
        plt.xlim(0, box_a.shape[1]-1)
        plt.ylim(0, box_a.shape[0]-1)
        plt.title("Residual")
        #plt.imshow(box_a - box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    else:

        residualimage = plt.imshow(box_a - box_b, origin='lower', interpolation="nearest", cmap="viridis")
        plt.xlim(0, box_a.shape[1]-1)
        plt.ylim(0, box_a.shape[0]-1)
        plt.title("Residual")
        plt.colorbar(residualimage, format="%.2f")

    # Set the main title
    if title is not None: plt.suptitle(title, size=16)

    plt.show()

# -----------------------------------------------------------------

def plot_difference_value(box, value, share_colorscale=False):

    """
    This function ...
    :param box:
    :param value:
    :param share_colorscale:
    :return:
    """

    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.nanmax(box)
    vmin = np.nanmin(box) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    #plt.imshow(box_a, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.imshow(box, origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")
    plt.subplot(1,3,2)
    #plt.imshow(box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    value_box = np.full(box.shape, value)

    #print np.median(box-value_box)

    plt.imshow(value_box, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
    plt.title("Constant value")
    plt.subplot(1,3,3)

    if share_colorscale:

        plt.imshow(box - value_box, origin='lower', interpolation="nearest", norm=norm, vmin=0.0, vmax=vmax, cmap="viridis")
        plt.title("Residual")
        #plt.imshow(box_a - box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    else:

        residualimage = plt.imshow(box - value_box, origin='lower', interpolation="nearest", cmap="viridis")
        plt.title("Residual")
        plt.colorbar(residualimage, format="%.2f")

    plt.show()

# -----------------------------------------------------------------

def plot_difference_model(box, model):

    """
    This function ...
    :param box:
    :param model:
    :return:
    """

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.nanmax(box)
    vmin = np.nanmin(box) if vmax <= 0 else 0.0

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    # Evaluate the model in the box
    model_box = model(x_plotvalues, y_plotvalues)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))

    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation="nearest", vmin=vmin, vmax=vmax, cmap="viridis")
    plt.title("Data")

    plt.subplot(1,3,2)
    plt.imshow(model_box, origin='lower', interpolation="nearest", vmin=vmin, vmax=vmax, cmap="viridis")
    plt.title("Model")

    plt.subplot(1,3,3)
    plt.imshow(box - model_box, origin='lower', interpolation="nearest", vmin=vmin, vmax=vmax, cmap="viridis")
    plt.title("Residual")

    plt.show()

# -----------------------------------------------------------------

def plot_table(filepath, column_x, column_y, output_path, x_log=False, y_log=False):

    """
    This function ...
    :param filepath:
    :param column_x:
    :param column_y:
    :param x_log:
    :param y_log:
    :param output_path:
    :return:
    """

    from ...core.tools import terminal

    # Set log flags
    xlog = "true" if x_log else "false"
    ylog = "true" if y_log else "false"

    # Make command
    command = "topcat -stilts plot2plane xlog=" + xlog + " ylog=" + ylog + " in=" + filepath + " ifmt=ASCII x=" + \
              column_x + " y=" + column_y + " legend=false layer_1=Mark layer_2=LinearFit out=" + output_path
    #print(command)

    # Execute the plotting command
    terminal.execute(command)

# -----------------------------------------------------------------

def plot_stilts(filepaths, xcolumn, ycolumn, xlabel, ylabel, path=None, title=None, xlimits=None, ylimits=None, legend=True,
                xlog=False, ylog=False):

    """
    This function ...
    :param filepaths:
    :param xcolumn:
    :param ycolumn:
    :param xlabel:
    :param ylabel:
    :param path:
    :param title:
    :param xlimits:
    :param ylimits:
    :param legend:
    :param xlog:
    :param ylog:
    :return:
    """

    from ...core.tools import terminal

    # xpix=932 ypix=371 \
    # xcrowd=0.9998301109057076 ycrowd=0.9998301109057076

    # Set log flags
    xlog = "true" if xlog else "false"
    ylog = "true" if ylog else "false"

    # Construct command
    command = "topcat -stilts plot2plane"

    # Set scales
    command += " xlog=" + xlog + " ylog=" + ylog

    # Add labels
    command += " xlabel='" + xlabel + "' ylabel='" + ylabel + "'"

    if xlimits is not None: command += " xmin=" + str(xlimits[0]) + " xmax=" + str(xlimits[1])
    if ylimits is not None: command += " ymin=" + str(ylimits[0]) + " ymax=" + str(ylimits[1])

    if legend: command += " legend=true"

    command += " ifmt=ASCII x=" + xcolumn + " y=" + ycolumn + " shading=density"

    # Loop over the table files
    colors = iter(["white", "green", "grey"])
    densemaps = iter(["greyscale", "viridis", "plasma"])
    for index, name in enumerate(filepaths.keys()):

        filepath = filepaths[name]

        color = colors.next()
        densemap = densemaps.next()

        line = " layer_" + str(index+1) + "=Mark in_" + str(index+1) + "='" + filepath + "' color_" + str(index+1) + "=" + color + " densemap_" + str(index+1) + "=" + densemap + " leglabel_" + str(index+1) + "='" + name + "'"
        command += line

    # Set output path
    if path is not None: command += " out='" + path + "'"

    # Debugging
    log.debug("Plotting command:")
    if log.is_debug: print(command)

    # Execute the plotting command
    terminal.execute(command)

# -----------------------------------------------------------------

def clean_xy_data(x, y, xlimits=None, ylimits=None, xlog=False, ylog=False, xpositive=False, ypositive=False,
                  xnonnegative=False, ynonnegative=False, xnonzero=False, ynonzero=False, adjust_limits=False,
                  apply_log=False, return_valid=False, get_limits=False):

    """
    This function ...
    :param x:
    :param y:
    :param xlimits:
    :param ylimits:
    :param xlog:
    :param ylog:
    :param xpositive:
    :param ypositive:
    :param xnonnegative:
    :param ynonnegative:
    :param xnonzero:
    :param ynonzero:
    :param adjust_limits:
    :param apply_log:
    :param return_valid:
    :param get_limits:
    :return:
    """

    # IF LOG, ONLY POSITIVE VALUES ARE ALLOWED!
    if xlog: xpositive = True
    if ylog: ypositive = True

    # Create mask of valid data points
    valid_x = np.isfinite(x)
    valid_y = np.isfinite(y)

    # Only positive values?
    if xpositive: valid_x *= (x > 0)
    if ypositive: valid_y *= (y > 0)

    # Only non-negative values?
    if xnonnegative: valid_x *= (x >= 0)
    if ynonnegative: valid_y *= (y >= 0)

    # No zero values
    if xnonzero: valid_x *= (x != 0)
    if ynonzero: valid_y *= (y != 0)

    # OUTSIDE LIMITS? -> ALREADY CLIP OUT TO POTENTIALLY SAVE RENDERING/CALCULATION TIME & MEMORY
    if xlimits is not None: valid_x *= (x >= xlimits[0]) * (x <= xlimits[1])
    if ylimits is not None: valid_y *= (y >= ylimits[0]) * (y <= ylimits[1])

    # Create combined mask
    valid = valid_x * valid_y

    # Keep only the valid data
    x = x[valid]
    y = y[valid]

    # Check length
    if len(x) == 0: raise ValueError("None of the data passes the conditions")

    # Make into log
    if apply_log:

        # X
        if xlog:
            x = np.log10(x)
            if xlimits is not None: xlimits = (np.log10(xlimits[0]), np.log10(xlimits[1]),)

        # Y
        if ylog:
            y = np.log10(y)
            if ylimits is not None: ylimits = (np.log10(ylimits[0]), np.log10(ylimits[1]),)

    # ADJUST LIMITS? -> SHRINK TO THE DATA
    if adjust_limits:
        xlimits = (np.min(x), np.max(x),)
        ylimits = (np.min(y), np.max(y),)
    elif get_limits:
        if xlimits is None: xlimits = (np.min(x), np.max(x),)
        if ylimits is None: ylimits = (np.min(y), np.max(y),)

    # Return cleaned data
    if return_valid: return x, y, xlimits, ylimits, valid
    else: return x, y, xlimits, ylimits

# -----------------------------------------------------------------

def align_marker(marker, halign='center', valign='middle',):
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """

    if isinstance(halign, (str, unicode)):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, (str, unicode)):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)

# -----------------------------------------------------------------

def create_marker_patch(label, marker, color, edgecolor="black", edgewidth=1, size=7):

    """
    This function ...
    :param marker:
    :param color
    :param edgecolor:
    :param edgewidth:
    :param size:
    :param label:
    :return:
    """

    from matplotlib.lines import Line2D

    # Create and return
    return Line2D([], [], marker=marker, markersize=size, label=label, linewidth=0,
                  markeredgecolor=edgecolor, markerfacecolor=color, markeredgewidth=edgewidth)

# -----------------------------------------------------------------

class Condition(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, lower=None, upper=None, colname=None):

        """
        The constructor ...
        :param lower:
        :param upper:
        :param colname:
        """

        # Set
        self.lower = lower
        self.upper = upper

        # Set colname
        self.colname = colname

# -----------------------------------------------------------------

class FilterCondition(Condition):

    """
    This class ...
    """

# -----------------------------------------------------------------

class RemoveCondition(Condition):

    """
    This class ...
    """

    def __init__(self, lower=None, upper=None, zeroes=False, nans=False, colname=None):
        """
        This function ...
        :param lower:
        :param upper:
        :param zeroes:
        :param nans:
        :param colname:
        """

        # Call the constructor of the base class
        super(RemoveCondition, self).__init__(lower=lower, upper=upper, colname=colname)

        # Set flags
        self.zeroes = zeroes
        self.nans = nans

# -----------------------------------------------------------------

class MultiRemoveCondition(object):

    """
    This class ...
    """

    def __init__(self, *conditions):
        """
        This function ...
        :param conditions:
        """

        # Set the conditions
        self.conditions = conditions

# -----------------------------------------------------------------

def get_mask_for_conditions(scatter, conditions):

    """
    This function ...
    :param scatter:
    :param conditions:
    :return:
    """

    # Initialize mask
    mask = np.ones_like(scatter.x_array, dtype=bool)

    # Loop over the columns in the conditions
    for key in conditions:

        # Get condition and actual colname (in case dictionary keys are just labels for the conditions)
        condition = conditions[key]

        # Filter condition (positive condition)
        if isinstance(condition, FilterCondition):

            # Get colname
            if condition.colname is not None: colname = condition.colname
            else: colname = key

            # Get column array (potentially with div/mul/add/sub symbol)
            array, array_unit = get_aux(scatter, colname)

            if condition.lower is not None:

                lower = condition.lower
                if hasattr(lower, "unit"):
                    if array_unit is None: raise ValueError("Condition for '" + colname + "' has unit but column has no unit")
                    mask *= array > lower.to(array_unit).value
                else: mask *= array > lower

            if condition.upper is not None:

                upper = condition.upper
                if hasattr(upper, "unit"):
                    if array_unit is None: raise ValueError("Condition for '" + colname + "' has unit but column has no unit")
                    mask *= array < upper.to(array_unit).value
                else: mask *= array < upper

        # Remove condition (negative condition)
        elif isinstance(condition, RemoveCondition):

            # Get colname
            if condition.colname is not None: colname = condition.colname
            else: colname = key

            # Get column array
            array, array_unit = get_aux(scatter, colname)

            # Initialize remove mask
            remove = np.ones_like(scatter.x_array, dtype=bool)

            if condition.lower is not None:

                lower = condition.lower
                if hasattr(lower, "unit"):
                    if array_unit is None: raise ValueError("Condition for '" + colname + "' has unit but column has no unit")
                    remove *= array > lower.to(array_unit).value
                else: remove *= array > lower

            if condition.upper is not None:

                upper = condition.upper
                if hasattr(upper, "unit"):
                    if array_unit is None: raise ValueError("Condition for '" + colname + "' has unit but column has no unit")
                    remove *= array < upper.to(array_unit).value
                else: remove *= array < upper

            # Add condition
            mask *= np.logical_not(remove)

        # Multi-remove condition
        elif isinstance(condition, MultiRemoveCondition):

            # Initialize remove mask
            remove = np.ones_like(scatter.x_array, dtype=bool)

            # Loop over the sub-conditions
            for cond in condition.conditions:

                # Get column name
                coln = cond.colname
                if coln is None: raise ValueError("Column name not defined for some sub-conditions")

                # Get column array
                array, array_unit = get_aux(scatter, coln)

                if cond.lower is not None:

                    lower = cond.lower
                    if hasattr(lower, "unit"):
                        if array_unit is None: raise ValueError("Condition for '" + coln + "' has unit but column has no unit")
                        remove *= array > lower.to(array_unit).value
                    else: remove *= array > lower

                if cond.upper is not None:

                    upper = cond.upper
                    if hasattr(upper, "unit"):
                        if array_unit is None: raise ValueError("Condition for '" + coln + "' has unit but column has no unit")
                        remove *= array < upper.to(array_unit).value
                    else: remove *= array < upper

            # Add condition
            mask *= np.logical_not(remove)

        # Invalid
        else: raise ValueError("Invalid type for condition")

    # Return mask
    return mask

# -----------------------------------------------------------------
