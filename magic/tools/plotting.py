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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Import astronomical modules
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import MinMaxInterval, ZScaleInterval
from photutils import CircularAperture

# Import the relevant PTS classes and modules


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

def color_dict(gradient):

  """
  Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on
    """

  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
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
    
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)

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

def plot_mask(mask, title=None, path=None, format=None):

    """
    This function ...
    :param mask:
    :param title:
    :param path:
    :param format:
    :return:
    """

    # Make the plot
    plt.figure(figsize=(7,7))
    plt.imshow(mask, origin="lower", interpolation="nearest", cmap='Greys')
    plt.xlim(0, mask.shape[1] - 1)
    plt.ylim(0, mask.shape[0] - 1)

    if title is not None: plt.title(title)
    else: plt.title("Black means True")

    if path is None: plt.show()
    else: plt.savefig(path, format=format)

    plt.close()

# -----------------------------------------------------------------

def plot_box(box, title=None, path=None, format=None, scale="log", interval="pts", cmap="viridis"):

    """
    This function ...
    :param box:
    :param title:
    :param path:
    :param format:
    :param scale:
    :param interval:
    :param cmap:
    :return:
    """

    # Other new colormaps: plasma, magma, inferno

    # Normalization
    if scale == "log": norm = ImageNormalize(stretch=LogStretch())
    elif scale == "sqrt": norm = ImageNormalize(stretch=SqrtStretch())
    #elif scale == "skimage": norm = exposure.equalize_hist
    else: raise ValueError("Invalid option for 'scale'")

    if interval == "zscale":

        vmin, vmax = ZScaleInterval().get_limits(box)

    elif interval == "pts":

        # Determine the maximum value in the box and the mimimum value for plotting
        vmin = max(np.nanmin(box), 0.)
        vmax = 0.5 * (np.nanmax(box) + vmin)

    elif interval == "minmax":

        vmin, vmax = MinMaxInterval().get_limits(box)

    elif isinstance(interval, tuple):

        vmin = interval[0]
        vmax = interval[1]

    else: raise ValueError("Invalid option for 'interval'")

    # Make the plot
    plt.figure(figsize=(7,7))
    plt.imshow(box, origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)

    if title is not None: plt.title(title)

    if path is None: plt.show()
    else: plt.savefig(path, format=format)

    plt.close()

    # Return vmin and vmax
    return vmin, vmax

# -----------------------------------------------------------------

def plot_peak_model(box, x_peak, y_peak, model, title=None, vmin=None, vmax=None):

    """
    This function ...
    :param box:
    :param x_peak:
    :param y_peak:
    :param model:
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
    :param x_peaks:
    :param y_peaks:
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
    :param interpolated_background:
    :param star:
    :param interpolated_background_star:
    :param fitted_star:
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
    :return:
    """

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
    plt.imshow(np.ma.masked_array(cutout, mask=mask), origin='lower', interpolation="nearest", norm=norm, vmin=vmin, vmax=vmax, cmap="viridis")
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
    :param background:
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
    :param peaks:
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
