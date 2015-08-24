#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import astronomical modules
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# *****************************************************************

def plot_box(box):

    """
    This function ...
    :param box:
    :return:
    """

    # Determine the maximum value in the box and the mimimum value for plotting
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Make the plot
    plt.figure(figsize=(8,2.5))
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")

    plt.show()

# *****************************************************************

def plot_peak_model(box, x_peak, y_peak, model):

    """
    This function ...
    :param box:
    :param x_peak:
    :param y_peak:
    :param model:
    :return:
    """

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    # Calcualte the pixel value at the peak for the data, model and residual
    x_peak_pixel = int(round(x_peak))
    y_peak_pixel = int(round(y_peak))
    peak_data_value = box[y_peak_pixel,x_peak_pixel]
    peak_model_value = model(x_plotvalues, y_plotvalues)[y_peak_pixel,x_peak_pixel]
    peak_residual_value = peak_data_value - peak_model_value

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.plot(x_peak, y_peak, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data " + str(peak_data_value))
    plt.subplot(1,3,2)
    plt.imshow(model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.title("Model " + str(peak_model_value))
    plt.subplot(1,3,3)
    plt.imshow(box - model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.title("Residual " + str(peak_residual_value))
    plt.show()

    # Grey-scale plot:
    #from astropy.visualization import SqrtStretch
    #from astropy.visualization.mpl_normalize import ImageNormalize
    #norm = ImageNormalize(stretch=SqrtStretch())
    #plt.imshow(square, cmap='Greys_r', origin='lower', norm=norm)

# *****************************************************************

def plot_peaks(box, x_peaks, y_peaks):

    """
    This function plots the data with peaks marked ...
    :param box:
    :param x_peaks:
    :param y_peaks:
    :return:
    """

    # Determine the maximum value in the box and the minium value for plotting
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Make the plot
    plt.figure(figsize=(8,2.5))
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.plot(x_peaks, y_peaks, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")
    plt.show()

# *****************************************************************

def plot_peaks_models(box, x_peaks, y_peaks, models):

    """
    This function plots the data with peaks marked and models subtracted
    :param box:
    :param x_peaks:
    :param y_peaks:
    :param models:
    :return:
    """

    # Determine the maximum value in the box and the minium value for plotting
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    # Calculate the sum of all models
    total_model = models[0]
    for i in range(1, len(models)): total_model += models[i]

    # Make the plot
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.plot(x_peaks, y_peaks, ls='none', color='white', marker='+', ms=40, lw=10, mew=4)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")
    plt.subplot(1,3,2)
    plt.imshow(total_model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.title("Model")
    plt.subplot(1,3,3)
    plt.imshow(box - total_model(x_plotvalues, y_plotvalues), origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.title("Residual")
    plt.show()

# *****************************************************************

def plot_star_model(background, background_clipped, est_background, star, est_background_star, fitted_star):

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
    vmax = np.max(background)
    vmin = np.min(background) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(20,3))
    plt.subplot(1,7,1)
    plt.imshow(background, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, background.shape[1]-1)
    plt.ylim(0, background.shape[0]-1)
    plt.title("Background")

    plt.subplot(1,7,2)
    plt.imshow(background_clipped, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, background_clipped.shape[1]-1)
    plt.ylim(0, background_clipped.shape[0]-1)
    plt.title("Sigma-clipped background")

    plt.subplot(1,7,3)
    plt.imshow(est_background, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, est_background.shape[1]-1)
    plt.ylim(0, est_background.shape[0]-1)
    plt.title("Estimated background")

    plt.subplot(1,7,4)
    plt.imshow(star, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star")

    plt.subplot(1,7,5)
    plt.imshow(star.data - est_background_star, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Star without background")

    plt.subplot(1,7,6)
    plt.imshow(fitted_star, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, fitted_star.shape[1]-1)
    plt.ylim(0, fitted_star.shape[0]-1)
    plt.title("Fitted star")

    plt.subplot(1,7,7)
    plt.imshow(star.data - fitted_star, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, star.shape[1]-1)
    plt.ylim(0, star.shape[0]-1)
    plt.title("Residual")

    plt.show()

# *****************************************************************

def plot_difference(box_a, box_b, share_colorscale=False):

    """
    This function ...
    :param box_a:
    :param box_b:
    :param share_colorscale:
    :return:
    """

    norm = ImageNormalize(stretch=SqrtStretch())

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.max(box_a)
    vmin = np.min(box_a) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    #plt.imshow(box_a, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.imshow(box_a, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, box_a.shape[1]-1)
    plt.ylim(0, box_a.shape[0]-1)
    plt.title("Data a")
    plt.subplot(1,3,2)
    #plt.imshow(box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)
    plt.imshow(box_b, origin='lower', interpolation='none', norm=norm, vmin=0.0, vmax=vmax)
    plt.xlim(0, box_a.shape[1]-1)
    plt.ylim(0, box_a.shape[0]-1)
    plt.title("Data b")
    plt.subplot(1,3,3)

    if share_colorscale:

        plt.imshow(box_a - box_b, origin='lower', interpolation='none', norm=norm, vmin=0.0, vmax=vmax)
        plt.xlim(0, box_a.shape[1]-1)
        plt.ylim(0, box_a.shape[0]-1)
        plt.title("Residual")
        #plt.imshow(box_a - box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    else:

        residualimage = plt.imshow(box_a - box_b, origin='lower', interpolation='none')
        plt.xlim(0, box_a.shape[1]-1)
        plt.ylim(0, box_a.shape[0]-1)
        plt.title("Residual")
        plt.colorbar(residualimage, format="%.2f")

    plt.show()

# *****************************************************************

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
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))
    plt.subplot(1,3,1)
    #plt.imshow(box_a, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.imshow(box, origin='lower', interpolation='none', norm=norm, vmin=vmin, vmax=vmax)
    plt.xlim(0, box.shape[1]-1)
    plt.ylim(0, box.shape[0]-1)
    plt.title("Data")
    plt.subplot(1,3,2)
    #plt.imshow(box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    value_box = np.full(box.shape, value)

    #print np.median(box-value_box)

    plt.imshow(value_box, origin='lower', interpolation='none', norm=norm, vmin=0.0, vmax=vmax)
    plt.title("Constant value")
    plt.subplot(1,3,3)

    if share_colorscale:

        plt.imshow(box - value_box, origin='lower', interpolation='none', norm=norm, vmin=0.0, vmax=vmax)
        plt.title("Residual")
        #plt.imshow(box_a - box_b, origin='lower', interpolation='nearest', vmin=0.0, vmax=vmax)

    else:

        residualimage = plt.imshow(box - value_box, origin='lower', interpolation='none')
        plt.title("Residual")
        plt.colorbar(residualimage, format="%.2f")

    plt.show()

# *****************************************************************

def plot_difference_model(box, model):

    """
    This function ...
    :param box:
    :param model:
    :return:
    """

    # Determine the maximum value in the box and the minimum value for plotting
    vmax = np.max(box)
    vmin = np.min(box) if vmax <= 0 else 0.0

    # Create x and y meshgrid for plotting
    y_plotvalues, x_plotvalues = np.mgrid[:box.shape[0], :box.shape[1]]

    # Evaluate the model in the box
    model_box = model(x_plotvalues, y_plotvalues)

    # Plot the data with the best-fit model
    plt.figure(figsize=(8,2.5))

    plt.subplot(1,3,1)
    plt.imshow(box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.title("Data")

    plt.subplot(1,3,2)
    plt.imshow(model_box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.title("Model")

    plt.subplot(1,3,3)
    plt.imshow(box - model_box, origin='lower', interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.title("Residual")

# *****************************************************************
