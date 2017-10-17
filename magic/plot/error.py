#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.plot.error Contains the ErrorPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
from scipy import ndimage
import copy
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.gridspec as gridspec
import glob
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import LogNorm
from collections import OrderedDict
from textwrap import wrap

from astropy.io import fits
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

#import wcsaxes
import matplotlib.colors as mpl_colors
import matplotlib.colorbar as mpl_colorbar

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...magic.tools.plotting import colours

# -----------------------------------------------------------------

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["dodgerblue", "r", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen",
                 "lightcoral", "crimson", "saddlebrown", "mediumslateblue", "dive", "lightslategrey", "firebrick",
                 "orange", "darkcyan", "hotpink", "indianred", "aqua"]
rainbow = cm.get_cmap("rainbow")

# -----------------------------------------------------------------

class ErrorPlotter(object):

    """
    This class ...
    """

    def __init__(self, title=None):

        """
        The constructor ...
        :param title:
        """

        # Set the title
        self.title = title

        # Figure and grid
        self._figure = None
        self._grid = None

        # Properties
        self.style = "dark" # "dark" or "light"
        self.transparent = True
        self.format = None
        self.colormap = "viridis"
        self.vmin = None

    # -----------------------------------------------------------------

    def add_image(self, image):

        """
        This function ...
        :return:
        """



    # -----------------------------------------------------------------

    def run(self, path):

        """
        This function ...
        :param path:
        :return:
        """

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        figure = plt.figure()

        # figure, (ax1, ax2, ax3) = plt.subplots(3)

        # fluxes = dict()
        # total_errors = dict()
        # poisson_errors = dict()
        # calibration_errors = dict()
        # sky_errors = dict()
        # colors = dict()

        # Loop over the various bands and get the flux values and error contributions for each pixel
        # counter = 0
        # for label in self.sorted_labels:
        # Determine the color
        #    colors[label] = pretty_colors[counter]
        #   not_nan = Mask.is_nan(self.images[label]).inverse()
        # fluxes[label] = self.images[label][not_nan].flatten()
        #     counter += 1

        ax1 = figure.add_subplot(4, 1, 1)  # Poisson
        ax2 = figure.add_subplot(4, 1, 2)  # Calibration
        ax3 = figure.add_subplot(4, 1, 3)  # Sky
        ax4 = figure.add_subplot(4, 1, 4)  # Total

        # ax1.set_xscale('log')
        ax2.set_yscale('log')

        counter = 0
        number_of_images = len(self.sorted_labels)
        for label in self.sorted_labels:

            color = rainbow(float(counter) / float(number_of_images))

            flux_notnan = Mask.is_nan(self.images[label]).inverse()
            error_notnan = Mask.is_nan(self.errors[label]).inverse()
            notnan = flux_notnan + error_notnan

            fluxes = self.images[label][notnan]
            errors = self.errors[label][notnan]
            rel_errors = errors / fluxes
            notinf = np.logical_not(np.isinf(rel_errors))

            fluxes = fluxes[notinf]
            rel_errors = rel_errors[notinf]

            if label in self.poisson_errors:
                rel_poisson_errors = self.poisson_errors[label][notnan][notinf] / fluxes
                notinf_poisson = np.logical_not(np.isinf(rel_poisson_errors))
                ax1.scatter(fluxes[notinf_poisson], rel_poisson_errors[notinf_poisson], color=color)

            if label in self.calibration_errors:
                rel_calibration_errors = self.calibration_errors[label][notnan][notinf] / fluxes

                notinf_calibration = np.logical_not(np.isinf(rel_calibration_errors))

                ax2.scatter(fluxes[notinf_calibration], rel_calibration_errors[notinf_calibration], color=color)

            if label in self.sky_errors:
                rel_sky_errors = self.sky_errors[label][notnan][notinf] / fluxes
                notinf_sky = np.logical_not(np.isinf(rel_sky_errors))
                ax3.scatter(fluxes[notinf_sky], rel_sky_errors[notinf_sky], color=color)

            ax4.scatter(fluxes, rel_errors)

            if counter == 5: break

            counter += 1

            # Loop over the various bands
            # for label in fluxes:

            #    ax1.scatter(fluxes[label], poisson_errors[label], color=colors[label], label=label)
            #    ax2.scatter(fluxes[label], calibration_errors[label], color=colors[label], label=label)
            #    ax3.scatter(fluxes[label], sky_errors[label], color=colors[label], label=label)
            #    ax4.scatter(fluxes[label], total_errors[label], color=colors[label], label=label)

        fig_a.hist2d(np.log10(old + yng + new), np.log10(tot), cmin=1, norm=LogNorm(), bins=200, cmap=plt.cm.BuGn_r)

# -----------------------------------------------------------------
