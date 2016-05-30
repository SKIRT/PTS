#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.scatter Contains the ScatterPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict
from textwrap import wrap
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["dodgerblue", "r", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class ScatterPlotter(object):
    
    """
    This class ...
    """

    def __init__(self, title=None, data=None):

        """
        This function ...
        :return:
        """

        # Set the title
        self.title = title

        if data is None:
            self.x = []
            self.y = []
            self.z = []
        else:
            self.x = data[0]
            self.y = data[1]
            self.z = data[2]

        # The axes labels
        self.x_label = None
        self.y_label = None
        self.z_label = None

        # The axes limits
        self.x_limits = [None, None]
        self.y_limits = [None, None]
        self.z_limits = [None, None]

        # Store the figure and its axes as references
        self._figure = None

        # Properties
        self.color_map = "viridis"
        self.format = None
        self.transparent = False
        self.density = True

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def add_point(self, x, y, z):

        """
        This function ...
        :param x:
        :param y:
        :param z:
        :return:
        """

        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

    # -----------------------------------------------------------------

    def set_x_limits(self, x_min, x_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :return:
        """

        self.x_limits[0] = x_min
        self.x_limits[1] = x_max

    # -----------------------------------------------------------------

    def set_y_limits(self, y_min, y_max):

        """
        This function ...
        :param y_min:
        :param y_max:
        :return:
        """

        self.y_limits[0] = y_min
        self.y_limits[1] = y_max

    # -----------------------------------------------------------------

    def set_z_limits(self, z_min, z_max):

        """
        This function ...
        :param z_min:
        :param z_max:
        :return:
        """

        self.z_limits[0] = z_min
        self.z_limits[1] = z_max

    # -----------------------------------------------------------------

    def set_x_label(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        self.x_label = label

    # -----------------------------------------------------------------

    def set_y_label(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        self.y_label = label

    # -----------------------------------------------------------------

    def set_z_label(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        self.z_label = label

    # -----------------------------------------------------------------

    def run(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Make the plot
        if self.density: self.plot_with_density(output_path)
        else: self.plot_simple(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the scatter plotter ...")

    # -----------------------------------------------------------------

    def plot_simple(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Making the scatter plot ...")

        # Create the figure
        self._figure = plt.figure(figsize=(10, 10))

        # Add first subplot
        ax = self._figure.add_subplot(1, 1, 1, projection='3d')

        ax.scatter(self.x, self.y, self.z)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.y_limits)
        ax.set_zlim(self.z_limits)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_zlabel(self.z_label)

        # Set the title
        if self.title is not None: plt.suptitle("\n".join(wrap(self.title, 60)))

        plt.tight_layout()

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the scatter plot to a buffer ...")
        else: log.debug("Saving the scatter plot to " + str(path) + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25, format=self.format, transparent=self.transparent)
        plt.close()

    # -----------------------------------------------------------------

    def plot_with_density(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Making the scatter plot ...")

        # Create the figure
        self._figure = plt.figure(figsize=(15, 15))

        # Add first subplot
        ax = self._figure.add_subplot(2, 2, 1, projection='3d')

        ax.scatter(self.x, self.y, self.z)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.y_limits)
        ax.set_zlim(self.z_limits)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_zlabel(self.z_label)

        # To draw projected points against the axis planes:
        # ax.plot(self.x, self.z, 'r+', zdir='y', zs=self.y_limits[1])
        # ax.plot(self.y, self.z, 'g+', zdir='x', zs=self.x_limits[0])
        # ax.plot(self.x, self.y, 'k+', zdir='z', zs=self.z_limits[0)

        # Add second subplot
        ax = self._figure.add_subplot(2, 2, 2)

        # Density plot of FUV young vs. FUV ionizing
        x = np.array(self.x)
        y = np.array(self.y)
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.y_limits)

        # Add third subplot
        ax = self._figure.add_subplot(2, 2, 3)

        # Density plot of FUV young vs. dust mass
        x = np.array(self.x)
        y = np.array(self.z)
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.z_label)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.z_limits)

        # Add fourth subplot
        ax = self._figure.add_subplot(2, 2, 4)

        # Density plot of FUV ionizing vs. dust mass
        x = np.array(self.y)
        y = np.array(self.z)
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)
        ax.set_xlabel(self.y_label)
        ax.set_ylabel(self.z_label)
        ax.set_xlim(self.y_limits)
        ax.set_ylim(self.z_limits)

        # Set the title
        if self.title is not None: plt.suptitle("\n".join(wrap(self.title, 60)))

        plt.tight_layout()

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the scatter plot to a buffer ...")
        else: log.debug("Saving the scatter plot to " + str(path) + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25, format=self.format, transparent=self.transparent)
        plt.close()

# -----------------------------------------------------------------
