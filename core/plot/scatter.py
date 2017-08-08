#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.scatter Contains the Scatter3DPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.scatter import Scatter
from ..tools import strings

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["dodgerblue", "r", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class Scatter3DPlotter(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Set the title
        self.title = None

        # The data (Scatter instances)
        self.data = OrderedDict()

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

    def add_data(self, label, scatter):

        """
        This function ...
        :param label:
        :param scatter:
        :return:
        """

        if label in self.data: raise ValueError("Already scatter data with the label '" + label + "'")
        self.data[label] = scatter

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

    def run(self, output_path, title=None):

        """
        This function ...
        :param output_path:
        :param title:
        :return:
        """

        if title is not None: self.title = title

        # Make the plot
        if self.density: self.plot_with_density(output_path)
        else: self.plot_simple(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function clears everything
        :return:
        """

        # Inform the user
        log.info("Clearing the scatter plotter ...")

    # -----------------------------------------------------------------

    def clear_figure(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the scatter plotter figure ...")

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

        # Loop over the different scatter data sets
        for label in self.data:

            x = self.data[label].x(asarray=True)
            y = self.data[label].y(asarray=True)
            z = self.data[label].z(asarray=True)
            ax.scatter(x, y, z, label=label)

        # Set axes limits
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.y_limits)
        ax.set_zlim(self.z_limits)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_zlabel(self.z_label)

        # Set the title
        if self.title is not None: plt.suptitle(strings.split_in_lines(self.title))

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

        # Loop over the different scatter data sets
        for label in self.data:

            x = self.data[label].x(asarray=True)
            y = self.data[label].y(asarray=True)
            z = self.data[label].z(asarray=True)
            ax.scatter(x, y, z)

        # Set axes limits
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

        for label in self.data:

            x = self.data[label].x(asarray=True)
            y = self.data[label].y(asarray=True)

            # Density plot of FUV young vs. FUV ionizing
            if len(x) > 4:

                x = np.array(x)
                y = np.array(y)
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)

            else: ax.scatter([],[])

        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.y_limits)

        # Add third subplot
        ax = self._figure.add_subplot(2, 2, 3)

        for label in self.data:

            x = self.data[label].x(asarray=True)
            z = self.data[label].z(asarray=True)

            # Density plot of FUV young vs. dust mass
            if len(x) > 4:

                x = np.array(x)
                y = np.array(z)
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)

            else: ax.scatter([],[])

        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.z_label)
        ax.set_xlim(self.x_limits)
        ax.set_ylim(self.z_limits)

        # Add fourth subplot
        ax = self._figure.add_subplot(2, 2, 4)

        for label in self.data:

            y = self.data[label].y(asarray=True)
            z = self.data[label].z(asarray=True)

            # Density plot of FUV ionizing vs. dust mass
            if len(y) > 4:

                x = np.array(y)
                y = np.array(z)
                xy = np.vstack([x, y])
                z = gaussian_kde(xy)(xy)

                # Sort the points by density, so that the densest points are plotted last
                idx = z.argsort()
                x, y, z = x[idx], y[idx], z[idx]
                ax.scatter(x, y, c=z, s=100, edgecolor='', cmap=self.color_map)

            else: ax.scatter([],[])

        ax.set_xlabel(self.y_label)
        ax.set_ylabel(self.z_label)
        ax.set_xlim(self.y_limits)
        ax.set_ylim(self.z_limits)

        # Set the title
        if self.title is not None: plt.suptitle(strings.split_in_lines(self.title))

        plt.tight_layout()

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the scatter plot to a buffer ...")
        else: log.debug("Saving the scatter plot to " + str(path) + " ...")

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25, format=self.format, transparent=self.transparent)
        plt.close()

# -----------------------------------------------------------------
