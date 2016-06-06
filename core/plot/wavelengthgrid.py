#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.wavelengthgrid Contains the WavelengthGridPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np
from textwrap import wrap
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class WavelengthGridPlotter(object):
    
    """
    This class ...
    """

    def __init__(self, title=None):

        """
        This function ...
        :return:
        """

        # Set the title
        self.title = title

        # The wavelength grids
        self.grids = OrderedDict()

        # SEDs to plot together with the wavelength grids
        self.seds = OrderedDict()

        # Emission lines
        self.emission_lines = []

        # The axis limits
        self.min_wavelength = None
        self.max_wavelength = None

        # The figure
        self._figure = None

        # Properties
        self.size = (20, 5)
        self.colormap = "rainbow"  # or "nipy_spectral"
        self.format = None
        self.transparent = False

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def add_wavelength_grid(self, grid, label):

        """
        This function ...
        :param grid:
        :param label:
        :return:
        """

        self.grids[label] = copy.deepcopy(grid)

    # -----------------------------------------------------------------

    def add_wavelength_point(self, label, wavelength):

        """
        This function ...
        :param label:
        :param wavelength:
        :return:
        """

        self.grids[label].add_point(wavelength)

    # -----------------------------------------------------------------

    def add_sed(self, sed, label):

        """
        This function ...
        :param sed:
        :param label:
        :return:
        """

        self.seds[label] = sed

    # -----------------------------------------------------------------

    def add_emission_line(self, line):

        """
        This function ...
        :param line:
        :return:
        """

        self.emission_lines.append(line)

    # -----------------------------------------------------------------

    def run(self, output_path, min_wavelength=0.019, max_wavelength=2050):

        """
        This function ...
        :param output_path:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        # Set the axis limits
        self.min_wavelength = min_wavelength
        self.max_wavelength = max_wavelength

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the wavelength grid plotter ...")

        # Set default values for all attributes
        self.title = None
        self.grids = OrderedDict()
        self.seds = OrderedDict()
        self.emission_lines = []
        self.min_wavelength = None
        self.max_wavelength = None
        self._figure = None
        self.colormap = "rainbow"
        self.format = None
        self.transparent = False

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grid ...")

        # Create the figure
        self._figure = plt.figure(figsize=self.size)
        plt.ylabel('$T_\lambda$', fontsize=28)
        plt.xlabel('$\lambda/\mu$m', fontsize=28)

        # -----------------------------------------------------------------

        plt.subplot2grid((4, 1), (0, 0), rowspan=3)

        # Setup x axis
        plt.xlim(self.min_wavelength, self.max_wavelength)
        plt.xscale('log')
        plt.gca().xaxis.set_ticklabels([])

        # Set up y axis
        plt.ylim(1, 4)
        plt.gca().yaxis.set_visible(False)
        plt.grid('on')

        # Loop over the SEDs
        for label in self.seds:

            # Get the wavelengths and fluxes array
            wavelengths = self.seds[label].wavelengths(asarray=True)
            fluxes = self.seds[label].fluxes(asarray=True)

            nonzero = fluxes != 0
            wavelengths = wavelengths[nonzero]
            fluxes = fluxes[nonzero]
            log_fluxes = 4. + 0.3 * np.log(fluxes / fluxes.max())

            # Plot the SED
            plt.plot(wavelengths, log_fluxes, color='c', lw=0.3)

        # Loop over the wavelength grids
        for label in self.grids:

            # Get the wavelength points
            grid = self.grids[label]
            wavelengths = grid.wavelengths(asarray=True)
            nwavelengths = grid.nwavelengths

            # Plot the grid points
            plt.scatter(wavelengths, [1.2 for _ in wavelengths], 10, marker='.', color='b', linewidths=0)
            t = plt.text(0.021, 1.3, "{} grid points".format(nwavelengths), horizontalalignment='left',
                         fontsize='xx-small', color='b', backgroundcolor='w')
            t.set_bbox(dict(color='w', alpha=0.5, edgecolor='w'))

            # Plot a vertical line for each grid point
            for w in grid: plt.vlines(w, 1, 4, color='b', lw=0.2, alpha=0.2)

            # Plot the deltas
            plt.subplot2grid((4, 1), (3, 0), rowspan=1)
            plt.plot(grid[1:], np.log10(grid[1:]) - np.log10(grid[:-1]), 'm-', lw=0.6)
            plt.xlim(0.019, 2050)
            plt.xscale('log')
            plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%g"))
            plt.xlabel(r"$\lambda\,(\mu \mathrm{m})$", fontsize='large')
            plt.ylim(0, 0.05)
            plt.ylabel(r"$\Delta\lambda\,(\mathrm{dex})$", fontsize='large')
            plt.grid('on')

        # plot a labeled vertical line for each emission line
        colors = ['m', 'r', 'y', 'g']
        positions = [1.3, 1.45, 1.6, 1.75]
        index = 0

        # Loop over the emission lines
        for label in self.emission_lines:

            # Get center wavelength and label
            emission_line = self.emission_lines[label]
            center = emission_line.center
            label = emission_line.label

            if len(label) > 0:

                plt.vlines(center, 1, 2, color=colors[index], lw=0.3)
                t = plt.text(center, positions[index], label, horizontalalignment='center',
                             fontsize='xx-small', color=colors[index], backgroundcolor='w')
                t.set_bbox(dict(color='w', alpha=0.5, edgecolor='w'))
                index = (index + 1) % 4

        # -----------------------------------------------------------------

        # Set the title
        if self.title is not None: self._figure.suptitle("\n".join(wrap(self.title, 60)))

        # Finish
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def finish_plot(self, path):

        """
        This function ...
        :return:
        """

        # Debugging
        if type(path).__name__ == "BytesIO": log.debug("Saving the SED plot to a buffer ...")
        elif path is None: log.debug("Showing the SED plot ...")
        else: log.debug("Saving the SED plot to " + str(path) + " ...")

        if path is not None:
            # Save the figure
            plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------
