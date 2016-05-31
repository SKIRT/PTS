#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.transmission Contains the TransmissionPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from textwrap import wrap
import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class TransmissionPlotter(object):
    
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

        # The different curves
        self.curves = OrderedDict()

        # The wavelengths
        self.wavelengths = []

        # The axis limits
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_transmission = None
        self.max_transmission = None

        # The figure
        self._figure = None

        # Properties
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

    def add_transmission_curve(self, transmission_curve, label):

        """
        This function ...
        :param transmission_curve:
        :param label:
        :return:
        """

        # Add the transmission curve
        self.curves[label] = transmission_curve

    # -----------------------------------------------------------------

    def add_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        # Add the wavelength
        self.wavelengths.append(wavelength)

    # -----------------------------------------------------------------

    def run(self, output_path, min_wavelength=None, max_wavelength=None, min_transmission=None, max_transmission=None):

        """
        This function ...
        :param output_path:
        :param min_wavelength:
        :param max_wavelength:
        :param min_transmission:
        :param max_transmission:
        :return:
        """

        # Set the axis limits
        self.min_wavelength = min_wavelength
        self.max_wavelength = max_wavelength
        self.min_transmission = min_transmission
        self.max_transmission = max_transmission

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the transmission plotter ...")

        # Set default values for all attributes
        self.title = None
        self.curves = OrderedDict()
        self.wavelengths = []
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_transmission = None
        self.max_transmission = None
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
        log.info("Plotting the transmission plot ...")

        self._figure = plt.figure(figsize=(17, 5))
        plt.ylabel('$T_\lambda$', fontsize=28)
        plt.xlabel('$\lambda/\mu$m', fontsize=28)

        # Set axes limits
        plt.xlim(self.min_wavelength, self.max_wavelength)
        plt.ylim(self.min_transmission, self.max_transmission)

        plt.xscale('log')

        plt.tick_params(labelsize=17)

        # Get the color map
        cm = plt.get_cmap(self.colormap)

        cNorm = colors.Normalize(vmin=0, vmax=len(self.curves) - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        # Sort the labels based on the peak wavelength
        sorted_labels = sorted(self.curves.keys(), key=lambda label: self.curves[label].peak_wavelength)

        # Plot the transmission curves
        counter = 0
        for label in sorted_labels:

            curve = self.curves[label]

            wavelengths = curve.wavelengths(unit="micron", add_unit=False)
            transmissions = curve.transmissions()

            colorVal = scalarMap.to_rgba(counter)

            # Plot the curve
            #plt.plot(wavelengths, transmissions, label=label, linewidth=2, color=colorVal)
            plt.fill(wavelengths, transmissions, label=label, linewidth=2, color=colorVal, alpha=0.5)

            counter += 1

        # Plot the wavelengths
        for wavelength in self.wavelengths:
            plt.axvline(wavelength.to("micron").value, color="0.8")

        # Set the title
        if self.title is not None: self._figure.suptitle("\n".join(wrap(self.title, 60)))

        # plt.tight_layout()

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
