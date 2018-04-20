#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.attenuation Contains the AttenuationPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ...magic.tools.plotting import pretty_colours, line_styles, filled_markers

# -----------------------------------------------------------------

def plot_attenuation_curve(curve, label):

    """
    This function ...
    :param curve:
    :param label:
    :return:
    """

    # Create the plotter
    plotter = AttenuationPlotter()

    # Add
    plotter.add_attenuation_curve(curve, label)

    # Run
    plotter.run()

# -----------------------------------------------------------------

def plot_attenuation_curves(curves):

    """
    This function ...
    :param curves:
    :return:
    """

    # Create the plotter
    plotter = AttenuationPlotter()

    # Add
    plotter.add_attenuation_curves(curves)

    # Run
    plotter.run()

# -----------------------------------------------------------------

class AttenuationPlotter(object):
    
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

        # The axis limits
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_attenuation = None
        self.max_attenuation = None

        # Properties
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

    def add_attenuation_curve(self, attenuation_curve, label):

        """
        This function ...
        :param attenuation_curve:
        :param label:
        :return:
        """

        self.curves[label] = attenuation_curve

    # -----------------------------------------------------------------

    def add_attenuation_curves(self, curves):

        """
        This function ...
        :param curves:
        :return:
        """

        for label in curves: self.add_attenuation_curve(curves[label], label)

    # -----------------------------------------------------------------

    def run(self, output_path=None, min_wavelength=None, max_wavelength=None, min_attenuation=None, max_attenuation=None):

        """
        This function ...
        :param output_path:
        :param min_wavelength:
        :param max_wavelength:
        :param min_attenuation:
        :param max_attenuation:
        :return:
        """

        # Set the axis limits
        self.min_wavelength = min_wavelength
        self.max_wavelength = max_wavelength
        self.min_attenuation = min_attenuation
        self.max_attenuation = max_attenuation

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the attenuation plotter ...")

        # Set default values for all attributes
        self.title = None
        self.curves = OrderedDict()
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_attenuation = None
        self.max_attenuation = None
        self.transparent = False

    # -----------------------------------------------------------------

    def plot(self, path=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Making the attenuation plot ...")

        # Plot the attenuation curves
        plt.figure(figsize=(10, 10))
        plt.ylabel('$A_\lambda/A_V$', fontsize=28)
        plt.xlabel('$\lambda/\mu$m', fontsize=28)
        plt.xlim(0.1, 2)
        plt.ylim(0, 8)
        #plt.ylim(0, 1)
        plt.xscale('log')
        x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2]
        plt.xticks(x, x)
        # plt.yscale('log')
        plt.tick_params(labelsize=17)
        # plt.subplots_adjust(bottom=0.1)

        #plt.plot(wl_MW, n_atts_MW, 'r-', label="MW", linewidth=2)
        #plt.plot(wl_SMC, n_atts_SMC, 'b--', label="SMC", linewidth=2)
        #plt.plot(wl_B16, Qfit_B16 + 1, 'g-.', label="Battisti+16", linewidth=3)
        ## plt.plot(wl_Calzetti,n_atts_Calzetti, 'g-', label="Cal+94",linewidth=2)

        #plt.plot(modWls, n_atts_tot, 'k-', label="M31 tot", linewidth=2)
        #plt.plot(modWls, n_atts_diff, 'k--', label="M31 diffuse", linewidth=2)
        #plt.plot(modWls, n_atts_mappings, 'k:', label="M31 SF regions", linewidth=2)

        colors = iter(pretty_colours)

        for label in self.curves:

            curve = self.curves[label]

            wavelengths = curve.wavelengths(unit="micron", add_unit=False)
            attenuations = curve.attenuations()

            # Plot the curve
            plt.plot(wavelengths, attenuations, label=label, linewidth=2, color=next(colors))

        #plt.errorbar(wl_M31, n_atts_M31, yerr=err_att_M31, color='m', fmt='ko', markersize=8, label="M31- Dong+14")
        ## plt.plot(1./wl_MWRv5,n_atts_MWRv5, 'm-', label="MW Rv=5")
        ## plt.plot([1./0.55,1./0.55],[0,15], 'k-')

        plt.tight_layout()
        plt.legend(loc='upper right', numpoints=1, markerscale=1.5, fontsize=24)

        # Debugging
        if path is not None: log.debug("Saving the attenuation plot to " + path + " ...")

        # Save the figure
        if path is not None: plt.savefig(path, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent)
        else: plt.show()
        plt.close()

# -----------------------------------------------------------------
