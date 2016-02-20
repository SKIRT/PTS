#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.sed Contains the SEDPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
line_styles = ['-','--','-.',':']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']

# -----------------------------------------------------------------

class SEDPlotter(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        self.models = dict()
        self.observations = dict()

    # -----------------------------------------------------------------

    def add_modeled_sed(self, sed, label):

        """
        This function ...
        :param sed:
        :return:
        """

        # Add the SED to the list
        self.models[label] = sed

    # -----------------------------------------------------------------

    def add_observed_sed(self, sed, label):

        """
        This function ...
        :param sed:
        :return:
        """

        # Add the SED to the list
        self.observations[label] = sed

    # -----------------------------------------------------------------

    def run(self, output_path, input=None):

        """
        This function ...
        :param output_path:
        :param input:
        :return:
        """

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        wavelengths = None

        figsize = (10,6)
        xlim = None
        #ylim = (-5,2)

        # Setup the figure
        figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        min_flux = None
        max_flux = None

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = self.observations[label].wavelengths
            fluxes = self.observations[label].fluxes
            instruments = self.observations[label].instruments
            bands = self.observations[label].bands
            errors = self.observations[label].errors

            # Filter descriptions
            #descriptions = [instrument + " " + band for instrument, band in zip(instruments, bands)]
            descriptions = bands

            # Set colours and markers for the flux points
            color = iter(cm.rainbow(np.linspace(0,1,len(wavelengths))))

            labels = list(set(list(self.observations[label].bands)))

            markers = filled_markers[:len(labels)]

            used_labels = []

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                if fluxes[k] <= 0.0:

                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                if min_flux is None or fluxes[k] < min_flux: min_flux = fluxes[k]
                if max_flux is None or fluxes[k] > max_flux: max_flux = fluxes[k]

                c = next(color)

                if descriptions[k] in labels:

                    # get marker
                    marker = markers[labels.index(descriptions[k])]

                    if descriptions[k] not in used_labels:

                        error_bar = np.array([[np.fabs(np.log10(fluxes[k]) - np.log10(fluxes[k] - errors[k])), np.fabs(np.log10(fluxes[k])-np.log10(fluxes[k] + errors[k]))]]).T
                        used_labels.append(labels[k])

                        ax1.errorbar(wavelengths[k], np.log10(fluxes[k]), yerr=error_bar, fmt=marker, markersize=7, color=c,markeredgecolor='black', ecolor=c, capthick=2)
                        ax1.plot(wavelengths[k], np.log10(fluxes[k]), marker, markersize=7, color=c, markeredgecolor='black', markerfacecolor=c, label=descriptions[k])
                        ax2.errorbar(wavelengths[k], 0., yerr=errors[k]/fluxes[k] * 100., fmt=marker, markersize=7, color=c, markeredgecolor='black', ecolor=c, capthick=2)

                    else:

                        error_bar = np.array([[np.fabs(np.log10(fluxes[k])-np.log10(fluxes[k]-errors[k])), np.fabs(np.log10(fluxes[k]) - np.log10(fluxes[k] + errors[k]))]]).T
                        ax1.errorbar(wavelengths[k], np.log10(fluxes[k]), yerr=error_bar, fmt=marker, markersize=7, color=c, markeredgecolor='black', ecolor=c, capthick=2)
                        ax2.errorbar(wavelengths[k], 0., yerr=errors[k]/fluxes[k] * 100., fmt=marker, markersize=7, color=c, markeredgecolor='black', ecolor=c, capthick=2)

            # Residuals
            model_labels = self.models.keys()
            for j in range(len(model_labels)):

                model_label = model_labels[j]

                log_model = np.log10(self.models[model_label].fluxes)

                f2 = interp1d(self.models[model_label].wavelengths, log_model, kind='cubic')
                ax2.plot(wavelengths, -(fluxes - f2(wavelengths))/fluxes * 100., line_styles[j], color='black', label='model')

        # Add model SEDs
        model_labels = self.models.keys()
        for j in range(len(model_labels)):

            log_model = np.log10(self.models[model_label].fluxes)
            ax1.plot(self.models[model_label].wavelengths, log_model, line_styles[j], color='black', label=label)

            if self.models[model_label].has_errors:

                bottom = []
                top = []
                for j in range(len(self.models[label].errors)):

                    value = self.models[model_label].fluxes[j]
                    bottom.append(value + self.models[model_label].errors[j][0])
                    top.append(value + self.models[model_label].errors[j][1])

                bottom = np.array(bottom)
                top = np.array(top)

                log_bottom = np.log10(bottom)
                log_top = np.log10(top)

                ax1.fill_between(self.models[model_label].wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor='cyan', edgecolor='cyan', interpolate=True, alpha=0.5)
                ax1.plot([], [], color='cyan', linewidth=10, label='spread')

        ax2.axhline(y=0., color='black', ls='-.')
        ax2.set_ylim(-95,95)

        # Set axis limits if requested
        print(min_flux)
        print(max_flux)
        log_min_flux = np.log10(min_flux)
        log_max_flux = np.log10(max_flux)
        ax1.set_ylim((log_min_flux, log_max_flux))

        figure.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)

        # Add axis labels and a legend
        ax2.set_xscale('log')
        ax1.set_xscale('log')
        ax2.set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')
        ax1.set_ylabel(r"Log $F_\nu$$[Jy]$", fontsize='large')
        ax2.set_ylabel(r"Residuals $[\%]$", fontsize='large')

        # Add the legend
        ax1.legend(numpoints=1, loc=4, frameon=True, ncol=2, fontsize=11)

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25)
        #plt.show()
        plt.close()

# -----------------------------------------------------------------
