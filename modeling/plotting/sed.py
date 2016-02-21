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
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
line_styles = ['-','--','-.',':']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'] * 6

distinguishable_colormaps = ["spring", "winter", "copper"]

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
        :param label:
        :return:
        """

        # Add the SED to the list
        self.models[label] = sed

    # -----------------------------------------------------------------

    def add_observed_sed(self, sed, label):

        """
        This function ...
        :param sed:
        :param label:
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

        if len(self.models) == 0: self.plot_no_models(path)
        else: self.plot_with_models(path)

        return

        wavelengths = None

        figsize = (10,6)

        # Setup the figure
        figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        ax1 = plt.subplot(gs[0])
        #ax2 = plt.subplot(gs[1], sharex=ax1)

        min_flux = None
        max_flux = None



        # Make iterable from color map names
        color_maps = iter(distinguishable_colormaps)

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = self.observations[label].wavelengths
            fluxes = self.observations[label].fluxes
            instruments = self.observations[label].instruments
            bands = self.observations[label].bands
            errors = self.observations[label].errors

            # Filter descriptions
            descriptions = []
            for j in range(len(instruments)):

                masked = False
                try: masked = instruments.mask[j]
                except AttributeError: pass

                if masked: descriptions.append(bands[j])
                else: descriptions.append(instruments[j])

            # Determine color map class
            if len(self.observations) == 1: colormap = plt.get_cmap("rainbow")
            else: colormap = plt.get_cmap(next(color_maps))

            # Create color range
            color_range = iter(colormap(np.linspace(0,1,len(wavelengths))))

            labels = list(set(list(descriptions)))

            markers = filled_markers[:len(labels)]

            used_labels = []

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                if fluxes[k] <= 0.0:

                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                if min_flux is None or fluxes[k] < min_flux: min_flux = fluxes[k]
                if max_flux is None or fluxes[k] > max_flux: max_flux = fluxes[k]

                color = next(color_range)

                if descriptions[k] in labels:

                    # get marker
                    marker = markers[labels.index(descriptions[k])]

                    if descriptions[k] not in used_labels:

                        error_bar = np.array([[np.fabs(np.log10(fluxes[k]) - np.log10(fluxes[k] - errors[k])), np.fabs(np.log10(fluxes[k])-np.log10(fluxes[k] + errors[k]))]]).T
                        used_labels.append(descriptions[k])

                        ax1.errorbar(wavelengths[k], np.log10(fluxes[k]), yerr=error_bar, fmt=marker, markersize=7, color=color,markeredgecolor='black', ecolor=color, capthick=2)
                        ax1.plot(wavelengths[k], np.log10(fluxes[k]), marker, markersize=7, color=color, markeredgecolor='black', markerfacecolor=color, label=descriptions[k])
                        ax2.errorbar(wavelengths[k], 0., yerr=errors[k]/fluxes[k] * 100., fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

                    else:

                        error_bar = np.array([[np.fabs(np.log10(fluxes[k])-np.log10(fluxes[k]-errors[k])), np.fabs(np.log10(fluxes[k]) - np.log10(fluxes[k] + errors[k]))]]).T
                        ax1.errorbar(wavelengths[k], np.log10(fluxes[k]), yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)
                        ax2.errorbar(wavelengths[k], 0., yerr=errors[k]/fluxes[k] * 100., fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

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
        log_min_flux = np.log10(min_flux)
        log_max_flux = np.log10(max_flux)
        plot_min = 0.9 * log_min_flux if log_min_flux > 0 else 1.1 * log_min_flux
        plot_max = 1.1 * log_max_flux
        ax1.set_ylim((plot_min, plot_max))

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
        plt.close()

    # -----------------------------------------------------------------

    def plot_no_models(self, path):

        """
        This function ...
        :return:
        """

        if len(self.observations) == 1: self.plot_one_observation(path)
        else: self.plot_more_observations(path)

    # -----------------------------------------------------------------

    def plot_one_observation(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        figsize = (10,6)

        # Setup the figure
        figure = plt.figure(figsize=figsize)

        # Keep track of the minimum and maximum flux encountered
        min_flux = None
        max_flux = None

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths
        fluxes = observation.fluxes
        instruments = observation.instruments
        bands = observation.bands
        errors = observation.errors

        # Get labels and descriptions
        labels, descriptions = get_labels_and_descriptions(instruments, bands)

        # Determine color map class
        colormap = plt.get_cmap("rainbow")

        # Create color range
        color_range = iter(colormap(np.linspace(0, 1, len(wavelengths))))

        # Unique labels
        unique_labels = list(set(list(labels)))

        # Markers for the unique labels
        markers = filled_markers[:len(unique_labels)]

        # Used labels
        used_labels = []

        # Loop over the wavelengths
        for k in range(len(wavelengths)):

            # Check validity of flux value
            if fluxes[k] <= 0.0:
                log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                continue

            # Keep track of minimum and maximum flux
            if min_flux is None or fluxes[k] < min_flux: min_flux = fluxes[k]
            if max_flux is None or fluxes[k] > max_flux: max_flux = fluxes[k]

            # Get next color
            color = next(color_range)

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            plot_wavelength(plt.gca(), labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

        # Set axis limits if requested
        log_min_flux = np.log10(min_flux)
        log_max_flux = np.log10(max_flux)
        plot_min = 0.9 * log_min_flux if log_min_flux > 0 else 1.1 * log_min_flux
        plot_max = 1.1 * log_max_flux
        plt.gca().set_ylim((plot_min, plot_max))

        # Hide ticks
        figure.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)

        # Add axis labels and a legend
        plt.gca().set_xscale('log')
        plt.gca().set_ylabel(r"Log $F_\nu$$[Jy]$", fontsize='large')

        # Add the legend
        plt.gca().legend(numpoints=1, loc=4, frameon=True, ncol=2, fontsize=11)

        # Save the figure
        plt.savefig(path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    def plot_more_observations(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        figsize = (10,6)

        # http://matplotlib.org/examples/color/colormaps_reference.html
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))

        # For each observed SED (different color map name):
        # Make a small new axis (like a legend)
        #    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))

        # Make smaller plot inside other axis:
        # http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
        #ax_number = add_subplot_axes(ax,rect)

        # Setup the figure
        figure = plt.figure(figsize=figsize)

        number_of_observations = len(self.observations)

        height_ratios = [4,1]
        for _ in range(number_of_observations): height_ratios.append(0.2)

        # Create axes
        gs = gridspec.GridSpec(2+number_of_observations, 1, height_ratios=height_ratios)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)
        colormap_axes = []
        for i in range(number_of_observations):
            ax = plt.subplot(gs[2+i])
            colormap_axes.append(ax)
        colormap_axes = iter(colormap_axes)

        # Keep track of the minimum and maximum flux encountered
        min_flux = None
        max_flux = None

        # Keep track of whether we are plotting the first observed SED
        first = True
        reference_sed = None

        # Make iterable from color map names
        color_maps = iter(distinguishable_colormaps)

        # Unique labels
        unique_labels = set()

        # Get all unique labels
        for label in self.observations:

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            instruments = observation.instruments
            bands = observation.bands

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            # Add labels
            for _label in labels: unique_labels.add(_label)

        # Convert unique labels (set) to list
        unique_labels = list(unique_labels)

        # Markers for the unique labels
        markers = filled_markers[:len(unique_labels)]

        # Used labels
        used_labels = []

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths
            fluxes = observation.fluxes
            instruments = observation.instruments
            bands = observation.bands
            errors = observation.errors

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            # Determine color map class
            colormap = plt.get_cmap(next(color_maps))

            ax3 = next(colormap_axes)
            ax3.imshow(gradient, aspect='auto', cmap=colormap)

            # Create color range
            color_range = iter(colormap(np.linspace(0, 1, len(wavelengths))))

            # Reference SED for comparison with other observed SEDs
            if first: reference_sed = observation

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                # Check validity of flux value
                if fluxes[k] <= 0.0:
                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                # Keep track of minimum and maximum flux
                if min_flux is None or fluxes[k] < min_flux: min_flux = fluxes[k]
                if max_flux is None or fluxes[k] > max_flux: max_flux = fluxes[k]

                # Get next color
                color = next(color_range)

                # Get marker
                marker = markers[unique_labels.index(labels[k])]

                # Plot the flux data point on the main axis with the specified marker and color
                plot_wavelength(ax1, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

                # Plot on axis 2
                if first:

                    value = 0.0
                    error = errors[k] / fluxes[k] * 100.

                else:

                    reference_flux = find_reference_flux(instruments[k], bands[k], wavelengths[k], reference_sed)

                    if reference_flux is None:

                        value = None
                        error = None

                    else:

                        value = (fluxes[k] - reference_flux) / reference_flux * 100.
                        error = errors[k] / reference_flux * 100.0

                if value is not None: ax2.errorbar(wavelengths[k], value, yerr=error, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

            # The next observation is not the first anymore
            first = False

        # Set linestyle and limit for axis2
        ax2.axhline(y=0., color='black', ls='-.')
        ax2.set_ylim(-95, 95)

        # Set axis limits if requested
        log_min_flux = np.log10(min_flux)
        log_max_flux = np.log10(max_flux)
        plot_min = 0.9 * log_min_flux if log_min_flux > 0 else 1.1 * log_min_flux
        plot_max = 1.1 * log_max_flux
        ax1.set_ylim((plot_min, plot_max))

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
        plt.close()

    # -----------------------------------------------------------------

    def plot_with_models(self, path):

        """
        This function ...
        :return:
        """

        figsize = (10,6)

# -----------------------------------------------------------------

def get_labels_and_descriptions(instruments, bands):

    """
    This function ...
    :param instruments:
    :param bands:
    :return:
    """

    # Filter descriptions
    labels = []
    descriptions = []
    for j in range(len(instruments)):

        masked = False
        try: masked = instruments.mask[j]
        except AttributeError: pass

        if masked:
            labels.append(bands[j])
            descriptions.append(bands[j])
        else:
            labels.append(instruments[j])
            descriptions.append(instruments[j] + " " + bands[j])

    return labels, descriptions

# -----------------------------------------------------------------

def plot_wavelength(axis, label, used_labels, wavelength, flux, error, marker, color):

    """
    This function ...
    :param axis:
    :param label:
    :param used_labels:
    :param wavelength:
    :param flux:
    :param error:
    :param marker:
    :param color:
    :return:
    """

    if label not in used_labels:

        error_bar = np.array([[np.fabs(np.log10(flux) - np.log10(flux - error)), np.fabs(np.log10(flux) - np.log10(flux + error))]]).T
        used_labels.append(label)

        axis.errorbar(wavelength, np.log10(flux), yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)
        axis.plot(wavelength, np.log10(flux), marker, markersize=7, color=color, markeredgecolor='black', markerfacecolor=color, label=label)

    else:

        error_bar = np.array([[np.fabs(np.log10(flux)-np.log10(flux - error)), np.fabs(np.log10(flux) - np.log10(flux + error))]]).T
        axis.errorbar(wavelength, np.log10(flux), yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

# -----------------------------------------------------------------

def find_reference_flux(instrument, band, wavelength, reference):

    """
    This function ...
    :param instrument:
    :param band:
    :param wavelength:
    :param reference:
    :return:
    """

    # Loop over the entries in the reference SED
    for i in range(len(reference.wavelengths)):

        if (reference.instruments[i] == instrument and reference.bands[i] == band) or np.isclose(reference.wavelengths[i], wavelength):
            return reference.fluxes[i]

    # No match
    return None

# -----------------------------------------------------------------

def add_subplot_axes(ax, rect, axisbg='w'):

    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

# -----------------------------------------------------------------
