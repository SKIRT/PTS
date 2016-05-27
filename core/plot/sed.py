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
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools.logging import log

# -----------------------------------------------------------------

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
line_styles = ['-', '--', '-.', ':']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'] * 6

distinguishable_colormaps = ["spring", "winter", "copper", "cool", "PRGn", "coolwarm"]

other_colormaps = ["YlGn", "YlGnBu", "YlOrRd", "Purples"]

# http://matplotlib.org/examples/color/named_colors.html
color_hex = colors.cnames

# Add the single letter colors
for name, rgb in colors.ColorConverter.colors.items():
    hex_ = colors.rgb2hex(rgb)
    color_hex[name] = hex_

color_plt_identifiers = color_hex.keys()

pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class SEDPlotter(object):
    
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

        # Create ordered dictionaries for the model and observed SEDs (the order of adding SEDs is remembered)
        self.models = OrderedDict()
        self.models_no_residuals = OrderedDict()
        self.observations = OrderedDict()

        # Keep track of the minimal and maximal wavelength and flux encountered during the plotting
        self._min_wavelength = None
        self._max_wavelength = None
        self._min_flux = None
        self._max_flux = None

        # The definite axes limits
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_flux = None
        self.max_flux = None

        # Store the figure and its axes as references
        self._figure = None
        self._main_axis = None
        self._residual_axes = []

        # Properties
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

    def add_modeled_sed(self, sed, label, residuals=True):

        """
        This function ...
        :param sed:
        :param label:
        :param residuals: whether plotting the residual curve makes sense for this SED
        (it does not when the SED is only for one component, or a contribution to the total SED)
        :return:
        """

        # Add the SED to the appropriate dictionary
        if residuals: self.models[label] = sed
        else: self.models_no_residuals[label] = sed

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

    def run(self, output_path, min_wavelength=None, max_wavelength=None, min_flux=None, max_flux=None):

        """
        This function ...
        :param output_path:
        :param min_wavelength:
        :param max_wavelength:
        :param min_flux:
        :param max_flux:
        :return:
        """

        # Set the axis limits
        self.min_wavelength = min_wavelength
        self.max_wavelength = max_wavelength
        self.min_flux = min_flux
        self.max_flux = max_flux

        # Make the plot
        self.plot(output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the SED plotter ...")

        # Set default values for all attributes
        self.models = OrderedDict()
        self.models_no_residuals = OrderedDict()
        self.observations = OrderedDict()
        self._min_wavelength = None
        self._max_wavelength = None
        self._min_flux = None
        self._max_flux = None
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_flux = None
        self.max_flux = None
        self._figure = None
        self._main_axis = None
        self._residual_axes = []

    # -----------------------------------------------------------------

    def plot(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Making the SED plot ...")

        if len(self.models) == 0 and len(self.models_no_residuals) == 0: self.plot_no_models(path)
        else: self.plot_with_models(path)

    # -----------------------------------------------------------------

    def plot_no_models(self, path):

        """
        This function ...
        :param path:
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

        # Determine color map class
        colormap = plt.get_cmap("rainbow")

        # The size of the figure
        figsize = (10,6)

        # Setup the figure
        self._figure = plt.figure(figsize=figsize)
        self._main_axis = self._figure.gca()

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit="micron", add_unit=False)
        fluxes = observation.fluxes(unit="Jy", add_unit=False)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit="Jy", add_unit=False)

        # Create colors
        colors = colormap(np.linspace(0, 1, len(wavelengths)))

        # Plot
        self.draw_observation(instruments, bands, wavelengths, fluxes, errors, colors)

        # Finish the plot
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def plot_more_observations(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # The size of the figure
        figsize = (10,6)

        # http://matplotlib.org/examples/color/colormaps_reference.html
        #gradient = np.linspace(0, 1, 256)
        #gradient = np.vstack((gradient, gradient))

        # For each observed SED (different color map name):
        # Make a small new axis (like a legend)
        #    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))

        # Make smaller plot inside other axis:
        # http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
        #ax_number = add_subplot_axes(ax,rect)

        # Setup the figure
        self._figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        self._main_axis = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=self._main_axis)
        self._residual_axes.append(ax2)

        # Count the number of observed SEDs
        number_of_observations = len(self.observations)

        # Keep track of whether we are plotting the first observed SED
        first = True
        reference_sed = None

        # Make iterable from color map names
        color_maps = iter(distinguishable_colormaps+other_colormaps)

        # Make iterable from distinct colors
        different_colors = iter(pretty_colors)

        # Unique labels
        unique_labels = get_unique_labels(self.observations)

        # Markers for the unique labels
        markers = filled_markers[:len(unique_labels)]

        # Used labels
        used_labels = []

        legend_patches = []
        legend_labels = []

        legend_rectangles = []
        rectangle_labels = []

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit="micron", add_unit=False)
            fluxes = observation.fluxes(unit="Jy", add_unit=False)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit="Jy", add_unit=False)

            # Create color range
            if number_of_observations <= 3:
                colormap = plt.get_cmap(next(color_maps))
                colors = colormap(np.linspace(0, 1, len(wavelengths)))
            else: colors = [color_hex[next(different_colors)]] * len(wavelengths)

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            # Reference SED for comparison with other observed SEDs
            if first: reference_sed = observation

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                # Check validity of flux value
                if fluxes[k] <= 0.0:
                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                # Get marker
                marker = markers[unique_labels.index(labels[k])]

                # Plot the flux data point on the main axis with the specified marker and color
                patch = self.plot_wavelength(self._main_axis, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, colors[k], return_patch=True)

                if patch is not None:
                    legend_patches.append(patch)
                    legend_labels.append(labels[k])

                # Plot on axis 2
                if first:

                    value = 0.0
                    error = errors[k] / fluxes[k] * 100.

                else:

                    reference_flux = find_reference_flux(instruments[k], bands[k], wavelengths[k], reference_sed.instruments(), reference_sed.bands(), reference_sed.wavelengths(unit="micron", add_unit=False), reference_sed.fluxes(unit="Jy", add_unit=False))

                    if reference_flux is None:

                        value = None
                        error = None

                    else:

                        value = (fluxes[k] - reference_flux) / reference_flux * 100.
                        error = errors[k] / reference_flux * 100.0

                if value is not None:

                    error_bar = np.array([[abs(error.lower), abs(error.upper)]]).T
                    ax2.errorbar(wavelengths[k], value, yerr=error_bar, fmt=marker, markersize=7, color=colors[k], markeredgecolor='black', ecolor=colors[k], capthick=2)

            # The next observation is not the first anymore
            first = False

            # If color gradient is used for observations instead of single colors ...
            #ax3 = next(colormap_axes)
            #ax3.imshow(gradient, aspect='auto', cmap=colormap)

            # Create rectangle for
            rectangle = patches.Rectangle((0, 0), 1, 1, fc=colors[0])
            legend_rectangles.append(rectangle)
            rectangle_labels.append(label)

        # Extra legend: the different observations
        # fancybox=True makes the legend corners rounded
        observations_legend = self._main_axis.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Finish the plot
        self.finish_plot(path, for_legend_patches=legend_patches, for_legend_parameters=legend_labels, extra_legend=observations_legend)

    # -----------------------------------------------------------------

    def plot_with_models(self, path):

        """
        This function ...
        :return:
        """

        if len(self.observations) == 0: self.plot_only_models(path)
        elif len(self.observations) == 1: self.plot_one_observation_with_models(path)
        else: self.plot_more_observations_with_models(path)

    # -----------------------------------------------------------------

    def plot_only_models(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # The size of the figure
        figsize = (10, 6)

        # Setup the figure
        self._figure = plt.figure(figsize=figsize)
        self._main_axis = self._figure.gca()

        # Add model SEDs
        model_labels = self.models.keys()
        for j in range(len(model_labels)):

            model_label = model_labels[j]

            # Get fluxes, wavelengths and errors
            fluxes = self.models[model_label].fluxes(unit="Jy", add_unit=False)
            wavelengths = self.models[model_label].wavelengths(unit="micron", add_unit=False)
            errors = self.models[model_label].errors(unit="Jy", add_unit=False) if self.models[
                model_label].has_errors else None

            # Plot the model SED as a line (with errors if present)
            self.draw_model(self._main_axis, wavelengths, fluxes, line_styles[j], model_label, errors=errors)

        if len(self.models_no_residuals) > 0:

            line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
            line_styles_models_no_residuals = ["-"] * len(self.models_no_residuals)

            # Add model SEDs in the models_no_residuals dictionary
            model_labels_no_residuals = self.models_no_residuals.keys()
            for j in range(len(model_labels_no_residuals)):

                model_label = model_labels_no_residuals[j]

                # Get fluxes, wavelengths and errors
                fluxes = self.models_no_residuals[model_label].fluxes(unit="Jy", add_unit=False)
                wavelengths = self.models_no_residuals[model_label].wavelengths(unit="micron", add_unit=False)
                errors = self.models_no_residuals[model_label].errors(unit="Jy", add_unit=False) if self.models_no_residuals[model_label].has_errors else None

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self._main_axis, wavelengths, fluxes, line_styles_models_no_residuals[j], model_label, errors=errors, linecolor=line_colors_models_no_residuals[j], adjust_extrema=False)

        # Finish the plot
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def plot_one_observation_with_models(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # The size of the figure
        figsize = (10,6)

        # Setup the figure
        self._figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        self._main_axis = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=self._main_axis)
        self._residual_axes.append(ax2)

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit="micron", add_unit=False)
        fluxes = observation.fluxes(unit="Jy", add_unit=False)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit="Jy", add_unit=False)

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

            # Get next color
            color = next(color_range)

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            self.plot_wavelength(self._main_axis, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

            # Plot point at y=0.0 with errorbar on axis 2
            value = 0.0
            error = errors[k] / fluxes[k] * 100.
            error_bar = np.array([[abs(error.lower), abs(error.upper)]]).T
            ax2.errorbar(wavelengths[k], value, yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

        # Model labels
        model_labels = self.models.keys()

        line_styles_models = line_styles
        line_colors_models = ['black'] * len(model_labels)

        # Residuals
        for j in range(len(model_labels)):

            model_label = model_labels[j]

            #log_model = np.log10(self.models[model_label].fluxes(unit="Jy", add_unit=False))

            model_fluxes= self.models[model_label].fluxes(unit="Jy", add_unit=False)
            f2 = interp1d(self.models[model_label].wavelengths(unit="micron", add_unit=False), model_fluxes, kind='cubic')
            residuals = -(fluxes - f2(wavelengths))/fluxes * 100.

            ax2.plot(wavelengths, residuals, line_styles_models[j], color=line_colors_models[j], label='model')

        # Add model SEDs
        for j in range(len(model_labels)):

            # Get the current model label
            model_label = model_labels[j]

            log_model = np.log10(self.models[model_label].fluxes(unit="Jy", add_unit=False))
            self._main_axis.plot(self.models[model_label].wavelengths(unit="micron", add_unit=False), log_model, line_styles_models[j], color=line_colors_models[j], label=model_label)

            if self.models[model_label].has_errors:

                bottom = []
                top = []
                for j in range(len(self.models[model_label].errors)):

                    value = self.models[model_label].fluxes[j]
                    bottom.append(value + self.models[model_label].errors[j][0])
                    top.append(value + self.models[model_label].errors[j][1])

                bottom = np.array(bottom)
                top = np.array(top)

                log_bottom = np.log10(bottom)
                log_top = np.log10(top)

                self._main_axis.fill_between(self.models[model_label].wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor='cyan', edgecolor='cyan', interpolate=True, alpha=0.5)
                self._main_axis.plot([], [], color='cyan', linewidth=10, label='spread')

        if len(self.models_no_residuals) > 0:

            line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
            line_styles_models_no_residuals = ["-"] * len(self.models_no_residuals)

            # Add model SEDs from the models_no_residuals dictionary
            model_labels_no_residuals = self.models_no_residuals.keys()
            for j in range(len(model_labels_no_residuals)):

                # Get the current model label
                model_label = model_labels_no_residuals[j]

                log_model = np.log10(self.models_no_residuals[model_label].fluxes(unit="Jy", add_unit=False))
                self._main_axis.plot(self.models_no_residuals[model_label].wavelengths(unit="micron", add_unit=False), log_model,
                                     line_styles_models_no_residuals[j], color=line_colors_models_no_residuals[j], label=model_label)

                if self.models_no_residuals[model_label].has_errors:

                    bottom = []
                    top = []
                    for j in range(len(self.models_no_residuals[model_label].errors)):
                        value = self.models_no_residuals[model_label].fluxes[j]
                        bottom.append(value + self.models_no_residuals[model_label].errors[j][0])
                        top.append(value + self.models_no_residuals[model_label].errors[j][1])

                    bottom = np.array(bottom)
                    top = np.array(top)

                    log_bottom = np.log10(bottom)
                    log_top = np.log10(top)

                    self._main_axis.fill_between(self.models_no_residuals[model_label].wavelengths, log_bottom, log_top,
                                                 where=log_top <= log_bottom, facecolor='cyan', edgecolor='cyan',
                                                 interpolate=True, alpha=0.5)
                    self._main_axis.plot([], [], color='cyan', linewidth=10, label='spread')

        # Finish the plot
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def plot_more_observations_with_models(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Count the number of observed SEDs
        number_of_observations = len(self.observations)

        figsize = [10,6]
        figsize[1] += (number_of_observations - 1) * 1.1

        height_ratios = [4]
        for _ in range(number_of_observations): height_ratios.append(1)

        # Setup the figure
        self._figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(1 + number_of_observations, 1, height_ratios=height_ratios)
        self._main_axis = plt.subplot(gs[0])
        for i in range(number_of_observations):
            ax = plt.subplot(gs[1+i], sharex=self._main_axis)
            self._residual_axes.append(ax)
        residual_axes = iter(self._residual_axes)

        # Make iterable from color map names
        color_maps = iter(distinguishable_colormaps+other_colormaps)

        # Make iterable from distinct colors
        colors = iter(pretty_colors)

        # Get unique labels
        unique_labels = get_unique_labels(self.observations)

        # Markers for the unique labels
        markers = filled_markers[:len(unique_labels)]

        # Used labels
        used_labels = []

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get the axis of the next residual plot
            ax2 = next(residual_axes)

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit="micron", add_unit=False)
            fluxes = observation.fluxes(unit="Jy", add_unit=False)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit="Jy", add_unit=False)

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            # Determine color map class
            colormap = plt.get_cmap(next(color_maps))

            # Create color range
            if number_of_observations <= 3: color_range = iter(colormap(np.linspace(0, 1, len(wavelengths))))
            else: color_range = iter([color_hex[next(colors)]] * len(wavelengths))

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                # Check validity of flux value
                if fluxes[k] <= 0.0:
                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                # Get next color
                color = next(color_range)

                # Get marker
                marker = markers[unique_labels.index(labels[k])]

                # Plot the flux data point on the main axis with the specified marker and color
                self.plot_wavelength(self._main_axis, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

                # Plot measurement points on residual plot
                value = 0.0
                error = errors[k] / fluxes[k] * 100.
                error_bar = np.array([[abs(error.lower), abs(error.upper)]]).T
                ax2.errorbar(wavelengths[k], value, yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

            # Residuals
            model_labels = self.models.keys()
            for j in range(len(model_labels)):

                model_label = model_labels[j]

                log_model = np.log10(self.models[model_label].fluxes(unit="Jy", add_unit=False))

                f2 = interp1d(self.models[model_label].wavelengths(unit="micron", add_unit=False), log_model, kind='cubic')
                ax2.plot(wavelengths, -(fluxes - f2(wavelengths))/fluxes * 100., line_styles[j], color='black', label='model')

        # Add model SEDs
        model_labels = self.models.keys()
        for j in range(len(model_labels)):

            model_label = model_labels[j]

            # Get fluxes, wavelengths and errors
            fluxes = self.models[model_label].fluxes(unit="Jy", add_unit=False)
            wavelengths = self.models[model_label].wavelengths(unit="micron", add_unit=False)
            errors = self.models[model_label].errors(unit="Jy", add_unit=False) if self.models[model_label].has_errors else None

            # Plot the model SED as a line (with errors if present)
            self.draw_model(self._main_axis, wavelengths, fluxes, line_styles[j], model_label, errors=errors)

        # Finish the plot
        self.finish_plot(path)

    # -----------------------------------------------------------------

    def plot_wavelength(self, axis, label, used_labels, wavelength, flux, error, marker, color, return_patch=False):

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

        # Keep track of minimum and maximum flux
        if self._min_flux is None or flux < self._min_flux: self._min_flux = flux
        if self._max_flux is None or flux > self._max_flux: self._max_flux = flux

        # Keep track of the minimal and maximal wavelength
        if self._min_wavelength is None or wavelength < self._min_wavelength: self._min_wavelength = wavelength
        if self._max_wavelength is None or wavelength > self._max_wavelength: self._max_wavelength = wavelength

        # If requested, the patch is returned
        patch = None

        # Check if a data point of this instrument has already been plotted
        if label not in used_labels:

            lower_flux = flux + error.lower
            upper_flux = flux + error.upper

            if lower_flux <= 0:
                flux_lower_flux = float("-inf")
            else: flux_lower_flux = flux/lower_flux

            flux_upper_flux = flux/upper_flux

            error_bar = np.array([[np.fabs(np.log10(flux_lower_flux)), np.fabs(np.log10(flux_upper_flux))]]).T
            used_labels.append(label)

            axis.errorbar(wavelength, np.log10(flux), yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)
            patch = axis.plot(wavelength, np.log10(flux), marker, markersize=7, color=color, markeredgecolor='black', markerfacecolor=color, label=label)

        # A data point of this instrument has already been plotted
        else:

            error_bar = np.array([[np.fabs(np.log10(flux)-np.log10(flux + error.lower)), np.fabs(np.log10(flux) - np.log10(flux + error.upper))]]).T
            axis.errorbar(wavelength, np.log10(flux), yerr=error_bar, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)

        # Return the patch if requested
        if return_patch: return patch

    # -----------------------------------------------------------------

    def draw_observation(self, instruments, bands, wavelengths, fluxes, errors, colors):

        """
        This function ...
        :param instruments:
        :param bands:
        :param wavelengths:
        :param fluxes:
        :param errors:
        :param colors:
        :return:
        """

        # Get labels and descriptions
        labels, descriptions = get_labels_and_descriptions(instruments, bands)

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

            # Get next color
            color = colors[k]

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            self.plot_wavelength(self._main_axis, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

    # -----------------------------------------------------------------

    def draw_model(self, axis, wavelengths, fluxes, linestyle, label, errors=None, linecolor='black', errorcolor='cyan', adjust_extrema=True):

        """
        This function ...
        :param axis:
        :param wavelengths:
        :param fluxes:
        :param linestyle:
        :param label:
        :param errors:
        :param linecolor:
        :param errorcolor:
        :param adjust_extrema:
        :return:
        """

        if adjust_extrema:

            # Keep track of the minimal and maximal flux
            min_flux_model = min(filter(None, fluxes)) # ignores zeros; see http://stackoverflow.com/questions/21084714/find-the-lowest-value-that-is-not-null-using-python
            max_flux_model = max(filter(None, fluxes)) # idem.
            if min_flux_model < self._min_flux or self._min_flux is None: self._min_flux = min_flux_model
            if max_flux_model > self._max_flux or self._max_flux is None: self._max_flux = max_flux_model

            # Keep track of the minimal and maximal wavelength
            min_wavelength_model = min(wavelengths)
            max_wavelength_model = max(wavelengths)
            if min_wavelength_model < self._min_wavelength or self._min_wavelength is None: self._min_wavelength = min_wavelength_model
            if max_wavelength_model > self._max_wavelength or self._max_wavelength is None: self._max_wavelength = max_wavelength_model

        # Plot the data points (as a line) on the axis
        log_model = np.log10(fluxes)
        axis.plot(wavelengths, log_model, linestyle, color=linecolor, label=label)

        # Plot errors
        if errors is not None:

            bottom = []
            top = []
            for j in range(len(errors)):

                value = fluxes[j]
                error = errors[j]

                bottom.append(value + error.lower)
                top.append(value + error.upper)

            bottom = np.array(bottom)
            top = np.array(top)

            log_bottom = np.log10(bottom)
            log_top = np.log10(top)

            axis.fill_between(wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor=errorcolor, edgecolor=errorcolor, interpolate=True, alpha=0.5)
            axis.plot([], [], color=errorcolor, linewidth=10, label='spread')

    # -----------------------------------------------------------------

    def finish_plot(self, path, for_legend_patches=None, for_legend_parameters=None, extra_legend=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Axis limits are now definite
        if self.min_flux is None: self.min_flux = self._min_flux
        if self.max_flux is None: self.max_flux = self._max_flux
        if self.min_wavelength is None: self.min_wavelength = self._min_wavelength
        if self.max_wavelength is None: self.max_wavelength = self._max_wavelength

        # Format residual axes
        for res_axis in self._residual_axes:

            # Set linestyle and limit for axis2
            res_axis.axhline(y=0., color='black', ls='-.')
            res_axis.set_ylim(-95, 95)

            # Set axis label
            res_axis.set_xscale('log')
            res_axis.set_ylabel(r"Res. $[\%]$", fontsize='large')

        # Set x label of the last residual plot
        if len(self._residual_axes) > 0: self._residual_axes[len(self._residual_axes)-1].set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')
        else: self._main_axis.set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')

        self._figure.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in self._figure.axes[:-1]], visible=False)

        # Add axis labels and a legend
        self._main_axis.set_xscale('log')
        self._main_axis.set_ylabel(r"Log $F_\nu$$[Jy]$", fontsize='large')

        # Enable grid
        self._main_axis.grid('on')

        # Set flux axis limits
        plot_min, plot_max = get_plot_flux_limits(self.min_flux, self.max_flux)
        self._main_axis.set_ylim((plot_min, plot_max))

        # Set wavelength axis limits
        plot_min_wavelength, plot_max_wavelength = get_plot_wavelength_limits(self.min_wavelength, self.max_wavelength)
        self._main_axis.set_xlim(plot_min_wavelength, plot_max_wavelength)

        # Add the legend
        if for_legend_patches is not None:

            # Set legend
            # fancybox=True makes the legend corners rounded
            self._main_axis.legend([l[0] for l in for_legend_patches], for_legend_parameters, numpoints=1, loc="lower right", frameon=True, ncol=2, fontsize=11, shadow=False)

            # Extra legend
            if extra_legend is not None: self._main_axis.add_artist(extra_legend)

        # No extra legend, no patches for the legend
        else: self._main_axis.legend(numpoints=1, loc="lower right", frameon=True, ncol=2, fontsize=11, shadow=False)

        # Add title if requested
        if self.title is not None: self._figure.suptitle(self.title, fontsize=14, fontweight='bold')

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

def find_reference_flux(instrument, band, wavelength, reference_instruments, reference_bands, reference_wavelengths, reference_fluxes):

    """
    This function ...
    :param instrument:
    :param band:
    :param wavelength:
    :param reference:
    :return:
    """

    # Loop over the entries in the reference SED
    for i in range(len(reference_wavelengths)):

        if (reference_instruments[i] == instrument and reference_bands == band) or np.isclose(reference_wavelengths[i], wavelength):
            return reference_fluxes[i]

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

def get_unique_labels(observations):

    """
    This function ...
    :param observations:
    :return:
    """

    # Unique labels
    unique_labels = set()

    # Get all unique labels
    for label in observations:

        # Get the observed SED
        observation = observations[label]

        # Get wavelengths, fluxes, instruments, bands, errors
        instruments = observation.instruments()
        bands = observation.bands()

        # Get labels and descriptions
        labels, descriptions = get_labels_and_descriptions(instruments, bands)

        # Add labels
        for _label in labels: unique_labels.add(_label)

    # Convert unique labels (set) to list
    unique_labels = list(unique_labels)

    # Return the list of unique labels
    return unique_labels

# -----------------------------------------------------------------

def get_plot_flux_limits(min_flux, max_flux):

    """
    This function ...
    :return:
    """

    log_min_flux = np.log10(min_flux)
    log_max_flux = np.log10(max_flux)
    plot_min = 0.9 * log_min_flux if log_min_flux > 0 else 1.1 * log_min_flux
    plot_max = 1.1 * log_max_flux

    return plot_min, plot_max

# -----------------------------------------------------------------

def get_plot_wavelength_limits(min_wavelength, max_wavelength):

    """
    This function ...
    :param min_wavelength:
    :param max_wavelength:
    :return:
    """

    log_min_wavelength = np.log10(min_wavelength)
    log_max_wavelength = np.log10(max_wavelength)

    plot_min_log_wavelength = math.floor(log_min_wavelength)
    plot_max_log_wavelength = math.ceil(log_max_wavelength)

    plot_min_wavelength = 10.**plot_min_log_wavelength
    plot_max_wavelength = 10.**plot_max_log_wavelength

    return plot_min_wavelength, plot_max_wavelength

# -----------------------------------------------------------------
