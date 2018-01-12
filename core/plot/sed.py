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
from matplotlib import colors
import matplotlib.patches as patches
from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
from scipy.interpolate import interp1d
from operator import itemgetter

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..data.sed import ObservedSED, SED, load_sed, load_multiple_seds
from ..basics.plot import MPLFigure, BokehFigure, BokehPlot, mpl, bokeh, dark_pretty_colors, pretty_colors, filled_markers
from ..tools import types
from ..filter.broad import BroadBandFilter
from ..basics.errorbar import ErrorBar
from ..tools import numbers
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

rc('text', usetex=True)

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']

distinguishable_colormaps = ["spring", "winter", "copper", "cool", "PRGn", "coolwarm"]
other_colormaps = ["YlGn", "YlGnBu", "YlOrRd", "Purples", "inferno", "plasma", "Spectral", "terrain"]

# http://matplotlib.org/examples/color/named_colors.html
color_hex = colors.cnames

# Add the single letter colors
for name, rgb in colors.ColorConverter.colors.items():
    hex_ = colors.rgb2hex(rgb)
    color_hex[name] = hex_

color_plt_identifiers = color_hex.keys()

# -----------------------------------------------------------------

def plot_sed(sed, label=None, path=None, title=None):

    """
    This function ...
    :param sed:
    :param label:
    :param path:
    :param title:
    :return:
    """

    # Create a new SEDPlotter instance
    plotter = SEDPlotter()

    # Determine label
    if label is None:
        if isinstance(sed, ObservedSED): label = "observation"
        elif isinstance(sed, SED): label = "model"
        else: raise ValueError("Invalid SED object")

    # Add the SED
    plotter.add_sed(sed, label)

    # Make the plot
    plotter.run(title=title, output=path)

# -----------------------------------------------------------------

class SEDPlotter(Configurable):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param *args
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SEDPlotter, self).__init__(*args, **kwargs)

        # The plot title
        self.title = None

        # Create ordered dictionaries for the model and observed SEDs (the order of adding SEDs is remembered)
        self.models = []
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

        # The output path
        self.out_path = None

        # Properties
        self.transparent = False

        # The figure
        self.figure = None

        # The main plot and the residual plots
        self.main_plot = None
        self.residual_plots = []

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        This function ...
        :return:
        """

        return len(self.models) + len(self.observations)

    # -----------------------------------------------------------------

    @property
    def no_seds(self):

        """
        This function ...
        :return:
        """

        return self.nseds == 0

    # -----------------------------------------------------------------

    @property
    def nmodels(self):

        """
        This function ...
        :return:
        """

        return len(self.models)

    # -----------------------------------------------------------------

    @property
    def has_models(self):

        """
        This function ...
        :return:
        """

        return self.nmodels > 0

    # -----------------------------------------------------------------

    @property
    def nobservations(self):

        """
        This function ...
        :return:
        """

        return len(self.observations)

    # -----------------------------------------------------------------

    @property
    def has_observations(self):

        """
        This function ...
        :return:
        """

        return self.nobservations > 0

    # -----------------------------------------------------------------

    def add_sed(self, sed, label, residuals=True, ghost=False):

        """
        This function ...
        :param sed:
        :param label:
        :param residuals: whether plotting the residual curve makes sense for this SED
        (it does not when the SED is only for one component, or a contribution to the total SED)
        :param ghost:
        :return:
        """

        # Add observed or modeled SED
        if isinstance(sed, ObservedSED):

            # Make a copy with only broad band filters
            only_broad = sed.only_broad_band()

            # Check whether there are any points in the SED
            if len(only_broad) == 0: log.warning("The '" + label + "' SED contains no (broad band) photometry points")

            # Set the SED
            self.observations[label] = only_broad

        elif isinstance(sed, SED): self.models.append((label, sed, residuals, ghost))
        else: raise ValueError("The SED must be an SED or ObservedSED instance")

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDPlotter, self).setup(**kwargs)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the axis limits
        self.min_wavelength = kwargs.pop("min_wavelength", None)
        self.max_wavelength = kwargs.pop("max_wavelength", None)
        self.min_flux = kwargs.pop("min_flux", None)
        self.max_flux = kwargs.pop("max_flux", None)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        if self.out_path is None and "output" in self.config and self.config.output is not None:
            full_output_path = fs.absolute_or_in(self.config.output, self.config.path)
            if fs.has_extension(full_output_path):
                directory_path = fs.directory_of(full_output_path)
                if not fs.is_directory(directory_path): fs.create_directory(directory_path)
            elif not fs.is_directory(full_output_path): fs.create_directory(full_output_path)
            self.out_path = full_output_path

        # Add SED files present in the current working directory (if nothing is added manually)
        if self.no_seds: self.load_seds()

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.config.plot.xsize

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        if self.nobservations > 1 and self.nmodels > 0:
            main_size = self.config.plot.ysize * self.config.plot.main_relsize
            other_size = float(self.config.plot.ysize) - main_size
            return main_size + (self.nobservations - 1) * other_size
        else: return self.config.plot.ysize

    # -----------------------------------------------------------------

    @property
    def figsize(self):

        """
        This function ...
        :return:
        """

        return (self.xsize, self.ysize)

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading SEDs ...")

        # From file paths
        if self.config.seds is not None: self.load_sed_files()

        # From files in working directory
        else: self.load_seds_cwd()

    # -----------------------------------------------------------------

    def load_sed_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading SED files ...")

        # Loop over the SED files
        for path in self.config.seds:

            # Get the SED name
            name = fs.strip_extension(fs.name(path))

            # Try getting multiple SEDs
            if self.config.multi:

                seds = load_multiple_seds(path, as_dict=True, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)
                for label in seds: self.add_sed(seds[label], label=label)

            # One SED per file
            else:

                # Load SED
                sed = load_sed(path, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)

                # Add the SED
                self.add_sed(sed, label=name)

    # -----------------------------------------------------------------

    def load_seds_cwd(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading SED files in the current working directory ...")

        # Loop over the files
        for path, name in fs.files_in_path(self.config.path, extension="dat", returns=["path", "name"]):

            # Skip emission lines
            if "Lines" in name: continue

            # Try getting multiple SEDs
            if self.config.multi:

                seds = load_multiple_seds(path, as_dict=True, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)
                for label in seds: self.add_sed(seds[label], label=label)

            # One SED per file
            else:

                # Load SED
                sed = load_sed(path, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)

                # Add the SED
                self.add_sed(sed, label=name)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the SED plotter ...")

        # Set default values for all attributes
        #self.models = OrderedDict()
        #self.models_no_residuals = OrderedDict()
        self.models = []
        self.observations = OrderedDict()
        self._min_wavelength = None
        self._max_wavelength = None
        self._min_flux = None
        self._max_flux = None
        self.min_wavelength = None
        self.max_wavelength = None
        self.min_flux = None
        self.max_flux = None
        self.main_plot = None
        self.residual_plots = []

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the SED plot ...")

        # No models
        if not self.has_models: self.plot_no_models()

        # With models
        else: self.plot_with_models()

    # -----------------------------------------------------------------

    def plot_no_models(self):

        """
        This function ...
        :return:
        """

        # One observation or more observations
        if self.nobservations == 1: self.plot_one_observation()
        else: self.plot_more_observations()

    # -----------------------------------------------------------------

    def plot_one_observation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting one observed SED ...")

        # Setup the figure
        self.main_plot = self.figure.create_one_plot()

        # Determine color map class
        colormap = plt.get_cmap("rainbow")

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
        fluxes = observation.photometry(unit=self.config.unit, add_unit=False)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit=self.config.unit, add_unit=False)

        # Create colors
        colors = colormap(np.linspace(0, 1, len(wavelengths)))

        # Plot
        self.draw_observation(instruments, bands, wavelengths, fluxes, errors, colors)

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    @property
    def has_any_observation_errors(self):

        """
        This function ...
        :return:
        """

        from ..tools import sequences

        # Loop over the observations
        for label in self.observations:

            # Get the errors
            observation = self.observations[label]
            errors = observation.errors(unit=self.config.unit, add_unit=False)

            # Has errors?
            has_errors = sequences.has_not_none(errors)

            # Certainly one observation with errors
            if has_errors: return True

        # No observation with errors
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_labels(self):

        """
        This function ...
        :return:
        """

        return get_unique_labels(self.observations)

    # -----------------------------------------------------------------

    @lazyproperty
    def most_npoints_observation_label(self):

        """
        This function ...
        :return:
        """

        most_npoints = None
        most_npoints_label = None

        for label in self.observations:
            observation = self.observations[label]
            if most_npoints_label is None or len(observation) > most_npoints:
                most_npoints = len(observation)
                most_npoints_label = label

        return most_npoints_label

    # -----------------------------------------------------------------

    @property
    def most_npoints_observation(self):

        """
        This function ...
        :return:
        """

        return self.observations[self.most_npoints_observation_label]

    # -----------------------------------------------------------------

    def plot_more_observations(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting multiple observed SEDs ...")

        nplots = 2
        height_ratios = [4, 1]

        # CREATE MAIN PLOT AND RESIDUAL PLOT(S)
        plots = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios)
        self.main_plot = plots[0]
        residual_plot = plots[1]
        self.residual_plots.append(residual_plot)

        # Keep track of whether we are plotting the first observed SED
        #first = True
        #reference_sed = None

        # Make iterable from distinct colors
        different_colors = iter(dark_pretty_colors)

        # Unique labels
        #unique_labels = get_unique_labels(self.observations)

        # Markers for the unique labels
        markers = filled_markers[:len(self.unique_labels)]

        # Used labels
        used_labels = []

        legend_patches = []
        legend_labels = []

        legend_rectangles = []
        rectangle_labels = []

        # Remember colors for each observation
        observation_colors = dict()

        # Determine the reference SED
        reference_sed_label = self.most_npoints_observation_label
        #print("reference", reference_sed_label)
        reference_sed = self.observations[reference_sed_label]

        # Loop over the different observed SEDs
        for label in self.observations:

            # Debugging
            log.debug("Plotting the data points of the '" + label + "' observed SED ...")

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            fluxes = observation.photometry(unit=self.config.unit, add_unit=False)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit=self.config.unit, add_unit=False)

            # Create color range
            observation_color = color_hex[next(different_colors)]
            #colors = [] * len(wavelengths)
            observation_colors[label] = observation_color

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                # Check filter
                instrument = instruments[k]
                band = bands[k]
                if self.config.ignore_filters is not None:
                    fltr = BroadBandFilter.from_instrument_and_band(instrument, band)
                    if fltr in self.config.ignore_filters: continue

                # Check validity of flux value
                if fluxes[k] <= 0.0:
                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                # Get marker
                marker = markers[self.unique_labels.index(labels[k])]

                # Set zero error if other observations also have errors (BECAUSE CALLING MATPLOTLIB'S ERRORBAR AND AFTER THAT PLOT DOESN'T WORK APPARENTLY!!)
                error = errors[k]
                #if sequences.has_any(with_errors) and error is None:
                if self.has_any_observation_errors and error is None: error = ErrorBar.zero()

                # Ignore upper limits?
                if self.config.ignore_upper and error.only_upper: continue

                #print(used_labels)

                # Plot the flux data point on the main axis with the specified marker and color
                patch, new_label = self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], error, marker, observation_color, return_new=True)

                # Add the patch and label
                if patch is not None and new_label:
                    legend_patches.append(patch)
                    legend_labels.append(labels[k])

                # Plot on axis 2
                if label == reference_sed_label:

                    if errors[k] is not None:

                        value = 0.0
                        error = errors[k] / fluxes[k] * 100.

                    else: value = error = None

                else:

                    reference_flux = find_reference_flux(instruments[k], bands[k], wavelengths[k], reference_sed.instruments(), reference_sed.bands(), reference_sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False), reference_sed.photometry(unit=self.config.unit, add_unit=False))

                    #print(fluxes[k], reference_flux)

                    if reference_flux is None:

                        value = None
                        error = None

                    else:

                        if reference_flux == 0: value = error = None
                        else:

                            if errors[k] is not None:

                                value = (fluxes[k] - reference_flux) / reference_flux * 100.
                                error = errors[k] / reference_flux * 100.0

                            #else: value = error = None
                            elif fluxes[k] is not None:

                                value = (fluxes[k] - reference_flux) / reference_flux * 100.
                                error = ErrorBar.zero()

                            else: value = error = None

                if value is not None and error is not None:

                    yerr, lolims, uplims = process_errorbar(error)
                    residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=observation_color, markeredgecolor='black', ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims)

            # The next observation is not the first anymore
            first = False

            # If color gradient is used for observations instead of single colors ...
            #ax3 = next(colormap_axes)
            #ax3.imshow(gradient, aspect='auto', cmap=colormap)

            # Create rectangle for this observation
            rectangle = patches.Rectangle((0, 0), 1, 1, fc=observation_color)
            legend_rectangles.append(rectangle)
            rectangle_labels.append(label)

        # Extra legend: the different observations
        # fancybox=True makes the legend corners rounded
        observations_legend = self.main_plot.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Finish the plot
        self.finish_plot(for_legend_patches=legend_patches, for_legend_parameters=legend_labels, extra_legend=observations_legend)

    # -----------------------------------------------------------------

    def plot_with_models(self):

        """
        This function ...
        :return:
        """

        # Only models
        if self.nobservations == 0: self.plot_only_models()

        # One observation
        elif self.nobservations == 1: self.plot_one_observation_with_models()

        # Multiple observations
        else: self.plot_more_observations_with_models()

    # -----------------------------------------------------------------

    def plot_only_models(self):

        """
        This function ...
        :param:
        :return:
        """

        # Debugging
        log.debug("Plotting only model SEDs ...")

        # Setup the figure
        #self.main_plot = self.figure.ax

        #plots = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios)
        #self.main_plot = plots[0]
        #self.residual_plots = plots[1:]

        # Create the main plot
        self.main_plot = self.figure.create_one_plot()

        # Add model SEDs

        counter_residuals = 0
        counter_no_residals = 0

        line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        line_styles_models_no_residuals = ["-"] * len(self.models)

        # Loop over the model SEDs
        for model_label, sed, plot_residuals, ghost in self.models:

            if ghost:

                # Get fluxes, wavelengths and errors
                fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                # Plot the model SED as a grey line (no errors)
                self.draw_model(self.main_plot, wavelengths, fluxes, "-", linecolor="lightgrey")

            elif plot_residuals:

                #print(self.config.unit.density)
                #print(sed.unit)

                # Get fluxes, wavelengths and errors
                fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                errors = sed.errors(unit=self.config.unit, add_unit=False) if sed.has_errors else None

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles[counter_residuals], model_label, errors=errors)

                counter_residuals += 1

            else:

                # Get fluxes, wavelengths and errors
                fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                errors = sed.errors(unit=self.config.unit, add_unit=False) if sed.has_errors else None

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles_models_no_residuals[counter_no_residals], model_label, errors=errors, linecolor=line_colors_models_no_residuals[counter_no_residals], adjust_extrema=False)

                counter_no_residals += 1

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    @property
    def nmodels_for_residuals(self):

        """
        Thisf unction ...
        :return:
        """

        count = 0

        # Loop over the models
        for model_label, sed, plot_residuals, ghost in self.models:
            if plot_residuals: count += 1
        return count

    # -----------------------------------------------------------------

    def plot_one_observation_with_models(self):

        """
        This function ...
        :param:
        :return:
        """

        # Debugging
        log.debug("Plotting one observed SED with model SEDs ...")

        # Setup the figure
        #gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        #self.main_plot = plt.subplot(gs[0])
        #ax2 = plt.subplot(gs[1], sharex=self.main_plot)
        #self.residual_plots.append(ax2)

        if self.config.residual_reference == "observations": nplots = 2
        elif self.config.residual_reference == "models": nplots = 1 + self.nmodels_for_residuals
        else: raise ValueError("")

        height_ratios = [4] + [1] * (nplots-1)

        # CREATE MAIN PLOT AND RESIDUAL PLOT(S)
        #self.main_plot, residual_plot = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios)
        plots = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios)
        #self.residual_plots.append(residual_plot)
        self.main_plot = plots[0]
        self.residual_plots = plots[1:]

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
        fluxes = observation.photometry(unit=self.config.unit, add_unit=False)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit=self.config.unit, add_unit=False)

        # Get labels and descriptions
        labels, descriptions = get_labels_and_descriptions(instruments, bands)

        # Determine color map class
        colormap = plt.get_cmap("rainbow")
        #print(colormap)

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

            # Check filter
            instrument = instruments[k]
            band = bands[k]
            if self.config.ignore_filters is not None:
                fltr = BroadBandFilter.from_instrument_and_band(instrument, band)
                #print(fltr)
                #print(self.config.ignore_filters)
                if fltr in self.config.ignore_filters: continue

            # Check validity of flux value
            if fluxes[k] <= 0.0:
                log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                continue

            #print(color_range)

            # Get next color
            color = next(color_range)

            #print(color)

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            # axis, label, used_labels, wavelength, flux, error, marker, color, return_patch=False
            #print(labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)
            # GALEX [] 0.153507951158 0.0686478012556 None o [ 0.5  0.   1.   1. ]
            self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

            # Observations as reference: plot at 0.0 (all in one panel)
            if self.config.residual_reference == "observations":

                # Get the (only) residual plot
                residual_plot = self.residual_plots[0]

                # Plot point at y=0.0 with errorbar on axis 2
                value = 0.0
                if errors[k] is not None:
                    error = errors[k] / fluxes[k] * 100.
                    yerr, lolims, uplims = process_errorbar(error)
                    residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims)

            # Models as reference: plot points at the relative difference with the models (one panel for each model)
            elif self.config.residual_reference == "models":

                from astropy.units import Unit

                residual_plot_index = 0

                # Loop over the models
                for model_label, sed, plot_residuals, ghost in self.models:

                    # Don't plot residuals for this model
                    if not plot_residuals: continue

                    #model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                    #sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                    actual_wavelength = wavelengths[k] * Unit(self.config.wavelength_unit)

                    # Wavelength out of bounds of model SED
                    if actual_wavelength < sed.min_wavelength or actual_wavelength > sed.max_wavelength:
                        residual_plot_index += 1
                        continue

                    model_value = sed.photometry_at(actual_wavelength, self.config.unit, add_unit=False)
                    rel_residual = (fluxes[k] - model_value) / model_value * 100.  # NORMALIZE TO MODEL VALUE

                    #print("MODEL VALUE:", model_value)
                    #print("FLUX:", fluxes[k])
                    #print("REL RESIDUAL:", rel_residual)

                    residual_plot = self.residual_plots[residual_plot_index]
                    #value = 0.0
                    if errors[k] is not None:
                        error = errors[k] / model_value * 100. # NORMALIZE TO MODEL VALUE
                        yerr, lolims, uplims = process_errorbar(error)
                        residual_plot.errorbar(wavelengths[k], rel_residual, yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims)

                    # Next residual plot
                    residual_plot_index += 1

        line_styles_models = line_styles
        line_colors_models = ['black'] * len(self.models)

        line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        line_styles_models_no_residuals = ["-"] * len(self.models)

        counter = 0

        # Plot models on the residual axis, WHEN OBSERVATION IS THE REFERENCE
        for model_label, sed, plot_residuals, ghost in self.models:

            if not plot_residuals: continue

            if ghost:

                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                min_sed_wavelength = min(sed_wavelengths)
                max_sed_wavelength = max(sed_wavelengths)
                wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                residual_plot.plot(wavelengths_residuals, residuals, linestyle="-", color="lightgrey")

            elif self.config.residual_reference == "observations":

                #log_model = np.log10(sed.fluxes(unit=self.config.unit, add_unit=False))
                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                min_sed_wavelength = min(sed_wavelengths)
                max_sed_wavelength = max(sed_wavelengths)
                wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                residual_plot.plot(wavelengths_residuals, residuals, linestyle=line_styles_models[counter], color=line_colors_models[counter], label='model')

                counter += 1

            # elif self.config.residual_reference == "models":
            #
            #     model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
            #     sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            #     f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')
            #
            #     #min_sed_wavelength = min(sed_wavelengths)
            #     #max_sed_wavelength = max(sed_wavelengths)

        counter = 0
        counter_no_residuals = 0

        # Add model SEDs
        for model_label, sed, plot_residuals, ghost in self.models:

            if self.nmodels == 1: model_label = "model"

            if ghost:

                log_model = np.log10(sed.photometry(unit=self.config.unit, add_unit=False))
                self.main_plot.plot(sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False), log_model, ls="-", color="lightgrey")

            elif plot_residuals:

                fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                log_model = np.log10(fluxes)
                self.main_plot.plot(wavelengths, log_model, linestyle=line_styles_models[counter], color=line_colors_models[counter], label=model_label)

                # Get lowest and highest flux and wavelength for this curve
                lowest = np.min(fluxes)
                highest = np.max(fluxes)
                #lowest_lambda = np.min(wavelengths)
                #highest_lambda = np.max(wavelengths)

                # Adapt flux range
                if self._min_flux is None or lowest < self._min_flux: self._min_flux = lowest
                if self._max_flux is None or highest > self._max_flux: self._max_flux = highest

                # Keep track of the minimal and maximal wavelength
                #if self._min_wavelength is None or lowest_lambda < self._min_wavelength: self._min_wavelength = lowest_lambda
                #if self._max_wavelength is None or highest_lambda > self._max_wavelength: self._max_wavelength = highest_lambda

                if sed.has_errors:

                    bottom = []
                    top = []
                    for j in range(len(sed.errors)):

                        value = sed.photometry[j]
                        bottom.append(value + sed.errors[j][0])
                        top.append(value + sed.errors[j][1])

                    bottom = np.array(bottom)
                    top = np.array(top)

                    lowest = np.min(bottom)
                    highest = np.max(top)

                    # Adapt flux range
                    if self._min_flux is None or lowest < self._min_flux: self._min_flux = lowest
                    if self._max_flux is None or highest > self._max_flux: self._max_flux = highest

                    log_bottom = np.log10(bottom)
                    log_top = np.log10(top)

                    self.main_plot.fill_between(sed.wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor='cyan', edgecolor='cyan', interpolate=True, alpha=0.5)
                    self.main_plot.plot([], [], color='cyan', linewidth=10, label='spread')

                counter += 1

            else:

                fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                log_model = np.log10(fluxes)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                self.main_plot.plot(wavelengths, log_model, linestyle=line_styles_models_no_residuals[counter_no_residuals], color=line_colors_models_no_residuals[counter_no_residuals], label=model_label)

                ## SAME AS ABOVE::
                # Get lowest and highest flux and wavelength for this curve
                lowest = np.min(fluxes)
                highest = np.max(fluxes)
                #lowest_lambda = np.min(wavelengths)
                #highest_lambda = np.max(wavelengths)

                # Adapt flux range
                if self._min_flux is None or lowest < self._min_flux: self._min_flux = lowest
                if self._max_flux is None or highest > self._max_flux: self._max_flux = highest

                # Keep track of the minimal and maximal wavelength
                #if self._min_wavelength is None or lowest_lambda < self._min_wavelength: self._min_wavelength = lowest_lambda
                #if self._max_wavelength is None or highest_lambda > self._max_wavelength: self._max_wavelength = highest_lambda
                ##

                if sed.has_errors:

                    bottom = []
                    top = []
                    for j in range(len(sed.errors)):
                        value = sed.photometry[j]
                        bottom.append(value + sed.errors[j][0])
                        top.append(value + sed.errors[j][1])

                    bottom = np.array(bottom)
                    top = np.array(top)

                    lowest = np.min(bottom)
                    highest = np.max(top)

                    # Adapt flux range
                    if self._min_flux is None or lowest < self._min_flux: self._min_flux = lowest
                    if self._max_flux is None or highest > self._max_flux: self._max_flux = highest

                    log_bottom = np.log10(bottom)
                    log_top = np.log10(top)

                    self.main_plot.fill_between(sed.wavelengths, log_bottom, log_top,
                                                 where=log_top <= log_bottom, facecolor='cyan', edgecolor='cyan',
                                                 interpolate=True, alpha=0.5)
                    self.main_plot.plot([], [], color='cyan', linewidth=10, label='spread')

                counter_no_residuals += 1

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    def plot_more_observations_with_models(self):

        """
        This function ...
        :return:
        """

        # Debuggging
        log.debug("Plotting multiple observed SEDs with models ...")

        # Determine number of plot rows
        if self.config.residual_reference == "observations": nplots = 1 + self.nobservations
        elif self.config.residual_reference == "models": nplots = 1 + self.nmodels_for_residuals
        else: raise ValueError("Invalid residual reference")

        #height_ratios = [4]
        #for _ in range(number_of_observations): height_ratios.append(1)

        #nplots = number_of_observations + 1
        #height_ratios = [4] + [1] * (nplots - 1)
        #height_ratios = [4] + [1] * number_of_observations
        height_ratios = [4] + [1] * (nplots - 1)

        # Setup the figure
        #gs = gridspec.GridSpec(1 + number_of_observations, 1, height_ratios=height_ratios)
        #self.main_plot = plt.subplot(gs[0])
        #for i in range(number_of_observations):
        #    ax = plt.subplot(gs[1+i], sharex=self.main_plot)
        #    self.residual_plots.append(ax)
        #residual_axes = iter(self.residual_plots)

        # CREATE MAIN PLOT AND RESIDUAL PLOT(S)
        plots = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios)
        self.main_plot = plots[0]
        self.residual_plots = plots[1:]

        # Make iterable from color map names
        #color_maps = iter(distinguishable_colormaps+other_colormaps)

        # Make iterable from distinct colors
        colors = iter(pretty_colors)

        # Get unique labels
        unique_labels = get_unique_labels(self.observations)

        # Markers for the unique labels
        markers = filled_markers[:len(unique_labels)]

        # Used labels
        used_labels = []

        line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        line_styles_models_no_residuals = ["-"] * len(self.models)

        legend_patches = []
        legend_labels = []

        legend_rectangles = []
        rectangle_labels = []

        # Loop over the different observed SEDs
        #index = 0
        observation_index = 0
        for label in self.observations:

            # Get the axis of the next residual plot
            #ax2 = next(self.residual_plots)
            #ax2 = self.residual_plots[index]
            #index += 1

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            fluxes = observation.photometry(unit=self.config.unit, add_unit=False)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit=self.config.unit, add_unit=False)

            # Get labels and descriptions
            labels, descriptions = get_labels_and_descriptions(instruments, bands)

            #print(label)
            #print(fluxes)

            # Create color range
            # if number_of_observations <= 3:
            #     # Determine color map class
            #     colormap = plt.get_cmap(next(color_maps))
            #     color_range = iter(colormap(np.linspace(0, 1, len(wavelengths))))
            # else: color_range = iter([color_hex[next(colors)]] * len(wavelengths))
            color_range = iter([color_hex[next(colors)]] * len(wavelengths))

            # Loop over the wavelengths
            for k in range(len(wavelengths)):

                #print(k)

                # Check filter
                instrument = instruments[k]
                band = bands[k]
                if self.config.ignore_filters is not None:
                    fltr = BroadBandFilter.from_instrument_and_band(instrument, band)
                    #print(fltr)
                    #print(self.config.ignore_filters)
                    if fltr in self.config.ignore_filters: continue

                # Check validity of flux value
                if fluxes[k] <= 0.0:
                    log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                    continue

                #wavelength = wavelengths[k]
                #flux = fluxes[k]
                #print(wavelength, flux)

                # Get next color
                color = next(color_range)

                # Get marker
                marker = markers[unique_labels.index(labels[k])]

                # Set zero error if other observations also have errors (BECAUSE CALLING MATPLOTLIB'S ERRORBAR AND AFTER THAT PLOT DOESN'T WORK APPARENTLY!!)
                error = errors[k]
                if self.has_any_observation_errors and error is None: error = ErrorBar.zero()

                #print(labels[k])

                # Plot the flux data point on the main axis with the specified marker and color
                patch = self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], error, marker, color)

                # Add patch and label
                if patch is not None:
                    legend_patches.append(patch)
                    legend_labels.append(labels[k])

                # Observations as reference: plot at 0.0 (all in one panel)
                if self.config.residual_reference == "observations":

                    # Get the residual plot
                    residual_plot = self.residual_plots[observation_index]

                    # Plot measurement points on residual plot
                    value = 0.0
                    if errors[k] is not None:
                        error = errors[k] / fluxes[k] * 100.
                        yerr, lolims, uplims = process_errorbar(error)
                        residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims)

                # Models as reference: plot points at the relative difference with the models (one panel for each model)
                elif self.config.residual_reference == "models":

                    from astropy.units import Unit

                    residual_plot_index = 0

                    # Loop over the models
                    for model_label, sed, plot_residuals, ghost in self.models:

                        # Don't plot residuals for this model
                        if not plot_residuals: continue

                        # model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                        # sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                        actual_wavelength = wavelengths[k] * Unit(self.config.wavelength_unit)

                        # Wavelength out of bounds of model SED
                        if actual_wavelength < sed.min_wavelength or actual_wavelength > sed.max_wavelength:
                            residual_plot_index += 1
                            continue

                        model_value = sed.photometry_at(actual_wavelength, self.config.unit, add_unit=False)
                        rel_residual = (fluxes[k] - model_value) / model_value * 100.  # NORMALIZE TO MODEL VALUE

                        # print("MODEL VALUE:", model_value)
                        # print("FLUX:", fluxes[k])
                        # print("REL RESIDUAL:", rel_residual)

                        # Get the residual plot for this model
                        residual_plot = self.residual_plots[residual_plot_index]

                        # value = 0.0
                        if errors[k] is None: errorbar = ErrorBar.zero()
                        else: errorbar = errors[k]

                        error = errorbar / model_value * 100.  # NORMALIZE TO MODEL VALUE

                        yerr, lolims, uplims = process_errorbar(error)

                        # Plot on residual axes
                        residual_plot.errorbar(wavelengths[k], rel_residual, yerr=yerr, fmt=marker,
                                               markersize=7, color=color, markeredgecolor='black', ecolor=color,
                                               capthick=2, lolims=lolims, uplims=uplims)

                        # Next residual plot
                        residual_plot_index += 1

                # Invalid
                else: raise RuntimeError("Invalid value for 'residual_reference'")

            # IF OBSERVATIONS ARE THE REFERENCE FOR THE RESIDUAL AXES
            if self.config.residual_reference == "observations":

                # Get the residual plot
                residual_plot = self.residual_plots[observation_index]

                # Residuals
                counter = 0
                for model_label, sed, plot_residuals, ghost in self.models:

                    if not plot_residuals: continue

                    model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                    sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                    f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                    min_sed_wavelength = min(sed_wavelengths)
                    max_sed_wavelength = max(sed_wavelengths)
                    wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                    wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                    fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                    residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                    if ghost: residual_plot.plot(wavelengths_residuals, residuals, ls="-", color='lightgrey')
                    else:
                        residual_plot.plot(wavelengths_residuals, residuals, ls=line_styles[counter], color='black', label=model_label)
                        counter += 1

            # Create rectangle for this observation
            rectangle = patches.Rectangle((0, 0), 1, 1, fc=color)
            legend_rectangles.append(rectangle)
            rectangle_labels.append(label)

            # Increment the observation index
            observation_index += 1

        ## END OF LOOP OVER OBSERVATIONS

        # Add model SEDs
        counter = 0
        counter_no_residuals = 0
        for model_label, sed, plot_residuals, ghost in self.models:

            # Get fluxes, wavelengths and errors
            fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
            wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            errors = sed.errors(unit=self.config.unit, add_unit=False) if sed.has_errors else None

            if ghost:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, "-", linecolor="lightgrey", adjust_extrema=False)

            elif plot_residuals:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles[counter], linecolor="black", label=model_label, adjust_extrema=False)
                counter += 1

            else:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles_models_no_residuals[counter_no_residuals], linecolor=line_colors_models_no_residuals[counter_no_residuals], label=model_label, adjust_extrema=False)
                counter_no_residuals += 1

        # Extra legend: the different observations
        # fancybox=True makes the legend corners rounded
        observations_legend = self.main_plot.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Finish the plot
        #self.finish_plot()

        # Finish the plot
        self.finish_plot(for_legend_patches=legend_patches, for_legend_parameters=legend_labels, extra_legend=observations_legend)

    # -----------------------------------------------------------------

    def plot_wavelength(self, axis, label, used_labels, wavelength, flux, error, marker, color, return_new=False):

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
        :param return_new:
        :return:
        """

        #print(label, wavelength, flux)

        # Keep track of minimum and maximum flux
        if error is not None:

            if error.lower == numbers.min_inf:
                if self._min_flux is None or flux < self._min_flux: self._min_flux = flux
            else:
                lower_flux = flux + error.lower
                if self._min_flux is None or lower_flux < self._min_flux: self._min_flux = lower_flux

            if error.upper == numbers.inf:
                if self._max_flux is None or flux > self._max_flux: self._max_flux = flux
            else:
                upper_flux = flux + error.upper
                if self._max_flux is None or upper_flux > self._max_flux: self._max_flux = upper_flux

        else:

            if self._min_flux is None or flux < self._min_flux: self._min_flux = flux
            if self._max_flux is None or flux > self._max_flux: self._max_flux = flux

        # Keep track of the minimal and maximal wavelength
        if self._min_wavelength is None or wavelength < self._min_wavelength: self._min_wavelength = wavelength
        if self._max_wavelength is None or wavelength > self._max_wavelength: self._max_wavelength = wavelength

        # Check if a data point of this instrument has already been plotted
        new_label = label not in used_labels
        if new_label:

            used_labels.append(label)

            if error is not None:

                lower_flux = flux + error.lower
                upper_flux = flux + error.upper

                lolims = False
                uplims = False

                if lower_flux <= 0:
                    flux_lower_flux = float("-inf")
                    uplims = [[True]]
                else: flux_lower_flux = flux / lower_flux

                if error.lower == numbers.min_inf: uplims = [[True]]
                if error.upper == numbers.inf: lolims = [[True]]

                flux_upper_flux = flux / upper_flux

                yerr = np.array([[np.fabs(np.log10(flux_lower_flux)), np.fabs(np.log10(flux_upper_flux))]]).T
                patch = axis.errorbar(wavelength, np.log10(flux), yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims)

            else: patch = axis.plot(wavelength, np.log10(flux), marker=marker, markersize=7, color=color, markeredgecolor='black') #markerfacecolor=color) #label=label)

        # A data point of this instrument has already been plotted
        else:

            # There is an error bar
            if error is not None:

                lolims = False
                uplims = False

                yerr = np.array([[np.fabs(np.log10(flux)-np.log10(flux + error.lower)), np.fabs(np.log10(flux) - np.log10(flux + error.upper))]]).T
                #yerr, lolims, uplims = process_errorbar(ErrorBar(np.log10(flux)-np.log10(flux + error.lower), np.log10(flux) - np.log10(flux + error.upper)))

                if error.lower == numbers.min_inf: uplims = [[True]]
                if error.upper == numbers.inf: lolims = [[True]]

                #print("c", wavelength, flux)
                patch = axis.errorbar(wavelength, np.log10(flux), yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims)

            #else: axis.plot(wavelength, np.log10(flux), fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)
            else: patch = axis.plot(wavelength, np.log10(flux), markersize=7, color=color, markeredgecolor='black')

        # Return the patch if requested
        if return_new: return patch, new_label
        else: return patch

    # -----------------------------------------------------------------

    # instruments, bands, wavelengths, fluxes, errors, colors
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

            # Check filter
            instrument = instruments[k]
            band = bands[k]
            if self.config.ignore_filters is not None:
                fltr = BroadBandFilter.from_instrument_and_band(instrument, band)
                if fltr in self.config.ignore_filters: continue

            # Check validity of flux value
            if fluxes[k] <= 0.0:
                log.warning("Negative flux encountered for " + str(descriptions[k]) + " band")
                continue

            # Get next color
            color = colors[k]

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            #print("fluxes", fluxes)
            #print("errors", errors)
            # axis, label, used_labels, wavelength, flux, error, marker, color, return_patch=False
            self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

    # -----------------------------------------------------------------

    def draw_model(self, axis, wavelengths, fluxes, linestyle, label=None, errors=None, linecolor='black', errorcolor='cyan', adjust_extrema=True):

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

            if errors is not None:
                zipped = zip(fluxes, errors)
                lower_flux_model = min(filter(None, [flux - error.lower for flux, error in zipped]))
                upper_flux_model = max(filter(None, [flux + error.upper for flux, error in zipped]))
                if lower_flux_model < self._min_flux or self._min_flux is None: self._min_flux = lower_flux_model
                if upper_flux_model > self._max_flux or self._max_flux is None: self._max_flux = upper_flux_model
            else:
                # Keep track of the minimal and maximal flux
                min_flux_model = min(filter(None, fluxes))  # ignores zeros; see http://stackoverflow.com/questions/21084714/find-the-lowest-value-that-is-not-null-using-python
                max_flux_model = max(filter(None, fluxes))  # idem.
                if min_flux_model < self._min_flux or self._min_flux is None: self._min_flux = min_flux_model
                if max_flux_model > self._max_flux or self._max_flux is None: self._max_flux = max_flux_model

            # Keep track of the minimal and maximal wavelength
            min_wavelength_model = min(wavelengths)
            max_wavelength_model = max(wavelengths)
            if min_wavelength_model < self._min_wavelength or self._min_wavelength is None: self._min_wavelength = min_wavelength_model
            if max_wavelength_model > self._max_wavelength or self._max_wavelength is None: self._max_wavelength = max_wavelength_model

        # Plot the data points (as a line) on the axis
        log_model = np.log10(fluxes)
        axis.plot(wavelengths, log_model, ls=linestyle, color=linecolor, label=label)

        # Plot errors
        if errors is not None:

            errors_wavelengths = []

            bottom = []
            top = []
            for j in range(len(errors)):

                if errors[j] is None: continue

                errors_wavelengths.append(wavelengths[j])

                value = fluxes[j]
                error = errors[j]

                bottom.append(value + error.lower)
                top.append(value + error.upper)

            bottom = np.array(bottom)
            top = np.array(top)

            log_bottom = np.log10(bottom)
            log_top = np.log10(top)

            axis.fill_between(errors_wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor=errorcolor, edgecolor=errorcolor, interpolate=True, alpha=0.5)
            axis.plot([], [], color=errorcolor, linewidth=10, label='spread')

    # -----------------------------------------------------------------

    def finish_plot(self, for_legend_patches=None, for_legend_parameters=None, extra_legend=None):

        """
        This function ...
        :param for_legend_patches:
        :param for_legend_parameters:
        :param extra_legend:
        :return:
        """

        # Debugging
        log.debug("Finishing the plot ...")

        had_min_flux = self.min_flux is not None
        had_max_flux = self.max_flux is not None
        had_min_wavelength = self.min_wavelength is not None
        had_max_wavelength = self.max_wavelength is not None

        # Axis limits are now definite
        if self.min_flux is None: self.min_flux = self._min_flux
        if self.max_flux is None: self.max_flux = self._max_flux
        if self.min_wavelength is None: self.min_wavelength = self._min_wavelength
        if self.max_wavelength is None: self.max_wavelength = self._max_wavelength

        # Adjust the limits by adding 10% at each side
        factor_x = 10 ** (0.1 * np.log10(self.max_wavelength / self.min_wavelength))
        factor_y = 10 ** (0.1 * np.log10(self.max_flux / self.min_flux))

        if not had_min_flux: self.min_flux /= factor_y
        if not had_max_flux: self.max_flux *= factor_y
        if not had_min_wavelength: self.min_wavelength /= factor_x
        if not had_max_wavelength: self.max_wavelength *= factor_x

        # Format residual axes
        for res_axis in self.residual_plots:

            # Set linestyle and limit for axis2
            res_axis.axhline(y=0., color='black', ls='-.')
            res_axis.set_ylim(-95, 95)

            # Set axis label
            res_axis.set_xscale('log')
            res_axis.set_ylabel(r"Res. $[\%]$", fontsize='large')

        # Set x label of the last residual plot
        if len(self.residual_plots) > 0: self.residual_plots[len(self.residual_plots)-1].set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')
        else: self.main_plot.set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')

        # Set log x scale
        self.main_plot.set_xscale('log')

        # Set ticks
        self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        self.main_plot.set_yticks(fontsize=self.config.plot.ticks_fontsize) # minor = False is possible

        # Add axis labels and a legend
        y_label = r"Log $" + self.config.unit.latex_string + r"$"
        self.main_plot.set_ylabel(y_label, fontsize='large')

        # Set grid
        #self.figure.set_grid(self.config.plot, which="both")
        self.main_plot.set_grid(self.config.plot, which="both")
        for plot in self.residual_plots: plot.set_grid(self.config.plot, which="both")

        # Set borders
        self.figure.set_borders(self.config.plot)

        # Set flux axis limits
        plot_min, plot_max = get_plot_flux_limits(self.min_flux, self.max_flux)
        self.main_plot.set_ylim((plot_min, plot_max))

        # Set wavelength axis limits
        #plot_min_wavelength, plot_max_wavelength = get_plot_wavelength_limits(self.min_wavelength, self.max_wavelength)
        #self._main_axis.set_xlim(plot_min_wavelength, plot_max_wavelength)

        # Add the plots to the figure
        if self.config.library == bokeh:
            self.figure.add_column(self.main_plot, *self.residual_plots)

        legends = []

        # Add the legend
        if for_legend_patches is not None:

            legend_properties = dict()
            legend_properties["loc"] = "lower center"
            #legend_properties["numpoints"] = 4
            legend_properties["scatterpoints"] = 4
            legend_properties["ncol"] = 2
            legend_properties["shadow"] = False
            legend_properties["frameon"] = True
            legend_properties["facecolor"] = None
            legend_properties["fontsize"] = self.config.plot.legend_fontsize

            # Set legend
            # fancybox=True makes the legend corners rounded
            legend = self.main_plot.legend([l[0] for l in for_legend_patches], for_legend_parameters, **legend_properties)
            legends.append(legend)

            # Extra legend
            if extra_legend is not None:

                self.main_plot.add_artist(extra_legend)
                legends.append(extra_legend)

        # No extra legend, no patches for the legend
        else:

            legend_properties = dict()
            legend_properties["loc"] = "lower center"
            legend_properties["numpoints"] = 1
            legend_properties["scatterpoints"] = 4
            legend_properties["ncol"] = 2
            legend_properties["shadow"] = False
            legend_properties["frameon"] = True
            legend_properties["facecolor"] = None
            legend_properties["fontsize"] = self.config.plot.legend_fontsize

            legend = self.main_plot.legend(**legend_properties)
            legends.append(legend)

        # Set legend properties
        #for legend in legends:
        #    if legend is None: continue # Bokeh
        #    set_legend(legend, self.config.plot)

        # Add title if requested
        if self.title is not None: self.figure.set_title(self.title) #self.figure.figure.suptitle("\n".join(wrap(self.title, 60)))

        # Save or show the plot
        #if self.out_path is None: self.figure.show()
        if self.config.show: self.figure.show()

        # Save the figure
        if self.out_path is not None: self.save_figure()

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Close the figure
        if self.figure is not None: self.figure.close()

    # -----------------------------------------------------------------

    def save_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Saving figure ...")

        # Determine path
        if types.is_string_type(self.out_path):
            if fs.is_directory(self.out_path): path = fs.join(self.out_path, "seds" + self.config.format)
            else: path = self.out_path
        else: path = self.out_path

        # Save
        self.figure.saveto(path)

# -----------------------------------------------------------------

def set_legend(legend, config):

    """
    This function ...
    :param legend:
    :param config:
    :return:
    """

    # if nlegends > 1: percentage = 25.
    # else: percentage = 10.
    #percentage = 20.

    # Shrink current axis's height by a certain percentage on the bottom
    #box = ax.get_position()
    #ax.set_position(
    #    [box.x0, box.y0 + box.height * percentage / 100., box.width, box.height * (100 - percentage) / 100.])

    # Plot legend
    #legend_title = r"\underline{" + legend_title + "}"
    #legend = plt.legend(loc="upper center", title=legend_title, bbox_to_anchor=(0.5, -0.25), fancybox=False,
    #                    shadow=False, ncol=4)

    frame = legend.get_frame()

    # Set legend frame color and line width
    if config.add_legend_border:
        frame.set_linewidth(config.legend_borderwidth)
    else:
        frame.set_linewidth(0)
    frame.set_edgecolor(config.legend_bordercolor)

    # Set background color
    frame.set_facecolor('0.85')
    legend.legendPatch.set_alpha(0.75)

    index = 0

    # Move to foreground
    legend.set_zorder(100 + index)

    # Set fontsize
    plt.setp(legend.get_title(), fontsize=str(config.legend_title_fontsize))

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
    :param reference_instruments:
    :param reference_bands:
    :param reference_wavelengths:
    :param reference_fluxes:
    :return:
    """

    # Loop over the entries in the reference SED
    for i in range(len(reference_wavelengths)):

        if (reference_instruments[i] == instrument and reference_bands[i] == band) or np.isclose(reference_wavelengths[i], wavelength):
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

def process_errorbar(error):

    """
    This function ...
    :param error:
    :return:
    """

    lolims = False
    uplims = False

    error_limits = [[abs(error.lower), abs(error.upper)]]

    if error.no_lower:

        error_limits[0][0] = error_limits[0][1]
        uplims = [[True]]

    if error.no_upper:

        error_limits[0][1] = error_limits[0][0]
        lolims = [[True]]

    yerr = np.array(error_limits).T
    return yerr, lolims, uplims

# -----------------------------------------------------------------
