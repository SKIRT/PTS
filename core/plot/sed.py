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
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
from collections import OrderedDict
from matplotlib import rc, rcdefaults
from scipy.interpolate import interp1d
from operator import itemgetter
from matplotlib import lines
from scipy.interpolate import spline

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..data.sed import ObservedSED, SED, load_sed, load_multiple_seds, is_sed
from ..basics.plot import MPLFigure, BokehFigure, BokehPlot, mpl, bokeh, dark_pretty_colors, pretty_colors, filled_markers
from ..tools import types
from ..filter.broad import BroadBandFilter
from ..basics.errorbar import ErrorBar
from ..tools import numbers
from ..tools.utils import lazyproperty
from ..tools import introspection
from ..tools import sequences
from ..basics.map import Map

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

models_reference = "models"
observations_reference = "observations"
residual_references = [models_reference, observations_reference]

# -----------------------------------------------------------------

def plot_sed(sed, label=None, path=None, title=None, show_file=False, format="pdf", unit=None,
             wavelength_unit=None, distance=None):

    """
    This function ...
    :param sed:
    :param label:
    :param path:
    :param title:
    :param show_file:
    :param format:
    :param unit:
    :param wavelength_unit:
    :parma distance:
    :return:
    """

    # Create a new SEDPlotter instance
    plotter = SEDPlotter()

    # Set units
    if wavelength_unit is None: wavelength_unit = sed.wavelength_unit
    if unit is None: unit = sed.unit
    plotter.config.wavelength_unit = wavelength_unit
    plotter.config.unit = unit
    #print(plotter.config.unit)

    # Set distance?
    if distance is not None: plotter.config.distance = distance

    # Determine label
    if label is None:
        if isinstance(sed, ObservedSED): label = "observation"
        elif isinstance(sed, SED): label = "model"
        else: raise ValueError("Invalid SED object")

    # Add the SED
    plotter.add_sed(sed, label)

    # Set filepath, if plot is to be shown as file
    if path is None and show_file: path = fs.join(introspection.pts_temp_dir, "seds." + format)

    # Make the plot
    plotter.run(title=title, output=path)

    # Show file?
    if show_file: fs.open_file(path)

    # Return the plotter
    return plotter

# -----------------------------------------------------------------

def plot_seds_quick(**seds):

    """
    This function ...
    :param seds:
    :return:
    """

    # Plot
    plot_seds(seds)

# -----------------------------------------------------------------

def plot_seds(seds, **kwargs):

    """
    This function ...
    :param seds:
    :param kwargs:
    :return:
    """

    # Get settings
    title = kwargs.pop("title", None)
    path = kwargs.pop("path", None)
    show_file = kwargs.pop("show_file", None)
    format = kwargs.pop("format", "pdf")
    residuals = kwargs.pop("residuals", True)
    ghost = kwargs.pop("ghost", False)
    wavelength_unit = kwargs.pop("wavelength_unit", None)
    unit = kwargs.pop("unit", None)
    distance = kwargs.pop("distance", None)
    options = kwargs.pop("options", {})
    tex = kwargs.pop("tex", True)

    # Figure and plots
    figure = kwargs.pop("figure", None)
    main_plot = kwargs.pop("main_plot", None)
    residual_plot = kwargs.pop("residual_plot", None)
    residual_plots = kwargs.pop("residual_plots", None)

    # Position of axis ticks and labels
    xaxis_position = kwargs.pop("xaxis_position", "bottom")
    yaxis_position = kwargs.pop("yaxis_position", "left")

    # Residual reference
    residual_reference = kwargs.pop("residual_reference", "models")

    # Smooth residuals
    smooth_residuals = kwargs.pop("smooth_residuals", False)

    # Only residuals on separate legend in the residuals panel
    only_residuals_legend = kwargs.pop("only_residuals_legend", False)

    # Show the figure after plotting
    show = kwargs.pop("show", None)

    # Other options
    show_models_residuals = kwargs.pop("show_models_residuals", False)
    show_smooth = kwargs.pop("show_smooth", False)
    smooth_residuals_ignore_close = kwargs.pop("smooth_residuals_ignore_close", False)

    # LEGENDS
    instruments_legend_ncols = kwargs.pop("instruments_legend_ncols", None)
    observations_legend_ncols = kwargs.pop("observations_legend_ncols", None)
    models_legend_ncols = kwargs.pop("models_legend_ncols", None)
    instruments_legend_location = kwargs.pop("instruments_legend_location", None)
    observations_legend_location = kwargs.pop("observations_legend_location", None)
    models_legend_location = kwargs.pop("models_legend_location", None)
    observations_residuals_legend_location = kwargs.pop("observations_residuals_legend_location", None)
    models_residuals_legend_location = kwargs.pop("models_residuals_legend_location", None)
    observations_residuals_legend_ncols = kwargs.pop("observations_residuals_legend_ncols", None)
    models_residuals_legend_ncols = kwargs.pop("models_residuals_legend_ncols", None)

    show_models_legend = kwargs.pop("show_models_legend", True)
    show_observations_legend = kwargs.pop("show_observations_legend", True)
    show_instruments_legend = kwargs.pop("show_instruments_legend", True)
    show_residuals_legends = kwargs.pop("show_residuals_legends", True)

    # Create SED plotter
    plotter = SEDPlotter(kwargs)

    # Legend props
    if instruments_legend_ncols is not None: plotter.config.legends.instruments_ncols = instruments_legend_ncols
    if observations_legend_ncols is not None: plotter.config.legends.observations_ncols = observations_legend_ncols
    if models_legend_ncols is not None: plotter.config.legends.models_ncols = models_legend_ncols
    if instruments_legend_location is not None: plotter.config.legends.instruments_location = instruments_legend_location
    if observations_legend_location is not None: plotter.config.legends.observations_location = observations_legend_location
    if models_legend_location is not None: plotter.config.legends.models_location = models_legend_location
    if observations_residuals_legend_ncols is not None: plotter.config.legends.observations_residuals_ncols = observations_residuals_legend_ncols
    if models_residuals_legend_ncols is not None: plotter.config.legends.models_residuals_ncols = models_residuals_legend_ncols
    if observations_residuals_legend_location is not None: plotter.config.legends.observations_residuals_location = observations_residuals_legend_location
    if models_residuals_legend_location is not None: plotter.config.legends.models_residuals_location = models_residuals_legend_location

    # Set other options
    plotter.config.legends.instruments = show_instruments_legend
    plotter.config.legends.observations = show_observations_legend
    plotter.config.legends.models = show_models_legend
    plotter.config.legends.residuals = show_residuals_legends
    plotter.config.models_residuals = show_models_residuals
    plotter.config.show_smooth = show_smooth
    plotter.config.smooth_residuals_ignore_close = smooth_residuals_ignore_close

    # Set TeX flag
    plotter.config.tex = tex

    # Set other flags
    plotter.config.smooth_residuals = smooth_residuals
    plotter.config.only_residuals_legend = only_residuals_legend

    # Set residual reference
    plotter.config.residual_reference = residual_reference

    # Set axis positions
    plotter.config.xaxis_position = xaxis_position
    plotter.config.yaxis_position = yaxis_position

    # Set show flag (None means automatically determined)
    plotter.config.show = show

    # Get the units
    if wavelength_unit is None:
        wavelength_units = [seds[label].wavelength_unit for label in seds]
        #print(wavelength_units)
        wavelength_unit = sequences.get_single(wavelength_units, none="error", method="first_not_none")
        #print(wavelength_unit)
    if unit is None:
        units = [seds[label].unit for label in seds]
        unit = sequences.get_single(units)

    # Set units
    #print("WAVELENGTH UNIT:", wavelength_unit)
    #print("UNIT:", unit)
    plotter.config.wavelength_unit = wavelength_unit
    plotter.config.unit = unit

    # Set distance
    if distance is not None: plotter.config.distance = distance

    # Add SEDs
    plotter.add_seds(seds, residuals=residuals, ghost=ghost, options=options)

    # Set filepath, if plot is to be shown as file
    if path is None and show_file: path = fs.join(introspection.pts_temp_dir, "seds." + format)

    # Set minima and maxima
    plotter.config.min_wavelength = kwargs.pop("min_wavelength", None)
    plotter.config.max_wavelength = kwargs.pop("max_wavelength", None)
    plotter.config.min_flux = kwargs.pop("min_flux", None)
    plotter.config.max_flux = kwargs.pop("max_flux", None)

    # Run the plotter
    plotter.run(title=title, output=path, figure=figure, main_plot=main_plot, residual_plot=residual_plot, residual_plots=residual_plots)

    # Show file
    if show_file: fs.open_file(path)

    # Return the plotter
    return plotter

# -----------------------------------------------------------------

def get_sed_template(name, **kwargs):

    """
    This function ...
    :param name:
    :return:
    """

    from ...modeling.core.mappings import create_mappings_sed
    from ...modeling.core.bruzualcharlot import create_bruzual_charlot_sed

    # Mappings
    if name == "mappings": return create_mappings_sed(**kwargs)

    # Bruzual-Charlot
    elif name == "bruzual_charlot": return create_bruzual_charlot_sed(**kwargs)

    # Not recognized
    else: raise ValueError("Template not recognized")

# -----------------------------------------------------------------

class Observation(object):

    """
    This class ...
    """

    def __init__(self):

        self.wavelengths = []
        self.fluxes = []
        self.errors = []
        self.instruments = []
        self.bands = []
        self.labels = []
        self.descriptions = []
        self.markers = []

# -----------------------------------------------------------------

class Model(object):

    """
    This class ...
    """

    def __init__(self):

        self.wavelengths = []
        self.fluxes = []
        self.errors = []
        self.above_wavelengths = []
        self.above_fluxes = []

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
        self.models = OrderedDict()
        self.observations = OrderedDict()

        # After processing
        self._observations = defaultdict(Observation)
        self._models = defaultdict(Model)

        # Options for models and observations
        self.model_options = defaultdict(Map)
        self.observation_options = defaultdict(Map)

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
        self._external_figure = False

        # The main plot and the residual plots
        self.main_plot = None
        self.residual_plots = []

        # Legends
        self.legends = []
        self.legends_residuals = []

        # FOR INSTRUMENTS LEGEND
        self.instruments_legend_patches = []
        self.instruments_legend_labels = []

        # FOR OBSERVATIONS LEGEND
        self.observations_legend_patches = []
        self.observations_legend_labels = []

        # FOR MODELS LEGEND
        self.models_legend_patches = []
        self.models_legend_labels = []

        # Used labels
        self.used_labels = []

        # Remember colors and linestyles for models
        self.model_colors = dict()
        self.model_styles = dict()

        # Remember colors for each observation
        self.observation_colors = dict()
        self.observation_markers = None

        # Dictionary that keeps data about filling between curves
        self.fills = dict()

    # -----------------------------------------------------------------

    @property
    def has_figure(self):
        return self.figure is not None

    # -----------------------------------------------------------------

    @property
    def has_plots(self):
        return self.has_main_plot

    # -----------------------------------------------------------------

    @property
    def has_main_plot(self):
        return self.main_plot is not None

    # -----------------------------------------------------------------

    @property
    def nresidual_plots(self):
        return len(self.residual_plots)

    # -----------------------------------------------------------------

    @property
    def has_residual_plots(self):
        return self.nresidual_plots > 0

    # -----------------------------------------------------------------

    @property
    def has_single_residual_plot(self):
        return self.nresidual_plots == 1

    # -----------------------------------------------------------------

    @property
    def residual_plot(self):
        if not self.has_single_residual_plot: raise RuntimeError("This setup does not have one residual plot pane or plots not yet created")
        return self.residual_plots[0]

    # -----------------------------------------------------------------

    @residual_plot.setter
    def residual_plot(self, plot):
        if self.has_residual_plots: raise ValueError("Cannot set residual plot when residual plots are already created")
        self.residual_plots[0] = plot

    # -----------------------------------------------------------------

    @property
    def nlegends(self):
        return len(self.legends)

    # -----------------------------------------------------------------

    @property
    def has_legends(self):
        return self.nlegends > 0

    # -----------------------------------------------------------------

    @property
    def nlegends_residuals(self):
        return len(self.legends_residuals)

    # -----------------------------------------------------------------

    @property
    def has_legends_residuals(self):
        return self.nlegends_residuals > 0

    # -----------------------------------------------------------------

    #@lazyproperty
    @property # HAS TO BE PROPERTY BECAUSE OF CHECKS IN SETUP
    def nseds(self):
        return self.nmodels + self.nobservations

    # -----------------------------------------------------------------

    #@lazyproperty
    @property # HAS TO BE PROPERTY BECAUSE OF CHECKS IN SETUP
    def no_seds(self):
        return self.nseds == 0

    # -----------------------------------------------------------------

    #@lazyproperty
    @property # HAS TO BE PROPERTY BECAUSE OF CHECKS IN SETUP
    def nmodels(self):
        return len(self.models)

    # -----------------------------------------------------------------

    @lazyproperty
    def one_model(self):
        return self.nmodels == 1

    # -----------------------------------------------------------------

    @lazyproperty
    def has_models(self):
        return self.nmodels > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def has_multiple_models(self):
        return self.nmodels > 1

    # -----------------------------------------------------------------

    @lazyproperty
    def no_models(self):
        return self.nmodels == 0

    # -----------------------------------------------------------------

    #@lazyproperty
    @property # HAS TO BE PROPERTY BECAUSE OF CHECKS IN SETUP
    def nobservations(self):
        return len(self.observations)

    # -----------------------------------------------------------------

    @lazyproperty
    def one_observation(self):
        return self.nobservations == 1

    # -----------------------------------------------------------------

    @lazyproperty
    def has_observations(self):
        return self.nobservations > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def has_multiple_observations(self):
        return self.nobservations > 1

    # -----------------------------------------------------------------

    @lazyproperty
    def no_observations(self):
        return self.nobservations == 0

    # -----------------------------------------------------------------

    def add_observation(self, sed, label, **kwargs):

        """
        This function ...
        :param sed:
        :param label:
        :param kwargs:
        :return:
        """

        # Check label
        if label in self.observation_labels: raise ValueError("Already an observation with the label '" + label + "'")

        # Make a copy with only broad band filters
        only_broad = sed.only_broad_band()

        # Check whether there are any points in the SED
        if len(only_broad) == 0: log.warning("The '" + label + "' SED contains no (broad band) photometry points")

        # Set the SED
        self.observations[label] = only_broad

        # Set options
        self.observation_options[label].only_residuals = kwargs.pop("only_residuals", False)
        self.observation_options[label].as_reference = kwargs.pop("as_reference", True) # will be used as reference if residual_reference = "observations"
        self.observation_options[label].color = kwargs.pop("color", None)
        self.observation_options[label].join_residuals = kwargs.pop("join_residuals", None)

    # -----------------------------------------------------------------

    def set_observation_options(self, label, **kwargs):
        self.observation_options[label].update(**kwargs)

    # -----------------------------------------------------------------

    def add_model(self, sed, label, **kwargs):

        """
        This function ...
        :param sed:
        :param label:
        :param kwargs:
        :return:
        """

        # Check label
        if label in self.model_labels: raise ValueError("Already a model with the label '" + label + "'")

        # Add the SED
        self.models[label] = sed

        # Set options
        self.model_options[label].only_residuals = kwargs.pop("only_residuals", False)
        self.model_options[label].residuals = kwargs.pop("residuals", True)
        self.model_options[label].ghost = kwargs.pop("ghost", False)
        self.model_options[label].above = kwargs.pop("above", None)
        self.model_options[label].above_name = kwargs.pop("above_name", None)
        self.model_options[label].fill = kwargs.pop("fill", True) # fill if when above other
        self.model_options[label].color = kwargs.pop("color", None) # None means automatic
        self.model_options[label].linestyle = kwargs.pop("linestyle", None) # None means automatic
        self.model_options[label].residual_color = kwargs.pop("residual_color", None) # None means automatic

    # -----------------------------------------------------------------

    def set_model_options(self, label, **kwargs):
        self.model_options[label].update(**kwargs)

    # -----------------------------------------------------------------

    def add_sed(self, sed, label, **kwargs):

        """
        This function ...
        :param sed:
        :param label:
        :param residuals: whether plotting the residual curve makes sense for this SED
        (it does not when the SED is only for one component, or a contribution to the total SED)
        :param ghost:
        :return:
        """

        # Debugging
        log.debug("Adding the '" + label + "' SED ...")

        # Add observed or modeled SED
        if isinstance(sed, ObservedSED): self.add_observation(sed, label, **kwargs)

        # Add model SED
        elif isinstance(sed, SED): self.add_model(sed, label, **kwargs)

        # Invalid
        else: raise ValueError("The SED must be an SED or ObservedSED instance")

    # -----------------------------------------------------------------

    def add_seds(self, seds, **kwargs):

        """
        This function ...
        :param seds:
        :param options:
        :return:
        """

        # Get options
        options = kwargs.pop("options", {})
        if options is None: options = {}

        # Loop over the SEDs
        for label in seds:

            # Are options given for this SED individually?
            if label in options: sed_options = options[label]
            else: sed_options = kwargs

            # Add
            self.add_sed(seds[label], label, **sed_options)

    # -----------------------------------------------------------------

    @property
    def do_create_plots(self):
        return not self.has_plots

    # -----------------------------------------------------------------

    @property
    def do_set_defaults(self):
        return not self._external_figure

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create the plots if necessary
        if self.do_create_plots: self.create_plots()

        # Plot the observations
        self.plot_observations()

        # Plot the models
        self.plot_models()

        # Plot the residusl
        self.plot_residuals()

        # Create the legends
        self.create_legends()

        # Finish the plot
        self.finish_plot()

        # Reset matplotlib defaults
        if self.do_set_defaults: self.set_defaults()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDPlotter, self).setup(**kwargs)

        # Use LaTeX rendering
        if self.config.tex: rc('text', usetex=True)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the axis limits
        self.min_wavelength = kwargs.pop("min_wavelength", None)
        self.max_wavelength = kwargs.pop("max_wavelength", None)
        self.min_flux = kwargs.pop("min_flux", None)
        self.max_flux = kwargs.pop("max_flux", None)

        # Set from config
        if self.min_wavelength is None: self.min_wavelength = self.config.min_wavelength
        if self.max_wavelength is None: self.max_wavelength = self.config.max_wavelength
        if self.min_flux is None: self.min_flux = self.config.min_flux
        if self.max_flux is None: self.max_flux = self.config.max_flux

        # Remove units from min and max wavelength
        if self.min_wavelength is not None and hasattr(self.min_wavelength, "unit"): self.min_wavelength = self.min_wavelength.to(self.config.wavelength_unit).value
        if self.max_wavelength is not None and hasattr(self.max_wavelength, "unit"): self.max_wavelength = self.max_wavelength.to(self.config.wavelength_unit).value

        # Remove units from min and max flux
        if self.min_flux is not None and hasattr(self.min_flux, "unit"):
            conv_info = dict(**self.conversion_info)
            if self.min_wavelength is not None: conv_info["wavelength"] = self.min_wavelength * self.config.wavelength_unit
            self.min_flux = self.min_flux.to(self.config.unit, **conv_info).value
        if self.max_flux is not None and hasattr(self.max_flux, "unit"):
            conv_info = dict(**self.conversion_info)
            if self.max_wavelength is not None: conv_info["wavelength"] = self.max_wavelength * self.config.wavelength_unit
            self.max_flux = self.max_flux.to(self.config.unit, **conv_info).value

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

        # STILL NO SEDS?
        if self.no_seds: raise RuntimeError("There are no SEDs loaded or found in the input directory")

        # Add additional errors
        if self.config.additional_error is not None: self.add_relative_errors()

        # Add SED templates
        if self.config.add_templates: self.create_templates()

        # Create the figure
        if kwargs.get("figure", None) is not None:
            self.figure = kwargs.pop("figure")
            self._external_figure = True
        else:
            if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
            elif self.config.library == bokeh: self.figure = BokehFigure()
            else: raise ValueError("Invalid libary: " + self.config.library)

        # Plots are passed
        if kwargs.get("main_plot", None) is not None: self.main_plot = kwargs.pop("main_plot")
        if kwargs.get("residual_plot", None) is not None: self.residual_plot = kwargs.pop("residual_plot")
        elif kwargs.get("residual_plots", None) is not None: self.residual_plots = kwargs.pop("residual_plots")

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

    # -----------------------------------------------------------------

    def create_templates(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Creating SED templates ...")

        # Loop over the template names
        for name in self.config.templates:

            # MAPPINGS
            if name == "mappings":

                # Debugging
                log.debug("Creating MAPPINGS SED template ...")

                properties = dict()
                properties["metallicity"] = self.config.metallicity
                properties["compactness"] = self.config.compactness
                properties["pressure"] = self.config.pressure
                properties["covering_factor"] = self.config.covering_factor

                # Set label
                label = "MAPPINGS"

                # Create the template SED
                sed = get_sed_template(name, **properties)

                # Add the SED
                self.add_sed(sed, label=label)

            # Stellar Bruzual Charlot
            elif name == "bruzual_charlot":

                # Debugging
                log.debug("Creating Bruzual-Charlot SED templates ...")

                # Loop over the ages
                for age in self.config.ages:

                    properties = dict()
                    properties["metallicity"] = self.config.metallicity
                    properties["age"] = age

                    # Create the SED template
                    label = "Bruzual-Charlot " + str(age)
                    sed = get_sed_template(name, **properties)

                    # Add the SED
                    self.add_sed(sed, label=label)

            # Invalid
            else: raise ValueError("Invalid SED template name")

    # -----------------------------------------------------------------

    @property
    def xsize(self):
        return self.config.plot.xsize

    # -----------------------------------------------------------------

    @lazyproperty
    def ysize(self):
        if self.nobservations > 1 and self.nmodels > 0:
            main_size = self.config.plot.ysize * self.config.plot.main_relsize
            other_size = float(self.config.plot.ysize) - main_size
            return main_size + (self.nobservations - 1) * other_size
        else: return self.config.plot.ysize

    # -----------------------------------------------------------------

    @lazyproperty
    def figsize(self):
        return (self.xsize, self.ysize,)

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

        # Get a dictionary
        if types.is_dictionary(self.config.seds): sed_paths = self.config.seds
        elif types.is_string_sequence(self.config.seds):
            sed_paths = OrderedDict() # initialize
            for filepath in self.config.seds:
                name = fs.strip_extension(fs.name(filepath))
                sed_paths[name] = filepath
        else: raise ValueError("Invalid type for 'seds': must be dictionary or sequence")

        # Loop over the SED files
        for name in sed_paths:

            # Get the SED filepath
            path = sed_paths[name]

            # Try getting multiple SEDs
            if self.config.multi:

                seds = load_multiple_seds(path, as_dict=True, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)
                for label in seds: self.add_sed(seds[label], label=name + " " + label)

            # One SED per file
            else:

                #print(path)

                # Load SED
                sed = load_sed(path, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)
                #print(sed)

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
        for path, name in fs.files_in_path(self.config.path, extension="dat", returns=["path", "name"],
                                           contains=self.config.contains, not_contains=self.config.not_contains,
                                           exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                           startswith=self.config.startswith, endswith=self.config.endswith):

            # Skip emission lines
            if "Lines" in name: continue

            # Check whether this is an SED file
            if not is_sed(path):
                log.warning("Ignoring file '" + path + "' because it is probably not an SED file ...")
                continue
            #print(path)

            # Try getting multiple SEDs
            if self.config.multi:

                seds = load_multiple_seds(path, as_dict=True, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)
                for label in seds: self.add_sed(seds[label], label=name + " " + label)

            # One SED per file
            else:

                # Load SED
                sed = load_sed(path, wavelength_unit=self.config.wavelength_unit_file, photometry_unit=self.config.unit_file)

                # Add the SED
                self.add_sed(sed, label=name)

    # -----------------------------------------------------------------

    def add_relative_errors(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adding relative errors ...")

        # Loop over the observed SEDs
        for label in self.observations:

            # Check
            if self.config.additional_for is not None and label not in self.config.additional_for: continue
            if self.config.additional_not_for is not None and label in self.config.additional_not_for: continue

            # Debugging
            log.debug("Adding additional relative errors to the data points of the '" + label + "' observed SED ...")

            # Get the observed SED, make copy
            observation = self.observations[label].copy()

            # Add the relative error
            observation.add_or_set_relative_error(self.config.additional_error)

            # Set
            self.observations[label] = observation

    # -----------------------------------------------------------------

    def clear(self):
        raise NotImplementedError("Clearing the SED plotter is no longer supported. Create a new instance for new plots.")

    # -----------------------------------------------------------------

    def create_plots(self):

        """
        This function ...
        :return:
        """

        # No models, just observed SEDs
        if self.no_models:

            # Only one observation: just create a main plot
            if self.one_observation: nrespanels = 0

            # Multiple observations
            elif self.config.observations_residuals: nrespanels = 1

            # Don't plot residuals
            else: nrespanels = 0

        # With models
        else:

            # Only models
            if self.no_observations:

                # Create
                if self.config.models_residuals and self.nmodels_for_residuals > 1: nrespanels = 1
                else: nrespanels = 0

            # One observation
            elif self.one_observation:

                # Determine number of residual panels
                if self.residual_reference == observations_reference: nrespanels = 1
                elif self.residual_reference == models_reference: nrespanels = self.nmodels_for_residuals
                else: raise ValueError("Invalid residual reference '" + self.residual_reference + "'")

            # Multiple observations
            else:

                # Determine number of plot rows
                if self.residual_reference == observations_reference: nrespanels = self.nobservations_as_reference
                elif self.residual_reference == models_reference: nrespanels = self.nmodels_for_residuals
                else: raise ValueError("Invalid residual reference '" + self.residual_reference + "'")

        # Create
        self.main_plot, self.residual_plots = self.figure.create_sed_plots(nresiduals=nrespanels)

    # -----------------------------------------------------------------

    @lazyproperty
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
            errors = observation.errors()

            # Has errors?
            has_errors = sequences.has_not_none(errors)

            # Certainly one observation with errors
            if has_errors: return True

        # No observation with errors
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_observation_point_labels(self):
        return get_unique_labels(self.observations)

    # -----------------------------------------------------------------

    @property
    def observation_labels(self):
        return self.observations.keys()

    # -----------------------------------------------------------------

    @property
    def model_labels(self):
        return self.models.keys()

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
        return self.observations[self.most_npoints_observation_label]

    # -----------------------------------------------------------------

    def plot_with_models(self):

        """
        This function ...
        :return:
        """

        # Only models
        if self.no_observations: pass #self.plot_only_models()

        # One observation
        elif self.one_observation: self.plot_one_observation_with_models()

        # Multiple observations
        else: self.plot_new()

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_labels_as_reference(self):
        return [obs_label for obs_label in self.observations if self.observation_options[obs_label].as_reference]

    # -----------------------------------------------------------------

    @lazyproperty
    def nobservations_as_reference(self):
        return len(self.observation_labels_as_reference)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_labels_for_residuals(self):
        return [model_label for model_label in self.models if self.model_options[model_label].residuals]

    # -----------------------------------------------------------------

    @lazyproperty
    def nmodels_for_residuals(self):
        return len(self.model_labels_for_residuals)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_labels_not_ghost(self):
        return [model_label for model_label in self.models if not self.model_options[model_label].ghost]

    # -----------------------------------------------------------------

    @lazyproperty
    def nmodels_not_ghost(self):
        return len(self.model_labels_not_ghost)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_labels_not_ghost_not_residuals(self):
        return [model_label for model_label in self.models if not self.model_options[model_label].ghost and not self.model_options[model_label].residuals]

    # -----------------------------------------------------------------

    @lazyproperty
    def nmodels_not_ghost_not_residuals(self):
        return len(self.model_labels_not_ghost_not_residuals)

    # -----------------------------------------------------------------

    @lazyproperty
    def conversion_info(self):
        info = dict()
        if self.config.distance is not None: info["distance"] = self.config.distance
        return info

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_labels_to_join(self):
        return [label for label in self.observation_labels if self.observation_options[label].join_residuals is not None]

    # -----------------------------------------------------------------

    @lazyproperty
    def join_observations_with(self):
        join_with = defaultdict(list)
        for label in self.observation_labels:
            for other_label in self.observation_labels_to_join:
                if self.observation_options[other_label].join_residuals == label: join_with[label].append(other_label)
        return join_with

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_labels_with_joined(self):
        return self.join_observations_with.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def joined_observations(self):
        from ..basics.containers import create_subdict
        joined = dict()
        for label in self.observation_labels_with_joined:
            joined_with_labels = self.join_observations_with[label]
            all_labels = [label] + joined_with_labels
            seds = create_subdict(self.observations, all_labels)
            sed = join_observed_seds(seds)
            joined[label] = sed
        return joined

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_labels_only_residuals(self):
        return [label for label in self.observation_labels if self.observation_options[label].only_residuals]

    # -----------------------------------------------------------------

    @lazyproperty
    def nobservations_only_residuals(self):
        return len(self.observation_labels_only_residuals)

    # -----------------------------------------------------------------

    @property
    def has_observations_only_residuals(self):
        return self.nobservations_only_residuals > 0

    # -----------------------------------------------------------------

    def create_marker_patch(self, label, marker, color, edgecolor="black", edgewidth=1, size=7):

        """
        This function ...
        :param marker:
        :param color
        :param edgecolor:
        :param edgewidth:
        :param size:
        :param label:
        :return:
        """

        # Create and return
        return lines.Line2D([], [], marker=marker, markersize=size, label=label, linewidth=0, markeredgecolor=edgecolor, markerfacecolor=color, markeredgewidth=edgewidth)

    # -----------------------------------------------------------------

    def create_line_patch(self, label, color, style="-"):

        """
        This function ...
        :param label:
        :param color:
        :param style:
        :return:
        """

        # Create and return
        return lines.Line2D([], [], label=label, color=color, linestyle=style)

    # -----------------------------------------------------------------

    def create_rectangle_patch(self, label, color, edgecolor="black"):

        """
        This function ...
        :param label:
        :param color:
        :param edgecolor:
        :return:
        """

        # Create and return
        return patches.Rectangle((0, 0), 1, 1, fc=color, edgecolor=edgecolor, label=label)

    # -----------------------------------------------------------------

    def create_marker_legend(self, specs, props=None, plot=None):

        """
        This function ...
        :param specs:
        :param props:
        :param plot:
        :return:
        """

        patches = []
        for label in specs:
            marker, color = specs[label]
            patch = self.create_marker_patch(label, marker, color)
            patches.append(patch)

        # Return the legend
        if props is None: props = {}
        if plot is None: plot = self.main_plot
        return plot.create_legend(lines, specs.keys(), **props)

    # -----------------------------------------------------------------

    def add_marker_legend(self, specs, residual_panel=None, props=None):

        """
        THis function ...
        :param specs:
        :param residual_panel:
        :param props:
        :return:
        """

        # Get the plot
        if residual_panel is not None: plot = self.residual_plots[residual_panel]
        else: plot = self.main_plot

        # Create legend
        legend = self.create_marker_legend(specs, props=props, plot=plot)

        # Add to the plot window
        plot.add_artist(legend)

    # -----------------------------------------------------------------

    def create_line_legend(self, specs, props=None, plot=None):

        """
        This function ...
        :param specs:
        :param props:
        :param props:
        :param plot:
        :return:
        """

        lines = []
        for label in specs:
            color = specs[label]
            line = self.create_line_patch(label, color)
            lines.append(line)

        # Return the legend
        if props is None: props = {}
        if plot is None: plot = self.main_plot
        return plot.create_legend(lines, specs.keys(), **props)

    # -----------------------------------------------------------------

    def add_line_legend(self, specs, residual_panel=None, props=None):

        """
        This function ...
        :param specs:
        :param residual_panel:
        :param props:
        :return:
        """

        # Get the plot
        if residual_panel is not None: plot = self.residual_plots[residual_panel]
        else: plot = self.main_plot

        # Create legend
        legend = self.create_line_legend(specs, props=props, plot=plot)

        # Add to the plot window
        plot.add_artist(legend)

    # -----------------------------------------------------------------

    def create_rectangle_legend(self, specs, props=None, plot=None):

        """
        This function ...
        :param specs:
        :param props:
        :param plot:
        :return:
        """

        rectangles = []
        for label in specs:
            color = specs[label]
            rectangle = self.create_rectangle_patch(label, color)
            rectangles.append(rectangle)

        # Return the legend
        if props is None: props = {}
        if plot is None: plot = self.main_plot
        return plot.create_legend(rectangles, specs.keys(), **props)

    # -----------------------------------------------------------------

    def add_rectangle_legend(self, specs, residual_panel=None, props=None):

        """
        This function ...
        :param specs:
        :param residual_panel:
        :param props:
        :return:
        """

        # Get the plot
        if residual_panel is not None: plot = self.residual_plots[residual_panel]
        else: plot = self.main_plot

        # Create legend
        legend = self.create_rectangle_legend(specs, props=props, plot=plot)

        # Add to the plot window
        plot.add_artist(legend)

    # -----------------------------------------------------------------

    @property
    def do_instruments_legend(self):
        return self.config.legends.instruments

    # -----------------------------------------------------------------

    @property
    def do_observations_legend(self):
        return self.config.legends.observations and self.has_multiple_observations

    # -----------------------------------------------------------------

    @property
    def do_models_legend(self):
        return self.config.legends.models and self.has_multiple_models

    # -----------------------------------------------------------------

    @property
    def do_residuals_legends(self):
        return self.config.legends.residuals and self.has_residual_plots

    # -----------------------------------------------------------------

    def plot_observations(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting the observed SEDs ...")

        # Markers for the unique labels
        self.observation_markers = filled_markers[:len(self.unique_observation_point_labels)]

        # Make iterable from distinct colors
        different_colors = iter(dark_pretty_colors)

        # Loop over the different observed SEDs
        for label in self.observations:

            # Get the observed SED
            observation = self.observations[label]

            # Get options
            only_residuals = self.observation_options[label].only_residuals
            #as_reference = self.observation_options[label].as_reference
            observation_color = self.observation_options[label].color

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

            # Create color range
            if observation_color is None:
                if self.one_observation:

                    # Determine color map class
                    colormap = plt.get_cmap(self.config.single_observation_cmap)

                    # Set colors
                    observation_color = None
                    observation_colors = colormap(np.linspace(0, 1, len(wavelengths)))
                    self.observation_colors[label] = observation_colors
                else:

                    # Set color
                    observation_color = color_hex[next(different_colors)]
                    observation_colors = None
                    self.observation_colors[label] = observation_color

            # Set user-defined color
            else:
                observation_colors = None
                self.observation_colors[label] = observation_color

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
                marker = self.observation_markers[self.unique_observation_point_labels.index(labels[k])]

                # Set zero error if other observations also have errors (BECAUSE CALLING MATPLOTLIB'S ERRORBAR AND AFTER THAT PLOT DOESN'T WORK APPARENTLY!!)
                error = errors[k]
                if self.has_any_observation_errors and error is None: error = ErrorBar.zero()

                # Add SED
                self._observations[label].wavelengths.append(wavelengths[k])
                self._observations[label].fluxes.append(fluxes[k])
                self._observations[label].errors.append(errors[k]) # None if no error
                self._observations[label].instruments.append(instrument)
                self._observations[label].bands.append(band)
                self._observations[label].labels.append(labels[k])
                self._observations[label].descriptions.append(descriptions[k])
                self._observations[label].markers.append(marker)

                # Plot the flux data point on the main axis with the specified marker and color
                # UNLESS ONLY_RESIDUALS FLAG IS ENABLED
                if not only_residuals:

                    # get color
                    if observation_colors is not None: color = observation_colors[k]
                    else: color = observation_color

                    # Plot the data point
                    patch, new_label = self.plot_wavelength(self.main_plot, labels[k], self.used_labels, wavelengths[k],
                                                            fluxes[k], error, marker, color, return_new=True)

                    # Add the patch and label
                    if patch is not None and new_label:

                        # Create black and white patch
                        patch = self.create_marker_patch(labels[k], marker, "white")

                        # Add the patch and label
                        self.instruments_legend_patches.append(patch)
                        self.instruments_legend_labels.append(labels[k])

            # Create rectangle for this observation
            if self.do_observations_legend:
                if not (only_residuals and self.config.only_residuals_legend):
                    rectangle = self.create_rectangle_patch(label, observation_color)
                    self.observations_legend_patches.append(rectangle)
                    self.observations_legend_labels.append(label.replace("_", "\_"))

    # -----------------------------------------------------------------

    @lazyproperty
    def model_labels_only_residuals(self):
        #return [label for label in self.model_labels if self.model_options[label].only_residuals]
        labels = []
        for label in self.model_labels:
            if self.model_options[label].only_residuals:
                if not self.model_options[label].residuals: raise ValueError("Cannot enable 'only_residuals' but disable plotting residuals")
                labels.append(label)
        return labels

    # -----------------------------------------------------------------

    @lazyproperty
    def nmodels_only_residuals(self):
        return len(self.model_labels_only_residuals)

    # -----------------------------------------------------------------

    @property
    def has_models_only_residuals(self):
        return self.nmodels_only_residuals > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def model_linestyle_iterator(self):
        return iter(line_styles)

    # -----------------------------------------------------------------

    def _get_next_model_linecolor_residuals(self):

        # If there are not enough linestyles, different colors
        if self.nmodels_for_residuals > len(line_styles): return sequences.find_first_not_in(["black"] + self.model_linecolors, self.model_colors.values())

        # If there are enough linestyles, only black color (line styles will be different)
        else: return "black"

    # -----------------------------------------------------------------

    def _get_next_model_linestyle_residuals(self):

        # If there are not enough linestyles, just normal lines
        if self.nmodels_for_residuals > len(line_styles): return "-"

        # If there are enough linestyles, return the next
        else: return self.model_linestyle_iterator.next()

    # -----------------------------------------------------------------

    def _get_next_model_linestyle(self):

        # line_styles_models = line_styles if self.nmodels_not_ghost <= len(line_styles) else ["-"] * len(line_colors_models)
        if self.nmodels_not_ghost > len(line_styles): return "-"
        else: return sequences.find_first_not_in(line_styles, self.model_styles.values()) # don't return already used styles (because specified by user for specific models)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_linecolors(self):
        # Make a nice sequence of colors, the most used ones first
        line_colors_models = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        for color in dark_pretty_colors:
            if color not in line_colors_models: line_colors_models.append(color)
        return line_colors_models

    # -----------------------------------------------------------------

    def _get_next_model_linecolor(self):

        # No: because if many models the linestyles are already the same
        #if self.model_labels_not_ghost_not_residuals > len(self.model_linecolors): return "black"
        #else: return sequences.find_first_not_in(self.model_linecolors, self.model_colors.values()) # don't return already used colors (because specified by user for specific models)
        return sequences.find_first_not_in(self.model_linecolors, self.model_colors.values())

    # -----------------------------------------------------------------

    def plot_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model SEDs ...")

        # Loop over the model SEDs
        for model_label in self.models:

            # Debugging
            log.debug("Plotting the '" + model_label + "' SED ...")

            # sed, plot_residuals, ghost = self.models[model_label]
            sed = self.models[model_label]
            plot_residuals = self.model_options[model_label].residuals
            only_residuals = self.model_options[model_label].only_residuals
            ghost = self.model_options[model_label].ghost
            above = self.model_options[model_label].above
            fill = self.model_options[model_label].fill
            model_color = self.model_options[model_label].color
            model_linestyle = self.model_options[model_label].linestyle
            #print(model_label, plot_residuals)

            # Get fluxes, wavelengths and errors
            # print(self.config.unit)
            fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info, asarray=True)
            wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False, asarray=True)
            errors = sed.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info) if sed.has_errors else None

            # Above another
            if above is not None:

                # Get wavelengths and fluxes of other SED
                above_wavelengths = self.models[above].wavelengths(unit=self.config.wavelength_unit, add_unit=False, asarray=True)
                above_fluxes = self.models[above].photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info, asarray=True)

                # Interpolate to same grid
                above_interpolation = interp1d(above_wavelengths, above_fluxes, kind='cubic')
                x_new = wavelengths
                below_bounds = x_new < above_interpolation.x[0]
                above_bounds = x_new > above_interpolation.x[-1]
                valid = np.logical_not(below_bounds) * np.logical_not(above_bounds)
                wavelengths = wavelengths[valid]
                fluxes = fluxes[valid]

                # Determine residuals
                above_fluxes = above_interpolation(wavelengths)
                fluxes = fluxes + above_fluxes

            # Not above another
            else: above_wavelengths = above_fluxes = None

            # Determine linecolor
            if model_color is not None: linecolor = model_color
            elif ghost: linecolor = "lightgrey"
            elif plot_residuals: linecolor = self._get_next_model_linecolor_residuals() #linecolor = "black"
            else: linecolor = self._get_next_model_linecolor()

            #print(model_label, plot_residuals)

            # Determine linestyle
            if model_linestyle is not None: linestyle = model_linestyle
            elif ghost: linestyle = "-"
            elif plot_residuals: linestyle = self._get_next_model_linestyle_residuals()
            else: linestyle = self._get_next_model_linestyle()

            # Set color and linestyle
            self.model_colors[model_label] = linecolor
            self.model_styles[model_label] = linestyle

            # Determine actual label
            if above is not None:
                above_name = self.model_options[model_label].above_name
                if above_name is not None: actual_label = above_name
                elif not fill: actual_label = model_label.replace("_", "\_")
                else: actual_label = None
            else: actual_label = model_label.replace("_", "\_")

            # Plot the model SED as a line (with errors if present)
            if not only_residuals: self.draw_model(self.main_plot, wavelengths, fluxes, linestyle, label=actual_label, errors=errors, linecolor=linecolor, adjust_extrema=not self.has_observations)

            # If residuals still have to be plotted
            #if plot_residuals:
                #self.models_residuals[model_label] = (fluxes, wavelengths, errors)

            # Add model (with fluxes added to another if 'above')
            self._models[model_label].wavelengths = wavelengths
            self._models[model_label].fluxes = fluxes
            self._models[model_label].errors = errors
            self._models[model_label].above_wavelengths = above_wavelengths
            self._models[model_label].above_fluxes = above_fluxes

            # Fill
            if above is not None and fill:

                # Get fill data
                fill_wavelengths = wavelengths
                fill_lower = np.log10(above_fluxes)
                fill_upper = np.log10(fluxes)
                fill_log_wavelengths = np.log10(fill_wavelengths)

                #print("fill")
                #print("wavelengths", fill_wavelengths)
                #print("lower", fill_lower)
                #print("upper", fill_upper)

                # Interpolate nans
                #fill_log_wavelengths = np.log10(fill_wavelengths)

                #nans_lower = np.isnan(fill_lower)
                #nans_upper = np.isnan(fill_upper)
                #nans = nans_lower + nans_upper
                #valid = np.isfinite(fill_lower) * np.isfinite(fill_upper)
                #invalid = np.logical_not(valid)
                #valid_fill_wavelengths = fill_wavelengths[valid]
                #valid_fill_lower = fill_lower[valid]
                #valid_fill_upper = fill_upper[valid]

                fill_lower = interpolate_and_extrapolate_nans(fill_log_wavelengths, fill_lower)
                fill_upper = interpolate_and_extrapolate_nans(fill_log_wavelengths, fill_upper)

                # Set fill data
                self.fills[model_label] = (fill_wavelengths, fill_lower, fill_upper)

                # Fill
                self.main_plot.axes.fill_between(fill_wavelengths, fill_lower, fill_upper, facecolor=self.model_colors[model_label], alpha=0.5, label=model_label.replace("_", "\_"))

            # Create for legend
            if not ghost and not (only_residuals and self.config.only_residuals_legend):

                # Create line patch
                patch = self.create_line_patch(model_label, linecolor, style=linestyle)

                # Add the patch and label
                self.models_legend_patches.append(patch)
                self.models_legend_labels.append(model_label.replace("_", "\_"))

    # -----------------------------------------------------------------

    @property
    def residual_reference(self):
        if self.no_models: return "observations"
        elif self.no_observations: return "models"
        else: return self.config.residual_reference

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals ...")

        # Observations as reference: plot at 0.0 (all in one panel)
        if self.residual_reference == observations_reference:

            # Plot observations each on their own residual axes
            self._plot_obs_residual_references()

            # Plot model to obs residuals
            if self.has_models: self._plot_model_to_obs_residuals()

            # Residuals between obs and obs
            # LOOP OVER OBSERVATIONS AGAIN IF WE HAVE TO PLOT THEIR RESIDUALS AGAINST THE REFERENCE (IF OBSERVATION IS REFERENCE)
            if self.config.observations_residuals: self._plot_obs_to_obs_residuals()

        # Models as reference: plot points at the relative difference with the models (one panel for each model)
        elif self.residual_reference == models_reference:

            # Plot models
            self._plot_obs_to_model_residuals()

            # Plot model to model
            if self.config.models_residuals: self._plot_model_to_model_residuals()

        # Invalid reference
        else: raise RuntimeError("Invalid value for 'residual_reference'")

    # -----------------------------------------------------------------

    def _plot_obs_residual_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the observation references on the residual panels ...")

        # Only observations?
        if self.no_models:

            # Determine the reference SED
            reference_sed_label = self.most_npoints_observation_label
            #reference_sed = self.observations[reference_sed_label]

            # Loop over the points
            obs = self._observations[reference_sed_label]
            for wavelength, flux, error, marker in zip(obs.wavelengths, obs.fluxes, obs.errors, obs.markers):

                if error is not None:
                    value = 0.0
                    error = error / flux * 100.
                else: value = error = None

                if value is not None and error is not None:

                    yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)
                    self.residual_plot.errorbar(wavelength, value, yerr=yerr, fmt=marker, markersize=7,
                                                color=self.observation_colors[reference_sed_label], markeredgecolor='black', ecolor=self.observation_colors[reference_sed_label],
                                                capthick=2, lolims=lolims, uplims=uplims, capsize=2)

        # THere are models
        else:

            # Loop over the observed SEDs
            observation_index = observation_reference_index = 0
            for label in self.observations:

                # Get options
                as_reference = self.observation_options[label].as_reference
                #observation_color = self.observation_options[label].color

                # Get color
                observation_color = self.observation_colors[label]

                # print(label, observation_index, wavelengths[k])
                if as_reference or self.observation_options[label].join_residuals:

                    # Get the residual plot
                    # print(observation_reference_index, len(self.residual_plots))
                    if as_reference: residual_plot = self.residual_plots[observation_reference_index]
                    else:
                        join_with_label = self.observation_options[label].join_residuals
                        residual_plot = self.residual_plots[self.observation_labels_as_reference.index(join_with_label)]

                    # Plot measurement points on residual plot
                    value = 0.0

                    # Loop over the points
                    obs = self._observations[label]
                    for wavelength, flux, error, marker in zip(obs.wavelengths, obs.fluxes, obs.errors, obs.markers):

                        # Error is defined
                        if error is not None:

                            # Process and plot errorbar
                            error = error / flux * 100.

                            # yerr, lolims, uplims = process_errorbar(error)
                            yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)

                            residual_plot.errorbar(wavelength, value, yerr=yerr, fmt=marker, markersize=7,
                                                   color=observation_color, markeredgecolor='black',
                                                   ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims,
                                                   capsize=2)

                        else:
                            yerr = [[0], [0]]
                            residual_plot.errorbar(wavelength, value, yerr=yerr, fmt=marker, markersize=7,
                                                   color=observation_color, markeredgecolor='black',
                                                   ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims,
                                                   capsize=2)

                # Increment the observation index
                observation_index += 1
                if as_reference: observation_reference_index += 1

    # -----------------------------------------------------------------

    def _plot_obs_to_model_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting observation residuals compared to model SEDs ...")

        from astropy.units import Unit

        residual_plot_index = 0

        # Loop over the observations
        for label in self.observations:

            # Get color
            observation_color = self.observation_colors[label]

            # Loop over the points
            obs = self._observations[label]
            for wavelength, flux, error, marker in zip(obs.wavelengths, obs.fluxes, obs.errors, obs.markers):

                # Loop over the models
                for model_label in self.models:

                    # Get entry
                    # sed, plot_residuals, ghost = self.models[model_label]
                    sed = self.models[model_label]
                    plot_residuals = self.model_options[model_label].residuals
                    ghost = self.model_options[model_label].ghost

                    # Don't plot residuals for this model
                    if not plot_residuals: continue

                    # model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                    # sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                    actual_wavelength = wavelength * Unit(self.config.wavelength_unit)

                    # Wavelength out of bounds of model SED
                    if actual_wavelength < sed.min_wavelength or actual_wavelength > sed.max_wavelength:
                        residual_plot_index += 1
                        continue

                    # print(model_label, sed.unit, sed.unit.physical_type)
                    model_value = sed.photometry_at(actual_wavelength, self.config.unit, add_unit=False,
                                                    conversion_info=self.conversion_info,
                                                    interpolate=self.config.interpolate_models_for_residuals)
                    try: rel_residual = (flux - model_value) / model_value * 100.  # NORMALIZE TO MODEL VALUE
                    except ZeroDivisionError:
                        residual_plot_index += 1
                        continue

                    # Get the residual plot for this model
                    residual_plot = self.residual_plots[residual_plot_index]

                    # value = 0.0
                    if error is None: errorbar = ErrorBar.zero()
                    else: errorbar = error

                    # NORMALIZE TO MODEL VALUE
                    error = errorbar / model_value * 100.

                    # yerr, lolims, uplims = process_errorbar(error)
                    yerr, lolims, uplims = process_errorbar_for_value(flux, error, lolim_abs_value=-75,
                                                                      uplim_value=75)

                    # Plot on residual axes
                    residual_plot.errorbar(wavelength, rel_residual, yerr=yerr, fmt=marker,
                                           markersize=7, color=observation_color, markeredgecolor='black',
                                           ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

                    # Next residual plot
                    residual_plot_index += 1

    # -----------------------------------------------------------------

    def _plot_model_to_model_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model residuals compared to other model SEDs ...")

        # Pair of models and no observations
        if self.no_observations and self.nmodels_for_residuals == 2:

            label_a = self.model_labels_for_residuals[0]
            label_b = self.model_labels_for_residuals[1]
            sed_a = self.models[label_a]
            sed_b = self.models[label_b]

            # Get wavelengths
            wavelengths_a = sed_a.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            wavelengths_b = sed_b.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

            # Get fluxes
            fluxes_a = sed_a.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
            fluxes_b = sed_b.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

            # Determine min and max wavelength
            min_wavelength = max(min(wavelengths_a), min(wavelengths_b))
            max_wavelength = min(max(wavelengths_a), max(wavelengths_b))

            # Filter
            wavelengths_fluxes_a = [(wavelength, flux) for wavelength, flux in zip(wavelengths_a, fluxes_a) if
                                    min_wavelength < wavelength < max_wavelength]
            wavelengths_fluxes_b = [(wavelength, flux) for wavelength, flux in zip(wavelengths_b, fluxes_b) if
                                    min_wavelength < wavelength < max_wavelength]

            # Wavelengths and fluxes again
            wavelengths_a = [wf[0] for wf in wavelengths_fluxes_a]
            fluxes_a = [wf[1] for wf in wavelengths_fluxes_a]
            wavelengths_b = [wf[0] for wf in wavelengths_fluxes_b]
            fluxes_b = [wf[1] for wf in wavelengths_fluxes_b]

            # Interpolate to same grid
            f2 = interp1d(wavelengths_b, fluxes_b, kind='cubic')

            # Determine residuals
            fluxes_b_wavelengths_a = f2(wavelengths_a)
            residuals = (fluxes_a - fluxes_b_wavelengths_a) / fluxes_a * 100.
            residuals_a = 0.5 * residuals
            residuals_b = - 0.5 * residuals

            # Plot both
            self.residual_plot.plot(wavelengths_a, residuals_a, linestyle=self.model_styles[label_a],
                                    color=self.model_colors[label_a], label=label_a)
            self.residual_plot.plot(wavelengths_a, residuals_b, linestyle=self.model_styles[label_b],
                                    color=self.model_colors[label_b], label=label_b)

        else:

            # Get reference
            reference_model_label = self.model_labels_for_residuals[0]
            reference_model_sed = self.models[reference_model_label]
            reference_wavelengths = reference_model_sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            reference_fluxes = reference_model_sed.photometry(unit=self.config.unit, add_unit=False,
                                                              conversion_info=self.conversion_info)

            # Plot reference
            residuals = [0.0] * len(reference_wavelengths)
            self.residual_plot.plot(reference_wavelengths, residuals, linestyle=self.model_styles[reference_model_label],
                                    color=self.model_colors[reference_model_label], label=reference_model_label)

            # Loop over the models to be plotted
            for model_label in self.model_labels_for_residuals:

                if model_label == reference_model_label: continue

                # Get color and style
                linecolor = self.model_colors[model_label]
                linestyle = self.model_styles[model_label]

                # Get the sed
                sed = self.models[model_label]

                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                min_sed_wavelength = min(sed_wavelengths)
                max_sed_wavelength = max(sed_wavelengths)

                wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(reference_wavelengths, reference_fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]

                fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                # Plot
                self.residual_plot.plot(wavelengths_residuals, residuals, linestyle=linestyle, color=linecolor, label=model_label.replace("_", "\_"))

    # -----------------------------------------------------------------

    def _plot_model_to_obs_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model residuals compared to observed SEDs ...")

        # Loop over the different observed SEDs
        observation_reference_index = 0
        for label in self.observations:

            # Get options
            #only_residuals = self.observation_options[label].only_residuals
            as_reference = self.observation_options[label].as_reference

            # Observation not as reference?
            if not as_reference: continue

            # JOINED WITH OTHER
            if label in self.observation_labels_with_joined:

                joined_observation = self.joined_observations[label]
                #print(label + " joined")
                # print(joined_observation)
                # print(joined_observation.wavelength_unit, joined_observation.unit)

                wavelengths = joined_observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                fluxes = joined_observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                bands = joined_observation.bands()

            # Not joined
            else:

                obs = self._observations[label]
                # for wavelength, flux, error, marker in zip(obs.wavelengths, obs.fluxes, obs.errors, obs.markers):
                wavelengths = obs.wavelengths
                fluxes = obs.fluxes
                bands = obs.bands

            # Get the residual plot
            residual_plot = self.residual_plots[observation_reference_index]
            observation_reference_index += 1

            #print(wavelengths)
            #print(fluxes)

            # PLOT RESIDUALS WITH MODELS
            for model_label in self.models:

                # Get entry
                # sed, plot_residuals, ghost = self.models[model_label]
                sed = self.models[model_label]
                plot_residuals = self.model_options[model_label].residuals
                ghost = self.model_options[model_label].ghost
                residual_color = self.model_options[model_label].residual_color

                if not plot_residuals: continue

                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info, asarray=True)
                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False, asarray=True)

                # FIRST SMOOTH THE OBSERVATION, THEN GET THE 'CONTINUOUS' RESIDUAL LINE
                if self.config.smooth_residuals:

                    logwavelengths = np.log10(wavelengths)
                    logfluxes = np.log10(fluxes)

                    # gradient = np.gradient(logwavelengths)
                    diffs = np.ediff1d(logwavelengths)
                    # print(len(gradient), len(wavelengths))
                    # print(diffs)
                    # print(len(diffs))
                    reldiffs = diffs / logwavelengths[:-1]
                    # for grad, band in zip(gradient, bands): print(band, grad)
                    ignore_bands = []
                    for first_band, first_wavelength, second_band, second_wavelength, reld in zip(bands[:-1], wavelengths[:-1], bands[1:], wavelengths[1:], reldiffs):

                        # print(first_band + " vs " + second_band + ": " + str(reld))
                        if abs(reld) < self.config.smooth_residuals_closeness_limit:
                            log.warning("The " + str(first_band) + " and " + second_band + " wavelengths are very close (" + str(first_wavelength) + " and " + str(second_wavelength) + ", closeness = " + str(abs(reld)) + "), but not necessarily their photometric value: this may give artefacts in the residual curve")
                            if self.config.smooth_residuals_ignore_close: ignore_bands.append(second_band)

                    # Ignore bands?
                    if len(ignore_bands) > 0:

                        ignore_indices = [bands.index(band) for band in ignore_bands]
                        ignore_mask = np.zeros(len(bands), dtype=bool)
                        ignore_mask[ignore_indices] = 1
                        pass_mask = np.logical_not(ignore_mask)
                        logwavelengths = logwavelengths[pass_mask]
                        logfluxes = logfluxes[pass_mask]

                    # Create spline (interpolate in log-log space)
                    logsed_wavelengths = np.log10(sed_wavelengths)
                    #logsmooth_model = interp1d(logwavelengths, logfluxes, kind="cubic")
                    #logsmooth_model = interp1d(logwavelengths, logfluxes, fill_value="extrapolate")
                    #logsmooth_fluxes = logsmooth_model(logsed_wavelengths)
                    logsmooth_fluxes = spline(logwavelengths, logfluxes, logsed_wavelengths)
                    smooth_fluxes = 10 ** logsmooth_fluxes

                    min_obs_wavelength = wavelengths[0]
                    max_obs_wavelength = wavelengths[-1]
                    within_observations_mask = (sed_wavelengths > min_obs_wavelength) * ( sed_wavelengths < max_obs_wavelength)
                    # print(within_observations_mask)

                    # Mask
                    wavelengths_residuals = sed_wavelengths[within_observations_mask]
                    smooth_fluxes = smooth_fluxes[within_observations_mask]
                    model_fluxes = model_fluxes[within_observations_mask]

                    # Plot the model SED as a line (with errors if present)
                    if self.config.show_smooth: self.draw_model(self.main_plot, wavelengths_residuals,
                                                                smooth_fluxes, "-", linecolor="lightgrey",
                                                                adjust_extrema=False)

                    # Calculate residuals
                    residuals = - (smooth_fluxes - model_fluxes) / smooth_fluxes * 100.

                # INTERPOLATE THE SED WHERE THE FLUX POINTS ARE TO COMPARE
                else:

                    f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                    min_sed_wavelength = min(sed_wavelengths)
                    max_sed_wavelength = max(sed_wavelengths)

                    wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                    wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                    fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                    residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                # SMOOTH?
                # THIS DOES NOT WORK WELL: RESIDUAL POINTS ARE NOT SMOOTH IN SE: THEY VARY BELOW AND ABOVE the expected:
                # THE INTERPOLATION GIVES MANY ARTEFACTS (over-fitting polynomial behaviour)
                # if self.config.smooth_residuals:
                #    residuals = spline(wavelengths_residuals, residuals, sed_wavelengths)
                #    wavelengths_residuals = sed_wavelengths

                # Get line color and style
                linecolor = self.model_colors[model_label] if residual_color is None else residual_color
                linestyle = self.model_styles[model_label]

                # Plot
                if ghost: actual_label = None
                else: actual_label = model_label.replace("_", "\_")
                residual_plot.plot(wavelengths_residuals, residuals, ls=linestyle, color=linecolor, label=actual_label)

    # -----------------------------------------------------------------

    def _plot_obs_to_obs_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the observation residuals compared to other observed SEDs ...")

        # Only observations: one panel for all observation residuals
        if self.no_models:

            # Determine the reference SED
            reference_sed_label = self.most_npoints_observation_label
            reference_sed = self.observations[reference_sed_label]

            # Loop over the other observations
            for label in self.observations:

                # Get the observed SED and its values
                #observation = self.observations[label]
                if label == reference_sed_label: continue
                obs = self._observations[label]

                # Get values
                # Get wavelengths, fluxes, instruments, bands, errors
                #wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                #fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                #instruments = observation.instruments()
                #bands = observation.bands()
                #errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

                # Loop over the points
                #for k in range(len(wavelengths)):
                for instrument, band, wavelength, flux, error, marker in zip(obs.instruments, obs.bands, obs.wavelengths, obs.fluxes, obs.errors, obs.markers):

                    # Find reference flux
                    reference_flux = find_reference_flux(instrument, band, wavelength, reference_sed.instruments(),
                                                         reference_sed.bands(), reference_sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False),
                                                         reference_sed.photometry(unit=self.config.unit, add_unit=False))

                    if reference_flux is None: continue
                    #else:
                    if reference_flux == 0: continue #value = error = None
                    #else:
                    if error is not None:
                        value = (error - reference_flux) / reference_flux * 100.
                        error = error / reference_flux * 100.0
                    # else: value = error = None
                    elif flux is not None:
                        value = (flux - reference_flux) / reference_flux * 100.
                        error = ErrorBar.zero()
                    else: continue #value = error = None

                    #if value is not None and error is not None:
                    # yerr, lolims, uplims = process_errorbar(error)
                    yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)
                    self.residual_plot.errorbar(wavelength, value, yerr=yerr, fmt=marker, markersize=7,
                                                color=self.observation_colors[label], markeredgecolor='black', ecolor=self.observation_colors[label],
                                                capthick=2, lolims=lolims, uplims=uplims, capsize=2)

        # There are models: each reference observation gets its own panel
        else:

            # Loop over the observations that are used as reference (that have a residual panel for them)
            observation_reference_index = 0
            for label in self.observations:

                # Get options
                # only_residuals = self.observation_options[label].only_residuals
                as_reference = self.observation_options[label].as_reference
                # print(label, as_reference, observation_reference_index)
                if not as_reference: continue

                # Get the observed SED and its values
                observation = self.observations[label]

                # JOINED WITH OTHER
                if label in self.observation_labels_with_joined:

                    joined_observation = self.joined_observations[label]
                    reference_wavelengths = joined_observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                    reference_fluxes = joined_observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                    reference_instruments = joined_observation.instruments()
                    reference_bands = joined_observation.bands()

                else:

                    reference_wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                    reference_fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                    reference_instruments = observation.instruments()
                    reference_bands = observation.bands()

                # Get the residual plot
                residual_plot = self.residual_plots[observation_reference_index]

                # Increment
                observation_reference_index += 1

                # Loop over the other observations
                for other_label in self.observations:
                    if other_label == label: continue

                    # Get the color defined for this observation
                    points_color = self.observation_colors[other_label]

                    # Get the other observed SED and its values
                    other_observation = self.observations[other_label]
                    wavelengths = other_observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                    fluxes = other_observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                    instruments = other_observation.instruments()
                    bands = other_observation.bands()
                    errors = other_observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

                    # Get labels and descriptions for the filters
                    labels, descriptions = get_labels_and_descriptions(instruments, bands)

                    # Loop over the points of the other observation
                    # Loop over the wavelengths
                    for k in range(len(wavelengths)):

                        # Get marker
                        marker = self.observation_markers[self.unique_observation_point_labels.index(labels[k])]

                        # Find reference flux
                        reference_flux = find_reference_flux(instruments[k], bands[k], wavelengths[k], reference_instruments, reference_bands, reference_wavelengths, reference_fluxes)
                        # print(fluxes[k], reference_flux)

                        if reference_flux is None: continue
                        if reference_flux == 0: continue

                        if errors[k] is not None:

                            value = (fluxes[k] - reference_flux) / reference_flux * 100.
                            error = errors[k] / reference_flux * 100.0

                        # else: value = error = None
                        elif fluxes[k] is not None:

                            value = (fluxes[k] - reference_flux) / reference_flux * 100.
                            error = ErrorBar.zero()

                        else: value = error = None

                        # Plot residual point
                        yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)
                        residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7,
                                               color=points_color, markeredgecolor='black', ecolor=points_color,
                                               capthick=2, lolims=lolims, uplims=uplims, capsize=2)

    # -----------------------------------------------------------------

    def create_instruments_legend(self):

        """
        This function ...
        :return:
        """

        # Create
        instruments_legend = self.main_plot.create_legend(self.instruments_legend_patches, self.instruments_legend_labels, **self.instruments_legend_properties)

        # Add the legend
        self.legends.append(instruments_legend)

    # -----------------------------------------------------------------

    def create_observations_legend(self):

        """
        This function ...
        :return:
        """

        # OBSERVATIONS
        observations_legend = self.main_plot.create_legend(self.observations_legend_patches,
                                                           self.observations_legend_labels,
                                                           **self.observations_legend_properties)
        # Extra legend: the different observations
        # observations_legend = self.main_plot.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Add the legend
        self.legends.append(observations_legend)

    # -----------------------------------------------------------------

    def create_models_legend(self):

        """
        This function ...
        :return:
        """

        # MODELS
        models_legend = self.main_plot.create_legend(self.models_legend_patches, self.models_legend_labels,
                                                     **self.models_legend_properties)

        # Add the legend
        self.legends.append(models_legend)

    # -----------------------------------------------------------------

    def create_residuals_legends(self):

        """
        This function ...
        :return:
        """

        # Initialize legends
        self.legends_residuals = [[] for _ in self.residual_plots]

        # Add legend for the only_residuals observations?
        if self.config.only_residuals_legend and self.has_observations_only_residuals:

            legend_patches = []
            legend_labels = []

            # Loop over the observations for only residuals
            for label in self.observation_labels_only_residuals:

                # Create rectangle
                rectangle = self.create_rectangle_patch(label, self.observation_colors[label])

                # Add
                legend_patches.append(rectangle)
                legend_labels.append(label.replace("_", "\_"))

            # Create legend for each residual axis
            for index, residual_plot in enumerate(self.residual_plots):
                # residual_plot = self.residual_plots[self.observation_labels_as_reference.index()]

                # Create legend
                legend = residual_plot.create_legend(legend_patches, legend_labels,
                                                     **self.observations_residuals_legend_properties)

                # Add legend
                self.legends_residuals[index].append(legend)

        # Add legend for the only_residuals models?
        if self.config.only_residuals_legend and self.has_models_only_residuals:

            legend_patches = []
            legend_labels = []

            # Loop over the models for only residuals
            for label in self.model_labels_only_residuals:
                # print("RES LEGEND MODEL", label)

                residual_color = self.model_options[label].residual_color
                linecolor = self.model_colors[label] if residual_color is None else residual_color

                # Create
                line = self.create_line_patch(label, linecolor, style=self.model_styles[label])

                # Add patch and label
                legend_patches.append(line)
                legend_labels.append(label.replace("_", "\_"))

            # print(legend_patches, legend_labels)

            # Create legend for each residual axis
            for index, residual_plot in enumerate(self.residual_plots):

                # Create legend
                legend = residual_plot.create_legend(legend_patches, legend_labels,
                                                     **self.models_residuals_legend_properties)

                #print(legend)

                # Add legend
                self.legends_residuals[index].append(legend)

    # -----------------------------------------------------------------

    def create_legends(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the legends ...")

        # Instruments
        if self.do_instruments_legend: self.create_instruments_legend()

        # Observations
        if self.do_observations_legend: self.create_observations_legend()

        # Models
        if self.do_models_legend: self.create_models_legend()

        # Residuals
        if self.do_residuals_legends: self.create_residuals_legends()

    # -----------------------------------------------------------------

    @lazyproperty
    def adjust_wavelength_minmax_observations(self):
        if self.config.minmax_wavelength_reference is None: return True # both observations and models
        else: return self.config.minmax_wavelength_reference == observations_reference # only observations

    # -----------------------------------------------------------------

    @lazyproperty
    def adjust_wavelength_minmax_models(self):
        if self.config.minmax_wavelength_reference is None: return True  # both observations and models
        else: return self.config.minmax_wavelength_reference == models_reference  # only models

    # -----------------------------------------------------------------

    @lazyproperty
    def adjust_photometry_minmax_observations(self):
        if self.config.minmax_photometry_reference is None: return True # both observations and models
        else: return self.config.minmax_photometry_reference == observations_reference # only observations

    # -----------------------------------------------------------------

    @lazyproperty
    def adjust_photometry_minmax_models(self):
        if self.config.minmax_photometry_reference is None: return True # both observations and models
        else: return self.config.minmax_photometry_reference == models_reference # only models

    # -----------------------------------------------------------------

    def adjust_photometry_minmax(self, flux, error=None):

        """
        This function ...
        :param flux:
        :param error:
        :return:
        """

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

    # -----------------------------------------------------------------

    def adjust_wavelength_minmax(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        if self._min_wavelength is None or wavelength < self._min_wavelength: self._min_wavelength = wavelength
        if self._max_wavelength is None or wavelength > self._max_wavelength: self._max_wavelength = wavelength

    # -----------------------------------------------------------------

    def plot_wavelength(self, axis, label, used_labels, wavelength, flux, error, marker, color, return_new=False,
                        markersize=7, markeredgecolor="black"):

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
        :param markersize:
        :param markeredgecolor:
        :return:
        """

        # Keep track of minimum and maximum flux
        if self.adjust_photometry_minmax_observations: self.adjust_photometry_minmax(flux, error=error)

        # Keep track of the minimal and maximal wavelength
        if self.adjust_wavelength_minmax_observations: self.adjust_wavelength_minmax(wavelength)

        # Check if a data point of this instrument has already been plotted
        new_label = label not in used_labels
        if new_label:

            used_labels.append(label)

            if error is not None:

                #print("a", label, new_label)

                #lower_flux = flux + error.lower
                #upper_flux = flux + error.upper

                #lolims = False
                #uplims = False

                #if lower_flux <= 0:
                #    flux_lower_flux = float("-inf")
                #    uplims = [[True]]
                #else: flux_lower_flux = flux / lower_flux

                #if error.lower == numbers.min_inf: uplims = [[True]]
                #if error.upper == numbers.inf: lolims = [[True]]

                #flux_upper_flux = flux / upper_flux
                #yerr = np.array([[np.fabs(np.log10(flux_lower_flux)), np.fabs(np.log10(flux_upper_flux))]]).T

                logflux = np.log10(flux)
                yerr, lolims, uplims = process_errorbar_for_value(flux, error, min_allowed_value=0, return_log=True, lolim_abs_value=0)
                #print(label, yerr)
                if yerr[0][0] == numbers.inf: yerr[0][0] = 20 # logflux - np.log10(self._min_flux) LARGE ENOUGH SO THAT IT GOES COMPLETELY DOWN
                if yerr[1][0] == numbers.inf: yerr[1][0] = np.log10(self._max_flux) - logflux
                patch = axis.errorbar(wavelength, logflux, yerr=yerr, fmt=marker, markersize=markersize, color=color, markeredgecolor=markeredgecolor, ecolor=color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

            else:
                #print("b", label, new_label)
                patch = axis.plot(wavelength, np.log10(flux), marker=marker, markersize=markersize, color=color, markeredgecolor=markeredgecolor) #markerfacecolor=color) #label=label)

        # A data point of this instrument has already been plotted
        else:

            # There is an error bar
            if error is not None:

                #print("c", label, new_label)

                lolims = False
                uplims = False

                #yerr = np.array([[np.fabs(np.log10(flux)-np.log10(flux + error.lower)), np.fabs(np.log10(flux) - np.log10(flux + error.upper))]]).T
                #yerr, lolims, uplims = process_errorbar(ErrorBar(np.log10(flux)-np.log10(flux + error.lower), np.log10(flux) - np.log10(flux + error.upper)))

                #if error.lower == numbers.min_inf: uplims = [[True]]
                #if error.upper == numbers.inf: lolims = [[True]]

                #print("c", wavelength, flux)
                logflux = np.log10(flux)
                yerr, lolims, uplims = process_errorbar_for_value(flux, error, min_allowed_value=0, return_log=True, lolim_abs_value=0)
                #print(label, yerr)
                if yerr[0][0] == numbers.inf: yerr[0][0] = 20 #logflux - np.log10(self._min_flux)  LARGE ENOUGH SO THAT IT GOES COMPLETELY DOWN
                if yerr[1][0] == numbers.inf: yerr[1][0] = np.log10(self._max_flux) - logflux
                patch = axis.errorbar(wavelength, logflux, yerr=yerr, fmt=marker, markersize=markersize, color=color, markeredgecolor=markeredgecolor, ecolor=color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

            #else: axis.plot(wavelength, np.log10(flux), fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2)
            else: patch = axis.plot(wavelength, np.log10(flux), marker=marker, markersize=markersize, color=color, markeredgecolor=markeredgecolor)

        # Return the patch if requested
        if return_new: return patch, new_label
        else: return patch

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

        #print(instruments, bands, wavelengths, fluxes, errors, colors)

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
            #print(color)

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            #print(labels[k], fluxes[k], errors[k])

            # Plot on the main axis with the specified marker and color
            #print("fluxes", fluxes)
            #print("errors", errors)
            # axis, label, used_labels, wavelength, flux, error, marker, color, return_patch=False
            self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

    # -----------------------------------------------------------------

    def adjust_photometry_minmax_fluxes(self, fluxes, errors=None):

        """
        This function ...
        :param fluxes:
        :param errors:
        :return:
        """

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

    # -----------------------------------------------------------------

    def adjust_wavelength_minmax_wavelengths(self, wavelengths):

        """
        This function ...
        :param wavelengths:
        :return:
        """

        # Keep track of the minimal and maximal wavelength
        min_wavelength_model = min(wavelengths)
        max_wavelength_model = max(wavelengths)

        if min_wavelength_model < self._min_wavelength or self._min_wavelength is None: self._min_wavelength = min_wavelength_model
        if max_wavelength_model > self._max_wavelength or self._max_wavelength is None: self._max_wavelength = max_wavelength_model

    # -----------------------------------------------------------------

    def draw_model(self, axis, wavelengths, fluxes, linestyle, label=None, errors=None, linecolor='black',
                   errorcolor='cyan', adjust_extrema=True):

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

            #print(label, fluxes)
            if self.adjust_photometry_minmax_models: self.adjust_photometry_minmax_fluxes(fluxes, errors=errors)
            if self.adjust_wavelength_minmax_models: self.adjust_wavelength_minmax_wavelengths(wavelengths)

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

            # In logspace
            log_bottom = np.log10(bottom)
            log_top = np.log10(top)

            # Fill between min and max of errorbar
            axis.fill_between(errors_wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor=errorcolor, edgecolor=errorcolor, interpolate=True, alpha=0.5)
            axis.plot([], [], color=errorcolor, linewidth=10, label='spread')

    # -----------------------------------------------------------------

    @lazyproperty
    def instruments_legend_properties(self):

        """
        This function ...
        :return:
        """

        # Initialize
        properties = dict()

        # Set properties
        properties["loc"] = self.config.legends.instruments_location
        properties["numpoints"] = 1
        properties["scatterpoints"] = 4
        properties["ncol"] = self.config.legends.instruments_ncols
        properties["shadow"] = False
        properties["frameon"] = True
        properties["facecolor"] = None
        properties["fontsize"] = self.config.plot.legend_fontsize

        # Return
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def observations_legend_properties(self):

        """
        This function ...
        :return:
        """

        # Initialize
        properties = dict()

        # Set properties
        properties["loc"] = self.config.legends.observations_location
        properties["numpoints"] = 1
        properties["scatterpoints"] = 4
        properties["ncol"] = self.config.legends.observations_ncols
        properties["shadow"] = False
        properties["frameon"] = True
        properties["facecolor"] = None
        properties["fontsize"] = self.config.plot.legend_fontsize

        # Return
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def observations_residuals_legend_properties(self):

        """
        This function ...
        :return:
        """

        # Initialize
        properties = dict()

        # Set properties
        properties["loc"] = self.config.legends.observations_residuals_location
        properties["numpoints"] = 1
        properties["scatterpoints"] = 4
        properties["ncol"] = self.config.legends.observations_residuals_ncols
        properties["shadow"] = False
        properties["frameon"] = True
        properties["facecolor"] = None
        properties["fontsize"] = self.config.plot.legend_fontsize

        # Return
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def models_legend_properties(self):

        """
        This function ...
        :return:
        """

        # Initialize
        properties = dict()

        # Set properties
        properties["loc"] = self.config.legends.models_location
        properties["ncol"] = self.config.legends.models_ncols
        properties["shadow"] = False
        properties["frameon"] = True
        properties["facecolor"] = None
        properties["fontsize"] = self.config.plot.legend_fontsize

        # Return
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def models_residuals_legend_properties(self):

        """
        This function ...
        :return:
        """

        # Initialize
        properties = dict()

        # Set properties
        properties["loc"] = self.config.legends.models_residuals_location
        properties["ncol"] = self.config.legends.models_residuals_ncols
        properties["shadow"] = False
        properties["frameon"] = True
        properties["facecolor"] = None
        properties["fontsize"] = self.config.plot.legend_fontsize

        # Return
        return properties

    # -----------------------------------------------------------------

    def finish_plot(self):

        """
        This function ...
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

        #print(self.min_flux, self.max_flux)

        if not had_min_flux: self.min_flux /= factor_y
        if not had_max_flux: self.max_flux *= factor_y
        if not had_min_wavelength: self.min_wavelength /= factor_x
        if not had_max_wavelength: self.max_wavelength *= factor_x

        #print(self.min_flux, self.max_flux)

        # Format residual axes
        for res_axis in self.residual_plots:

            # Set linestyle and limit for axis2
            res_axis.axhline(y=0., color='black', ls='-.')
            res_axis.set_ylim(-95, 95)

            # Set axis label
            res_axis.set_xscale('log')
            res_axis.set_ylabel(r"Res. $[\%]$", fontsize='large')

        # Hide x ticks of all but last residual plot
        if self.has_residual_plots:

            # Hide other tick labels but last residual plot
            self.main_plot.hide_xtick_labels(shared_axis=True)
            for j in range(len(self.residual_plots)-1): self.residual_plots[j].hide_xtick_labels(shared_axis=True)

            # Set xticks for last residual plot
            last_residual_plot_index = len(self.residual_plots)-1
            self.residual_plots[last_residual_plot_index].set_xticks(fontsize=self.config.plot.ticks_fontsize)

        # Set x label of the last residual plot
        #xlabel = r"Wavelength $\lambda\,[\mu \mathrm{m}]$"
        xlabel = r"Wavelength $\lambda\,[\mathrm{" + str(self.config.wavelength_unit) + "}]$".replace("\mathrm{um}", "\mu \mathrm{m}")
        #print(xlabel)
        if len(self.residual_plots) > 0: self.residual_plots[len(self.residual_plots)-1].set_xlabel(xlabel, fontsize='large')
        else: self.main_plot.set_xlabel(xlabel, fontsize='large')

        # Set log x scale
        self.main_plot.set_xscale('log')

        # Set ticks
        #if len(self.residual_plots) == 0: self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        if not self.has_residual_plots: self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        else: self.residual_plots[-1].set_xticks(fontsize=self.config.plot.ticks_fontsize)
        self.main_plot.set_yticks(fontsize=self.config.plot.ticks_fontsize) # minor = False is possible

        # Add axis labels and a legend
        y_label = r"log$_{10} " + self.config.unit.latex_string + r"$"
        #print(y_label)
        self.main_plot.set_ylabel(y_label, fontsize='large')

        # SHOW TICK LABELS
        #print(self.main_plot.xtick_labels, self.main_plot.ytick_labels)
        #for res in self.residual_plots: print(res.xtick_labels, res.ytick_labels)

        # Set position
        # X: only for main plot
        if self.config.xaxis_position == "top": self.main_plot.set_xaxis_position("top")

        # Y: for main and residual plots
        if self.config.yaxis_position == "right":
            for plot in [self.main_plot] + self.residual_plots: plot.set_yaxis_position("right")

        # Set grid
        #self.figure.set_grid(self.config.plot, which="both")
        self.main_plot.set_grid(self.config.plot, which="both")
        for plot in self.residual_plots: plot.set_grid(self.config.plot, which="both")

        # FAILED TEST TO HIDE TICK LABELS ON MAIN BUT SHOW ON LAST RESIDUAL PLOT
        #if self.has_residual_plots:
        #    self.residual_plots[-1].set_xticks(fontsize=self.config.plot.ticks_fontsize)
        #    #self.main_plot.hide_xtick_labels()

        # Set borders
        if not self._external_figure: self.figure.set_borders(self.config.plot)

        # Set flux axis limits
        #print(self.min_flux, self.max_flux)
        plot_min, plot_max = get_plot_flux_limits(self.min_flux, self.max_flux)
        #print(plot_min, plot_max)

        # Set Y limits
        if plot_min is None:
            if plot_max is None: pass
            else: self.main_plot.set_ylim(top=plot_max)
        else:
            if plot_max is None: self.main_plot.set_ylim(bottom=plot_min)
            else: self.main_plot.set_ylim((plot_min, plot_max)) # NORMAL

        # Set wavelength axis limits
        self.main_plot.set_xlim(self.min_wavelength, self.max_wavelength)

        # Add the plots to the figure
        if self.config.library == bokeh and not self._external_figure: self.figure.add_column(self.main_plot, *self.residual_plots)

        # Add legends
        self.add_legends()

        # Add title if requested
        #print(self.title)
        if self.title is not None and not self._external_figure: self.figure.set_title(self.title) #self.figure.figure.suptitle("\n".join(wrap(self.title, 60)))

        # Save or show the plot
        #if self.out_path is None: self.figure.show()
        if self.config.show: self.figure.show()

        # Save the figure
        if self.out_path is not None: self.save_figure()

    # -----------------------------------------------------------------

    def add_legends(self):
        
        """
        Thisf unction ...
        :return: 
        """

        # Inform the user
        log.info("Adding the legends to the plot ...")

        # Loop over the legends
        for legend in self.legends:

            # Add to the plot window
            self.main_plot.add_artist(legend)

        # Legends for residual panels?
        if self.has_legends_residuals:

            # Loop over the residual legends
            for residual_plot, legends in zip(self.residual_plots, self.legends_residuals):

                # Add to the plot
                for legend in legends: residual_plot.add_artist(legend)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Close the figure, if it was saved or told to show: in other cases, the user probably wants to use the figure later
        if self.has_figure and (not self._external_figure) and (self.out_path is not None or self.config.show): self.figure.close()

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

    def set_defaults(self):

        """
        This function ...
        :return:
        """

        rcdefaults()

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

    # Add range
    plot_min = 0.9 * log_min_flux if log_min_flux > 0 else 1.1 * log_min_flux
    plot_max = 1.1 * log_max_flux if log_max_flux > 0 else 0.95 * log_max_flux

    #print(log_min_flux, log_max_flux)

    if numbers.is_invalid(plot_min): plot_min = None
    if numbers.is_invalid(plot_max): plot_max = None

    # Return
    return plot_min, plot_max

# -----------------------------------------------------------------

def process_errorbar(error, lolim_value=None, uplim_value=None, default_error=1.):

    """
    This function ...
    :param error:
    :param lolim_value:
    :param uplim_value:
    :param default_error:
    :return:
    """

    lolim = False
    uplim = False

    error_limits = [[abs(error.lower), abs(error.upper)]]

    # No lower limit
    if error.no_lower:

        if lolim_value is not None: low_value = lolim_value
        elif error.has_upper: low_value = abs(error_limits[0][1])
        else: low_value = default_error

        # Set
        error_limits[0][0] = low_value
        uplim = True

    # No upper value
    if error.no_upper:

        if uplim_value is not None: up_value = uplim_value
        elif error.has_lower: up_value = abs(error_limits[0][0])
        else: up_value = default_error

        # Set
        error_limits[0][1] = up_value
        lolim = True

    # Set variables
    lolims = [[lolim]]
    uplims = [[uplim]]

    # Return
    yerr = np.array(error_limits).T
    return yerr, lolims, uplims

# -----------------------------------------------------------------

def process_errorbar_for_value(value, error, lolim_value=None, uplim_value=None, lolim_abs_value=None,
                               uplim_abs_value=None, default_rel_error=None, default_error=None, min_allowed_value=None,
                               max_allowed_value=None, return_log=False):

    """
    Thisf unction ...
    :param value:
    :param error:
    :param lolim_value:
    :param uplim_value:
    :param lolim_abs_value:
    :param uplim_abs_value:
    :param default_rel_error:
    :param default_error:
    :param min_allowed_value:
    :param max_allowed_value:
    :param return_log:
    :return:
    """

    lowerr = abs(error.lower)
    uperr = abs(error.upper)

    value_min_lower = value - lowerr
    value_plus_upper = value + uperr

    # Initialize
    low_value = None
    up_value = None

    #if lower_flux <= 0:
    #    flux_lower_flux = float("-inf")
    #    uplims = [[True]]
    #else: flux_lower_flux = flux / lower_flux

    # Is upper limit?
    if min_allowed_value is not None and value_min_lower < min_allowed_value:
        uplim = True
        low_value = abs(min_allowed_value - value)
    elif error.no_lower: uplim = True
    else: uplim = False

    # Is lower limit?
    if max_allowed_value is not None and value_plus_upper > max_allowed_value:
        lolim = True
        up_value = max_allowed_value - value
    elif error.no_upper: lolim = True
    else: lolim = False

    # Set limits
    error_limits = [[lowerr, uperr]]

    # No lower limit
    if uplim:

        if low_value is None:

            if lolim_value is not None:
                if lolim_abs_value is not None: raise ValueError("Cannot specify both lolim value and lolum abs value")
                low_value = lolim_value
            elif lolim_abs_value is not None: low_value = abs(lolim_abs_value - value)
            elif error.has_upper: low_value = abs(error_limits[0][1])
            elif default_rel_error is not None:
                if default_error is not None: raise ValueError("Cannot specify both relative default error and default error")
                low_value = default_rel_error * value
            elif default_error is not None: low_value = default_error
            else: raise ValueError("Default error is not given")

        # Set
        error_limits[0][0] = low_value
        #uplim = True

    # No upper limit
    if lolim:

        if up_value is None:

            if uplim_value is not None:
                if uplim_abs_value is not None: raise ValueError("Cannot specify both uplim value and uplim abs value")
                up_value = uplim_value
            elif uplim_abs_value is not None: up_value = uplim_abs_value - value
            elif error.has_lower: up_value = abs(error_limits[0][0])
            elif default_rel_error is not None:
                if default_error is not None: raise ValueError("Cannot specify both relative default error and default error")
                up_value = default_rel_error * value
            elif default_error is not None: up_value = default_error
            else: raise ValueError("Default error is not given")

        # Set
        error_limits[0][1] = up_value
        #lolim = True

    # Set variables
    lolims = [[lolim]]
    uplims = [[uplim]]

    # Return for log scale
    if return_log:

        logvalue = np.log10(value)

        lowerr = error_limits[0][0]
        uperr = error_limits[0][1]

        value_min_lower = value - lowerr
        value_plus_upper = value + uperr

        logvalue_min_lower = np.log10(value_min_lower)
        logvalue_plus_upper = np.log10(value_plus_upper)

        # Set errobar
        error_limits = [[logvalue - logvalue_min_lower, logvalue_plus_upper - logvalue]]

    # Return
    yerr = np.array(error_limits).T
    return yerr, lolims, uplims

# -----------------------------------------------------------------

def join_observed_seds(seds):

    """
    This function ...
    :param seds:
    :return:
    """

    # Get SED objects and original labels
    if isinstance(seds, dict):
        labels = seds.keys()
        seds = seds.values()
    else: labels = None
    nseds = len(seds)

    # Create new SED
    unit = seds[0].unit
    if labels is not None: joined_sed = ObservedSED(photometry_unit=unit, extra_columns=[("Label", str, None, "original SED label")])
    else: joined_sed = ObservedSED(photometry_unit=unit)

    # Add the data points
    if labels is None: labels = [None] * nseds
    for label, sed in zip(labels, seds):

        # Loop over the data points
        for index in range(len(sed)):

            # Get filters and photometric values
            fltr = sed.get_filter(index)
            phot = sed.get_photometry(index, add_unit=True)

            # Add points
            if label is not None: extra = [label]
            else: extra = None
            joined_sed.add_point(fltr, phot, extra=extra)

    # Return
    return joined_sed

# -----------------------------------------------------------------

def interpolate_and_extrapolate_nans(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    # Get nans
    nan_indices = np.argwhere(np.isnan(y))
    if len(nan_indices) == 0: return y
    nan_indices = nan_indices[0]

    valid = np.isfinite(y)
    valid_x = x[valid]
    valid_y = y[valid]
    #print(valid_x)
    #print(valid_y)

    # Do the interpolation
    interpolation = interp1d(valid_x, valid_y, kind="linear", fill_value="extrapolate")

    # Get array of new x
    new_y = y.copy()

    # Loop over the nan indices
    for index in nan_indices:

        xi = x[index]
        interpolated = interpolation(xi)
        #original = y[index]
        #print(original, interpolated)
        new_y[index] = interpolated

    # Return the new data
    return new_y

# -----------------------------------------------------------------
