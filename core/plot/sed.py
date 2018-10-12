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

    # Show the figure after plotting
    show = kwargs.pop("show", None)

    # Create SED plotter
    plotter = SEDPlotter(kwargs)

    # Set TeX flag
    plotter.config.tex = tex

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
        wavelength_unit = sequences.get_single(wavelength_units)
    if unit is None:
        units = [seds[label].unit for label in seds]
        unit = sequences.get_single(units)

    # Set units
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

        # Legend
        self.for_legend_patches = None
        self.for_legend_parameters = None
        self.extra_legend = None

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
        self.model_options[label].residuals = kwargs.pop("residuals", True)
        self.model_options[label].ghost = kwargs.pop("ghost", False)
        self.model_options[label].above = kwargs.pop("above", None)
        self.model_options[label].above_name = kwargs.pop("above_name", None)

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

        # Plot the SEDs
        self.plot()

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
            else: nrespanels = 1

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
                if self.config.residual_reference == observations_reference: nrespanels = 1
                elif self.config.residual_reference == models_reference: nrespanels = self.nmodels_for_residuals
                else: raise ValueError("Invalid residual reference '" + self.config.residual_reference + "'")

            # Multiple observations
            else:

                # Determine number of plot rows
                if self.config.residual_reference == observations_reference: nrespanels = self.nobservations_as_reference
                elif self.config.residual_reference == models_reference: nrespanels = self.nmodels_for_residuals
                else: raise ValueError("Invalid residual reference")

        # Create
        self.main_plot, self.residual_plots = self.figure.create_sed_plots(nresiduals=nrespanels)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the SED plot ...")

        # No models
        if self.no_models: self.plot_no_models()

        # With models
        else: self.plot_with_models()

    # -----------------------------------------------------------------

    def plot_no_models(self):

        """
        This function ...
        :return:
        """

        # One observation
        if self.one_observation: self.plot_one_observation()

        # More observations
        else: self.plot_more_observations()

    # -----------------------------------------------------------------

    def plot_one_observation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting one observed SED ...")

        # Determine color map class
        colormap = plt.get_cmap("rainbow")

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
        fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

        # Create colors
        colors = colormap(np.linspace(0, 1, len(wavelengths)))

        # Plot
        self.draw_observation(instruments, bands, wavelengths, fluxes, errors, colors)

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

    def plot_more_observations(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Plotting multiple observed SEDs ...")

        # Make iterable from distinct colors
        different_colors = iter(dark_pretty_colors)

        # Markers for the unique labels
        markers = filled_markers[:len(self.unique_observation_point_labels)]

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
        reference_sed = self.observations[reference_sed_label]

        # Loop over the different observed SEDs
        for label in self.observations:

            # Debugging
            log.debug("Plotting the data points of the '" + label + "' observed SED ...")

            # Get the observed SED
            observation = self.observations[label]

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

            # Create color range
            observation_color = color_hex[next(different_colors)]
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
                marker = markers[self.unique_observation_point_labels.index(labels[k])]

                # Set zero error if other observations also have errors (BECAUSE CALLING MATPLOTLIB'S ERRORBAR AND AFTER THAT PLOT DOESN'T WORK APPARENTLY!!)
                error = errors[k]
                if self.has_any_observation_errors and error is None: error = ErrorBar.zero()

                # Ignore upper limits?
                if self.config.ignore_upper and error.only_upper: continue

                # Plot the flux data point on the main axis with the specified marker and color
                patch, new_label = self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], error, marker, observation_color, return_new=True)

                # Add the patch and label
                if patch is not None and new_label:

                    # Remove the color -> make a new patch?
                    patch = lines.Line2D([], [], marker=marker, markersize=7, label=labels[k], linewidth=0, markeredgecolor="black", markerfacecolor="white", markeredgewidth=1)

                    # Add patch
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

                    #yerr, lolims, uplims = process_errorbar(error)
                    yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)
                    self.residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=observation_color, markeredgecolor='black', ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

            # The next observation is not the first anymore
            first = False

            # If color gradient is used for observations instead of single colors ...
            #ax3 = next(colormap_axes)
            #ax3.imshow(gradient, aspect='auto', cmap=colormap)

            # Create rectangle for this observation
            rectangle = patches.Rectangle((0, 0), 1, 1, fc=observation_color)
            legend_rectangles.append(rectangle)
            rectangle_labels.append(label.replace("_", "\_"))

        # Extra legend: the different observations
        # fancybox=True makes the legend corners rounded
        observations_legend = self.main_plot.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Set
        self.for_legend_patches = legend_patches
        self.for_legend_parameters = legend_labels
        self.extra_legend = observations_legend

    # -----------------------------------------------------------------

    def plot_with_models(self):

        """
        This function ...
        :return:
        """

        # Only models
        if self.no_observations: self.plot_only_models()

        # One observation
        elif self.one_observation: self.plot_one_observation_with_models()

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

        # Keep
        counter = 0
        model_colors = dict()
        model_styles = dict()

        #line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        #line_colors_models = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        #line_styles_models_no_residuals = ["-"] * len(self.models)
        line_colors_models = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        for color in dark_pretty_colors:
            if color not in line_colors_models: line_colors_models.append(color)
        #line_colors_models = dark_pretty_colors

        line_styles_models = line_styles if self.nmodels_not_ghost <= len(line_styles) else ["-"] * len(line_colors_models)

        # Loop over the model SEDs
        for model_label in self.models:

            # Debugging
            log.debug("Plotting the '" + model_label + "' SED ...")

            #sed, plot_residuals, ghost = self.models[model_label]
            sed = self.models[model_label]
            #plot_residuals = self.model_options[model_label].residuals
            ghost = self.model_options[model_label].ghost
            above = self.model_options[model_label].above

            # Get fluxes, wavelengths and errors
            fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info, asarray=True)
            wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False, asarray=True)

            if above is not None:

                above_wavelengths = self.models[above].wavelengths(unit=self.config.wavelength_unit, add_unit=False, asarray=True)
                above_fluxes = self.models[above].photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info, asarray=True)

                #fluxes = fluxes + above_fluxes

                # Interpolate to same grid
                above_interpolation = interp1d(above_wavelengths, above_fluxes, kind='cubic')

                #print("min", above_interpolation.x[0])
                #print("max", above_interpolation.x[-1])
                x_new = wavelengths
                below_bounds = x_new < above_interpolation.x[0]
                above_bounds = x_new > above_interpolation.x[-1]
                #print(below_bounds)
                #print(above_bounds)

                valid = np.logical_not(below_bounds) * np.logical_not(above_bounds)
                #valid_indices = np.where(valid)
                wavelengths = wavelengths[valid]
                fluxes = fluxes[valid]

                # Determine residuals
                #print(above, above_wavelengths)
                #print(model_label, wavelengths)
                above_fluxes = above_interpolation(wavelengths)
                fluxes = fluxes + above_fluxes

            else: above_fluxes = None

            # Ghost
            if ghost:

                # Plot the model SED as a grey line (no errors)
                self.draw_model(self.main_plot, wavelengths, fluxes, "-", linecolor="lightgrey")
                model_colors[model_label] = "lightgrey"
                model_styles[model_label] = "-"

            # No
            # PLOT_RESIDUALS IS IRRELEVANT, RIGHT?
            else:

                # Get fluxes, wavelengths and errors
                errors = sed.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info) if sed.has_errors else None

                # Determine actual label
                if above is not None:

                    above_name = self.model_options[model_label].above_name
                    if above_name is not None: actual_label = above_name
                    else: actual_label = None

                else: actual_label = model_label.replace("_", "\_")

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles_models[counter], label=actual_label, errors=errors, linecolor=line_colors_models[counter])

                # Set colour
                model_colors[model_label] = line_colors_models[counter]
                model_styles[model_label] = line_styles_models[counter]

                counter += 1

            # Fill
            #print(wavelengths)
            #print(above_fluxes)
            #print(fluxes)
            if above is not None: self.main_plot.axes.fill_between(wavelengths, np.log10(above_fluxes), np.log10(fluxes), facecolor=model_colors[model_label], alpha=0.5, label=model_label.replace("_", "\_"))

        # Plot residual axis
        if self.has_residual_plots:

            # Pair of models
            if self.nmodels_for_residuals == 2:

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
                wavelengths_fluxes_a = [(wavelength, flux) for wavelength, flux in zip(wavelengths_a, fluxes_a) if min_wavelength < wavelength < max_wavelength]
                wavelengths_fluxes_b = [(wavelength, flux) for wavelength, flux in zip(wavelengths_b, fluxes_b) if min_wavelength < wavelength < max_wavelength]

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
                self.residual_plot.plot(wavelengths_a, residuals_a, linestyle=model_styles[label_a], color=model_colors[label_a], label=label_a)
                self.residual_plot.plot(wavelengths_a, residuals_b, linestyle=model_styles[label_b], color=model_colors[label_b], label=label_b)

            # More models
            elif self.nmodels_for_residuals > 2:

                # Get reference
                reference_model_label = self.model_labels_for_residuals[0]
                reference_model_sed = self.models[reference_model_label]
                reference_wavelengths = reference_model_sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                reference_fluxes = reference_model_sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

                # Plot reference
                residuals = [0.0] * len(reference_wavelengths)
                self.residual_plot.plot(reference_wavelengths, residuals, linestyle=model_styles[reference_model_label], color=model_colors[reference_model_label], label=reference_model_label)

                # Loop over the models to be plotted
                for model_label in self.model_labels_for_residuals:
                    if model_label == reference_model_label: continue

                    # Get color and style
                    linecolor = model_colors[model_label]
                    linestyle = model_styles[model_label]

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
    def conversion_info(self):
        info = dict()
        if self.config.distance is not None: info["distance"] = self.config.distance
        return info

    # -----------------------------------------------------------------

    def plot_one_observation_with_models(self):

        """
        This function ...
        :param:
        :return:
        """

        # Debugging
        log.debug("Plotting one observed SED with model SEDs ...")

        # Get the first (only) observation
        observation = self.observations[self.observations.keys()[0]]

        # Get wavelengths, fluxes, instruments, bands, errors
        wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
        fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
        instruments = observation.instruments()
        bands = observation.bands()
        errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

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

            # Get next color
            color = next(color_range)

            # Get marker
            marker = markers[unique_labels.index(labels[k])]

            # Plot on the main axis with the specified marker and color
            # axis, label, used_labels, wavelength, flux, error, marker, color, return_patch=False
            #print(labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)
            # GALEX [] 0.153507951158 0.0686478012556 None o [ 0.5  0.   1.   1. ]
            self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], errors[k], marker, color)

            # Observations as reference: plot at 0.0 (all in one panel)
            if self.config.residual_reference == observations_reference:

                # Get the (only) residual plot
                residual_plot = self.residual_plots[0]

                # Plot point at y=0.0 with errorbar on axis 2
                value = 0.0
                if errors[k] is not None:
                    error = errors[k] / fluxes[k] * 100.

                    #yerr, lolims, uplims = process_errorbar(error)
                    yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)
                    residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

            # Models as reference: plot points at the relative difference with the models (one panel for each model)
            elif self.config.residual_reference == models_reference:

                from astropy.units import Unit

                residual_plot_index = 0

                # Loop over the models
                for model_label in self.models:

                    #sed, plot_residuals, ghost = self.models[model_label]
                    sed = self.models[model_label]
                    plot_residuals = self.model_options[model_label].residuals
                    ghost = self.model_options[model_label].ghost

                    # Don't plot residuals for this model
                    if not plot_residuals: continue

                    #model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                    #sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                    actual_wavelength = wavelengths[k] * Unit(self.config.wavelength_unit)

                    # Wavelength out of bounds of model SED
                    if actual_wavelength < sed.min_wavelength or actual_wavelength > sed.max_wavelength:
                        residual_plot_index += 1
                        continue

                    model_value = sed.photometry_at(actual_wavelength, self.config.unit, add_unit=False, conversion_info=self.conversion_info, interpolate=self.config.interpolate_models_for_residuals)
                    rel_residual = (fluxes[k] - model_value) / model_value * 100.  # NORMALIZE TO MODEL VALUE

                    #print("MODEL VALUE:", model_value)
                    #print("FLUX:", fluxes[k])
                    #print("REL RESIDUAL:", rel_residual)

                    residual_plot = self.residual_plots[residual_plot_index]
                    #value = 0.0
                    if errors[k] is not None:

                        error = errors[k] / model_value * 100. # NORMALIZE TO MODEL VALUE

                        #yerr, lolims, uplims = process_errorbar(error)
                        yerr, lolims, uplims = process_errorbar_for_value(fluxes[k], error, lolim_abs_value=-75, uplim_value=75)
                        residual_plot.errorbar(wavelengths[k], rel_residual, yerr=yerr, fmt=marker, markersize=7, color=color, markeredgecolor='black', ecolor=color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

                    # Next residual plot
                    residual_plot_index += 1

        line_styles_models = line_styles
        line_colors_models = ['black'] * len(self.models)

        line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        for color in dark_pretty_colors:
            if color not in line_colors_models_no_residuals: line_colors_models_no_residuals.append(color)

        line_styles_models_no_residuals = ["-"] * len(self.models)

        counter = 0

        # Plot models on the residual axis, WHEN OBSERVATION IS THE REFERENCE
        for model_label in self.models:

            # Get entry
            #sed, plot_residuals, ghost = self.models[model_label]
            sed = self.models[model_label]
            plot_residuals = self.model_options[model_label].residuals
            ghost = self.model_options[model_label].ghost

            if not plot_residuals: continue

            if ghost:

                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                min_sed_wavelength = min(sed_wavelengths)
                max_sed_wavelength = max(sed_wavelengths)
                wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                self.residual_plot.plot(wavelengths_residuals, residuals, linestyle="-", color="lightgrey")

            elif self.config.residual_reference == "observations":

                #log_model = np.log10(sed.fluxes(unit=self.config.unit, add_unit=False))
                model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                f2 = interp1d(sed_wavelengths, model_fluxes, kind='cubic')

                min_sed_wavelength = min(sed_wavelengths)
                max_sed_wavelength = max(sed_wavelengths)
                wavelengths_fluxes_residuals = sorted([(wavelength, flux) for wavelength, flux in zip(wavelengths, fluxes) if min_sed_wavelength < wavelength < max_sed_wavelength], key=itemgetter(0))
                wavelengths_residuals = [item[0] for item in wavelengths_fluxes_residuals]
                fluxes_residuals = [item[1] for item in wavelengths_fluxes_residuals]
                residuals = -(fluxes_residuals - f2(wavelengths_residuals)) / fluxes_residuals * 100.

                self.residual_plot.plot(wavelengths_residuals, residuals, linestyle=line_styles_models[counter], color=line_colors_models[counter], label='model')

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
        for model_label in self.models:

            # Get entry
            #sed, plot_residuals, ghost = self.models[model_label]
            sed = self.models[model_label]
            plot_residuals = self.model_options[model_label].residuals
            ghost = self.model_options[model_label].ghost

            if self.nmodels == 1: model_label = "model"

            if ghost:

                log_model = np.log10(sed.photometry(unit=self.config.unit, add_unit=False))
                self.main_plot.plot(sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False), log_model, ls="-", color="lightgrey")

            elif plot_residuals:

                fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                log_model = np.log10(fluxes)

                # Adjust extrema
                if self.adjust_wavelength_minmax_models: self.adjust_wavelength_minmax_wavelengths(wavelengths)
                if self.adjust_photometry_minmax_models: self.adjust_photometry_minmax_fluxes(fluxes)

                # Plot
                self.main_plot.plot(wavelengths, log_model, linestyle=line_styles_models[counter], color=line_colors_models[counter], label=model_label.replace("_", "\_"))

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

                fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
                log_model = np.log10(fluxes)
                wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                # Adjust extrema
                if self.adjust_wavelength_minmax_models: self.adjust_wavelength_minmax_wavelengths(wavelengths)
                if self.adjust_photometry_minmax_models: self.adjust_photometry_minmax_fluxes(fluxes)

                # Plot
                self.main_plot.plot(wavelengths, log_model, linestyle=line_styles_models_no_residuals[counter_no_residuals], color=line_colors_models_no_residuals[counter_no_residuals], label=model_label.replace("_", "\_"))

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

    # -----------------------------------------------------------------

    def plot_more_observations_with_models(self):

        """
        This function ...
        :return:
        """

        # Debuggging
        log.debug("Plotting multiple observed SEDs with models ...")

        # Make iterable from distinct colors
        different_colors = iter(dark_pretty_colors)

        # Markers for the unique labels
        markers = filled_markers[:len(self.unique_observation_point_labels)]

        # Used labels
        used_labels = []

        line_colors_models_no_residuals = ["r", "lawngreen", "blueviolet", "deepskyblue", "orange"]
        line_styles_models_no_residuals = ["-"] * len(self.models)

        legend_patches = []
        legend_labels = []

        legend_rectangles = []
        rectangle_labels = []

        # Remember colors for each observation
        observation_colors = dict()

        # Loop over the different observed SEDs
        observation_index = 0
        observation_reference_index = 0
        for label in self.observations:

            # Get the observed SED
            observation = self.observations[label]

            # Get options
            as_reference = self.observation_options[label].as_reference

            # Get wavelengths, fluxes, instruments, bands, errors
            wavelengths = observation.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            fluxes = observation.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
            instruments = observation.instruments()
            bands = observation.bands()
            errors = observation.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)

            # Create color range
            observation_color = color_hex[next(different_colors)]
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
                marker = markers[self.unique_observation_point_labels.index(labels[k])]

                # Set zero error if other observations also have errors (BECAUSE CALLING MATPLOTLIB'S ERRORBAR AND AFTER THAT PLOT DOESN'T WORK APPARENTLY!!)
                error = errors[k]
                if self.has_any_observation_errors and error is None: error = ErrorBar.zero()

                # Plot the flux data point on the main axis with the specified marker and color
                patch, new_label = self.plot_wavelength(self.main_plot, labels[k], used_labels, wavelengths[k], fluxes[k], error, marker, observation_color, return_new=True)

                # Add the patch and label
                if patch is not None and new_label:

                    # Remove color
                    # Remove the color -> make a new patch?
                    patch = lines.Line2D([], [], marker=marker, markersize=7, label=labels[k], linewidth=0,
                                         markeredgecolor="black", markerfacecolor="white", markeredgewidth=1)

                    # Add the patch and label
                    legend_patches.append(patch)
                    legend_labels.append(labels[k])

                # Observations as reference: plot at 0.0 (all in one panel)
                if self.config.residual_reference == observations_reference:

                    #print(label, observation_index, wavelengths[k])
                    if as_reference:

                        # Get the residual plot
                        #print(observation_reference_index, len(self.residual_plots))
                        residual_plot = self.residual_plots[observation_reference_index]
                        #print(residual_plot)

                        # Plot measurement points on residual plot
                        value = 0.0
                        #if errors[k] is not None:

                        if errors[k] is not None:

                            # Process and plot errorbar
                            error = errors[k] / fluxes[k] * 100.

                            #print(label, wavelengths[k], observation_color, observation_index, residual_plot._plot)

                            #yerr, lolims, uplims = process_errorbar(error)
                            yerr, lolims, uplims = process_errorbar_for_value(value, error, lolim_abs_value=-75, uplim_value=75)

                            residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=observation_color, markeredgecolor='black', ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

                        else:
                            #residual_plot.scatter(wavelengths[k], value, marker=marker, s=7, c=observation_color, edgecolors="black")
                            yerr = [[0], [0]]
                            residual_plot.errorbar(wavelengths[k], value, yerr=yerr, fmt=marker, markersize=7, color=observation_color, markeredgecolor='black', ecolor=observation_color, capthick=2, lolims=lolims, uplims=uplims, capsize=2)

                # Models as reference: plot points at the relative difference with the models (one panel for each model)
                elif self.config.residual_reference == models_reference:

                    from astropy.units import Unit

                    residual_plot_index = 0

                    # Loop over the models
                    for model_label in self.models:

                        # Get entry
                        #sed, plot_residuals, ghost = self.models[model_label]
                        sed = self.models[model_label]
                        plot_residuals = self.model_options[model_label].residuals
                        ghost = self.model_options[model_label].ghost

                        # Don't plot residuals for this model
                        if not plot_residuals: continue

                        # model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False)
                        # sed_wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)

                        actual_wavelength = wavelengths[k] * Unit(self.config.wavelength_unit)

                        # Wavelength out of bounds of model SED
                        if actual_wavelength < sed.min_wavelength or actual_wavelength > sed.max_wavelength:
                            residual_plot_index += 1
                            continue

                        #print(model_label, sed.unit, sed.unit.physical_type)
                        model_value = sed.photometry_at(actual_wavelength, self.config.unit, add_unit=False, conversion_info=self.conversion_info, interpolate=self.config.interpolate_models_for_residuals)
                        rel_residual = (fluxes[k] - model_value) / model_value * 100.  # NORMALIZE TO MODEL VALUE

                        # Get the residual plot for this model
                        residual_plot = self.residual_plots[residual_plot_index]

                        # value = 0.0
                        if errors[k] is None: errorbar = ErrorBar.zero()
                        else: errorbar = errors[k]

                        error = errorbar / model_value * 100.  # NORMALIZE TO MODEL VALUE

                        #yerr, lolims, uplims = process_errorbar(error)
                        yerr, lolims, uplims = process_errorbar_for_value(fluxes[k], error, lolim_abs_value=-75, uplim_value=75)

                        # Plot on residual axes
                        residual_plot.errorbar(wavelengths[k], rel_residual, yerr=yerr, fmt=marker,
                                               markersize=7, color=observation_color, markeredgecolor='black', ecolor=observation_color,
                                               capthick=2, lolims=lolims, uplims=uplims, capsize=2)

                        # Next residual plot
                        residual_plot_index += 1

                # Invalid
                else: raise RuntimeError("Invalid value for 'residual_reference'")

            # IF OBSERVATIONS ARE THE REFERENCE FOR THE RESIDUAL AXES
            if self.config.residual_reference == observations_reference:

                if as_reference:

                    # Get the residual plot
                    residual_plot = self.residual_plots[observation_reference_index]

                    # Residuals
                    counter = 0
                    for model_label in self.models:

                        # Get entry
                        #sed, plot_residuals, ghost = self.models[model_label]
                        sed = self.models[model_label]
                        plot_residuals = self.model_options[model_label].residuals
                        ghost = self.model_options[model_label].ghost

                        if not plot_residuals: continue

                        model_fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
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
                            residual_plot.plot(wavelengths_residuals, residuals, ls=line_styles[counter], color='black', label=model_label.replace("_", "\_"))
                            counter += 1

            # Create rectangle for this observation
            rectangle = patches.Rectangle((0, 0), 1, 1, fc=observation_color)
            legend_rectangles.append(rectangle)
            rectangle_labels.append(label.replace("_", "\_"))

            # Increment the observation index
            observation_index += 1
            if as_reference: observation_reference_index += 1

        ## END OF LOOP OVER OBSERVATIONS

        # Add model SEDs
        counter = 0
        counter_no_residuals = 0
        for model_label in self.models:

            # Get entry
            #sed, plot_residuals, ghost = self.models[model_label]
            sed = self.models[model_label]
            plot_residuals = self.model_options[model_label].residuals
            ghost = self.model_options[model_label].ghost

            # Get fluxes, wavelengths and errors
            fluxes = sed.photometry(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info)
            wavelengths = sed.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
            errors = sed.errors(unit=self.config.unit, add_unit=False, conversion_info=self.conversion_info) if sed.has_errors else None

            if ghost:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, "-", linecolor="lightgrey", adjust_extrema=False, errors=errors)

            elif plot_residuals:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles[counter], linecolor="black", label=model_label.replace("_", "\_"), adjust_extrema=False, errors=errors)
                counter += 1

            else:

                # Plot the model SED as a line (with errors if present)
                self.draw_model(self.main_plot, wavelengths, fluxes, line_styles_models_no_residuals[counter_no_residuals], linecolor=line_colors_models_no_residuals[counter_no_residuals], label=model_label, adjust_extrema=False, errors=errors)
                counter_no_residuals += 1

        # Extra legend: the different observations
        observations_legend = self.main_plot.legend(legend_rectangles, rectangle_labels, loc='upper left', shadow=False, fontsize=11, ncol=3)

        # Set
        self.for_legend_patches = legend_patches
        self.for_legend_parameters = legend_labels
        self.extra_legend = observations_legend

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
            else:
                #print("d", label, new_label)
                patch = axis.plot(wavelength, np.log10(flux), marker=marker, markersize=markersize, color=color, markeredgecolor=markeredgecolor)

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

            log_bottom = np.log10(bottom)
            log_top = np.log10(top)

            axis.fill_between(errors_wavelengths, log_bottom, log_top, where=log_top<=log_bottom, facecolor=errorcolor, edgecolor=errorcolor, interpolate=True, alpha=0.5)
            axis.plot([], [], color=errorcolor, linewidth=10, label='spread')

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
        if len(self.residual_plots) > 0:
            #for j in range(len(self.residual_plots)-1):
                #self.residual_plots[j].hide_xticks()
            last_residual_plot_index = len(self.residual_plots)-1
            self.residual_plots[last_residual_plot_index].set_xticks(fontsize=self.config.plot.ticks_fontsize)

        # Set x label of the last residual plot
        xlabel = r"Wavelength $\lambda\,[\mu \mathrm{m}]$"
        #print(xlabel)
        if len(self.residual_plots) > 0: self.residual_plots[len(self.residual_plots)-1].set_xlabel(xlabel, fontsize='large')
        else: self.main_plot.set_xlabel(xlabel, fontsize='large')

        # Set log x scale
        self.main_plot.set_xscale('log')

        # Set ticks
        #if len(self.residual_plots) == 0: self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        self.main_plot.set_xticks(fontsize=self.config.plot.ticks_fontsize)
        self.main_plot.set_yticks(fontsize=self.config.plot.ticks_fontsize) # minor = False is possible

        # Add axis labels and a legend
        y_label = r"Log $" + self.config.unit.latex_string + r"$"
        #print(y_label)
        self.main_plot.set_ylabel(y_label, fontsize='large')

        # SHOW TICK LABELS
        #print(self.main_plot.xtick_labels, self.main_plot.ytick_labels)
        #for res in self.residual_plots: print(res.xtick_labels, res.ytick_labels)

        # Set position
        # X: only for main plot
        if self.config.xaxis_position == "top":
            self.main_plot.xaxis.set_ticks_position("top")
            self.main_plot.xaxis.set_label_position("top")

        # Y: for main and residual plots
        if self.config.yaxis_position == "right":
            for plot in [self.main_plot] + self.residual_plots:
                plot.yaxis.set_ticks_position("right")
                plot.yaxis.set_label_position("right")

        # Set grid
        #self.figure.set_grid(self.config.plot, which="both")
        self.main_plot.set_grid(self.config.plot, which="both")
        for plot in self.residual_plots: plot.set_grid(self.config.plot, which="both")

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

        legends = []

        # Add the legend
        if self.for_legend_patches is not None:

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
            #legend = self.main_plot.legend([l[0] for l in for_legend_patches], for_legend_parameters, **legend_properties)
            legend = self.main_plot.legend(self.for_legend_patches, self.for_legend_parameters, **legend_properties)
            legends.append(legend)

            # Extra legend
            if self.extra_legend is not None:

                self.main_plot.add_artist(self.extra_legend)
                legends.append(self.extra_legend)

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

        # Add title if requested
        #print(self.title)
        if self.title is not None and not self._external_figure: self.figure.set_title(self.title) #self.figure.figure.suptitle("\n".join(wrap(self.title, 60)))

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
