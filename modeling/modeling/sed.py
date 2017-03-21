#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.sed Contains the SEDModeler class.
#  Perform radiative transfer modeling using genetic algorithms on a simple SED

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.sed import SEDFittingInitializer
from .base import ModelerBase
from ..component.sed import get_ski_template, get_observed_sed, get_sed_plot_path
from ...core.basics.range import IntegerRange, QuantityRange
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter
from ...core.plot.sed import SEDPlotter
from ..build.sed import SEDModelBuilder
from ..build.sedrepresentation import SEDRepresentationBuilder
from ..core.environment import get_ski_input_paths

# -----------------------------------------------------------------

class SEDModeler(ModelerBase):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(SEDModeler, self).__init__(config, interactive)

        # Optional configs for the fitting configurer
        #self.parameters_config = None
        self.descriptions_config = None
        self.types_config = None
        self.units_config = None
        self.ranges_config = None
        self.filters_config = None
        self.genetic_config = None

        # Configuration for the fitting initializer
        self.initialize_config = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the data
        self.load_data()

        # 3. Build the model
        self.build_model()

        # Build the model representation
        self.build_representation()

        # 3. Do the fitting
        self.fit()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SEDModeler, self).setup(**kwargs)

        # Set configs for the fitting configurer
        #if "parameters_config" in kwargs: self.parameters_config = kwargs.pop("parameters_config")
        if "descriptions_config" in kwargs: self.descriptions_config = kwargs.pop("descriptions_config")
        if "types_config" in kwargs: self.types_config = kwargs.pop("types_config")
        if "units_config" in kwargs: self.units_config = kwargs.pop("units_config")
        if "ranges_config" in kwargs: self.ranges_config = kwargs.pop("ranges_config")
        if "filters_config" in kwargs: self.filters_config = kwargs.pop("filters_config")
        if "genetic_config" in kwargs: self.genetic_config = kwargs.pop("genetic_config")

        # Config for the fitting initializer
        if "initialize_config" in kwargs: self.initialize_config = kwargs.pop("initialize_config")

        # Set ranges dict
        #if self.ranges_config is not None and self.parameters_config is not None:
        #    self.parameter_ranges = dict()
        #    for parameter_name in self.parameters_config.free_parameters:
        #        range = self.ranges_config[parameter_name + "_range"]
        #        self.parameter_ranges[parameter_name] = range

        # Set ranges dict

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input data ...")

        # Plot SED
        if "plot_sed" not in self.history: self.plot_sed()

    # -----------------------------------------------------------------

    def build_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the model ...")

        # Create configuration
        config = dict()
        config["name"] = self.model_name

        # Create the builder
        builder = SEDModelBuilder(config)

        # Run
        builder.run()

    # -----------------------------------------------------------------

    def build_representation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the representation ...")

        # Create configuration
        config = dict()
        config["name"] = self.representation_name
        config["model_name"] = self.model_name

        # Create the representation
        builder = SEDRepresentationBuilder(config)

        # Run
        builder.run()

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Add an entry to the history
        self.history.add_entry("plot_sed")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add the observed SED
        sed = get_observed_sed(self.modeling_path)
        plotter.add_sed(sed, "Observations")

        # Run the plotter
        plotter.run(output=get_sed_plot_path(self.modeling_path))

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def configure_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the fitting ...")

        # Create configuration for the FittingConfigurer
        config = dict()

        # Load the ski template, get the free parameters
        ski = get_ski_template(self.config.path)
        free_parameter_names = ski.labels

        # Set free parameters
        config["parameters"] = free_parameter_names

        # Load the SED, get the fitting filters
        if self.filters_config is not None:
            sed = get_observed_sed(self.config.path)
            fitting_filter_names = sed.filter_names()
            config["filters"] = fitting_filter_names

        # Set fitting run name and model name
        config["name"] = self.fitting_run_name
        config["model_name"] = self.model_name

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Add an entry to the history
        self.history.add_entry(FittingConfigurer.command_name())

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Run the fitting configurer
        #configurer.run(parameters_config=self.parameters_config, descriptions_config=self.descriptions_config,
        #               types_config=self.types_config, units_config=self.units_config, ranges_config=self.ranges_config,
        #               filters_config=self.filters_config, genetic_config=self.genetic_config, settings=self.config.fitting_settings)
        configurer.run(descriptions_config=self.descriptions_config, types_config=self.types_config,
                       units_config=self.units_config, ranges_config=self.ranges_config, filters_config=self.filters_config,
                       genetic_config=self.genetic_config, settings=self.config.fitting_settings)

        # Get the parameter ranges
        # Set ranges dict
        if self.ranges_config is not None:
            self.parameter_ranges = dict()
            for parameter_name in free_parameter_names:

                # Get the range
                range = self.ranges_config[parameter_name + "_range"]

                # Debugging
                log.debug("Setting the range of the '" + parameter_name + "' parameter to '" + str(range) + "' for the parameter exploration ...")

                # Set the range
                self.parameter_ranges[parameter_name] = range

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the fitting ...")

        # Create configuration
        config = dict()
        config["name"] = self.fitting_run_name

        # Create the fitting initializer
        initializer = SEDFittingInitializer()

        # Add an entry to the history
        self.history.add_entry(SEDFittingInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Load the current ski template
        ski = get_ski_template(self.modeling_path)

        # Load the ski input paths
        ski_input_paths = get_ski_input_paths(self.modeling_path)

        # Create a definition
        definition = ConfigurationDefinition()
        default_npackages = max(int(1e4), int(ski.packages() / 10))
        definition.add_optional("npackages", "positive_integer", "the number of photon packages per wavelength for the initial generation", default=default_npackages)
        definition.add_flag("selfabsorption", "enable dust self-absorption", default=ski.dustselfabsorption())
        definition.add_flag("transient_heating", "enable transient heating", default=ski.transientheating())

        # Add option for the range of the number of wavelengths
        nwavelengths = ski.get_nwavelengths(ski_input_paths)
        min_nwavelengths = max(int(0.1 * nwavelengths), 45)
        max_nwavelengths = max(5 * nwavelengths, 5 * min_nwavelengths)
        default_nwavelengths_range = IntegerRange(min_nwavelengths, max_nwavelengths)
        definition.add_optional("nwavelengths_range", "integer_range", "range for the number of wavelengths to vary over the generations", default=default_nwavelengths_range)
        definition.add_optional("ngrids", "positive_integer", "number of wavelength grids to be generated", default=10)
        definition.add_flag("add_emission_lines", "add additional points to the wavelength grids to sample important dust/gas emission lines", default=False)
        default_wavelength_range = QuantityRange(ski.get_min_wavelength(ski_input_paths), ski.get_max_wavelength(ski_input_paths))
        definition.add_optional("wavelength_range", "quantity_range", "wavelength range for all wavelength grids", default=default_wavelength_range)

        # Create the setter
        if self.initialize_config is None:
            setter = InteractiveConfigurationSetter("Initialization of the ski template", add_cwd=False, add_logging=False)
            config = setter.run(definition, prompt_optional=True)
        else:
            setter = DictConfigurationSetter(self.initialize_config, "Initialization of the ski template")
            config = setter.run(definition)

        # Set fixed settings for the ski model
        initializer.config.npackages = config.npackages
        initializer.config.selfabsorption = config.selfabsorption
        initializer.config.transient_heating = config.transient_heating

        # Set options for the wavelength grids
        initializer.config.wg.npoints_range = config.nwavelengths_range
        initializer.config.wg.ngrids = config.ngrids
        initializer.config.wg.add_emission_lines = config.add_emission_lines
        initializer.config.wg.min_wavelength = config.wavelength_range.min
        initializer.config.wg.max_wavelength = config.wavelength_range.max

        # OPTIONS FOR THE DUST GRID NOT RELEVANT FOR SED MODELING (YET)

        # Run the fitting initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
