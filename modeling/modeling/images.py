#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.sed Contains the ImagesModeler class.
#  Performs radiative transfer modeling using genetic algorithms based on observed images

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.images import ImagesFittingInitializer
from .base import ModelerBase
from ..component.sed import get_ski_template
from ...core.basics.range import IntegerRange, QuantityRange
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter
from ..core.environment import ImagesModelingEnvironment
from ..component.images import get_ski_input_path
from ..build.models.images import ImagesModelBuilder
from ..build.representations.images import ImagesRepresentationBuilder
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class ImagesModeler(ModelerBase):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImagesModeler, self).__init__(*args, **kwargs)

        # Optional configs for the fitting configurer
        self.descriptions_config = None
        self.types_config = None
        self.units_config = None
        self.ranges_config = None
        self.filters_config = None
        self.genetic_config = None
        self.grid_config = None

        # Configuration for the fitting initializer
        self.initialize_config = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the data
        self.load_data()

        # Build model
        self.build_model()

        # Build representation
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
        super(ImagesModeler, self).setup(**kwargs)

        # Load the modeling environment
        self.environment = ImagesModelingEnvironment(self.modeling_path)

        # Set configs for the fitting configurer
        if "descriptions_config" in kwargs: self.descriptions_config = kwargs.pop("descriptions_config")
        if "types_config" in kwargs: self.types_config = kwargs.pop("types_config")
        if "units_config" in kwargs: self.units_config = kwargs.pop("units_config")
        if "ranges_config" in kwargs: self.ranges_config = kwargs.pop("ranges_config")
        if "filters_config" in kwargs: self.filters_config = kwargs.pop("filters_config")
        if "genetic_config" in kwargs: self.genetic_config = kwargs.pop("genetic_config")
        if "grid_config" in kwargs: self.grid_config = kwargs.pop("grid_config")

        # Config for the fitting initializer
        if "initialize_config" in kwargs: self.initialize_config = kwargs.pop("initialize_config")

        # CHECK IF RERUN IS DEFINED, IF SO, REMOVE COMMANDS FROM THE MODELLING HISTORY
        if self.config.rerun is not None: self.set_rerun()

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input data ...")

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
        builder = ImagesModelBuilder(config)

        # Run
        with self.register(builder): builder.run()

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
        builder = ImagesRepresentationBuilder(config)

        # Run
        with self.write_log(builder), self.register(builder), self.write_config(builder): builder.run()

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

        # Set the free parameters
        config["parameters"] = free_parameter_names

        # Load the filters
        if self.filters_config is None: config["filters"] = self.environment.image_filter_names

        # Set fitting run name and model name
        config["name"] = self.fitting_run_name
        config["model_name"] = self.model_name

        # Ask for and set the fitting method
        config["fitting_method"] = self.fitting_method

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Set the input dictionary
        input_dict = dict()
        input_dict["descriptions_config"] = self.descriptions_config
        input_dict["types_config"] = self.types_config
        input_dict["units_config"] = self.units_config
        input_dict["ranges_config"] = self.ranges_config
        input_dict["filters_config"] = self.filters_config
        input_dict["genetic_config"] = self.genetic_config
        input_dict["grid_config"] = self.grid_config
        input_dict["settings"] = self.config.fitting_settings

        # Run the fitting configurer
        with self.write_log(configurer), self.register(configurer), self.write_config(configurer), self.write_input(configurer, **input_dict): configurer.run(**input_dict)

            #configurer.run(descriptions_config=self.descriptions_config, types_config=self.types_config,
            #               units_config=self.units_config, ranges_config=self.ranges_config, filters_config=self.filters_config,
            #               genetic_config=self.genetic_config, grid_config=self.grid_config, settings=self.config.fitting_settings)

        # Set the parameter ranges
        if self.ranges_config is not None:
            self.parameter_ranges = dict()
            for parameter_name in free_parameter_names:
                # Get the range
                parameter_range = self.ranges_config[parameter_name + "_range"]
                # Debugging
                log.debug("Setting the range of the '" + parameter_name + "' parameter to '" + tostr(parameter_range, fancy=True) + "' for the parameter exploration ...")
                # Set the range
                self.parameter_ranges[parameter_name] = parameter_range

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
        initializer = ImagesFittingInitializer()

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Load the current ski template
        ski = get_ski_template(self.modeling_path)

        # Get the original ski input path
        ski_input_path = get_ski_input_path(self.modeling_path)

        # Create a definition
        definition = ConfigurationDefinition()
        default_npackages = max(int(1e4), int(ski.packages() / 10))
        definition.add_optional("npackages", "positive_integer", "the number of photon packages per wavelength for the initial generation", default=default_npackages)
        definition.add_flag("selfabsorption", "enable dust self-absorption", default=ski.dustselfabsorption())
        definition.add_flag("transient_heating", "enable transient heating", default=ski.transientheating())

        # Add option for the range of the number of wavelengths
        nwavelengths = ski.get_nwavelengths(ski_input_path)
        min_nwavelengths = max(int(0.1 * nwavelengths), 45)
        max_nwavelengths = max(5 * nwavelengths, 5 * min_nwavelengths)
        default_nwavelengths_range = IntegerRange(min_nwavelengths, max_nwavelengths)
        definition.add_optional("nwavelengths_range", "integer_range", "range for the number of wavelengths to vary over the generations", default=default_nwavelengths_range)
        definition.add_optional("ngrids", "positive_integer", "number of wavelength grids to be generated", default=10)
        definition.add_flag("add_emission_lines", "add additional points to the wavelength grids to sample important dust/gas emission lines", default=False)
        default_wavelength_range = QuantityRange(ski.get_min_wavelength(ski_input_path), ski.get_max_wavelength(ski_input_path))
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
        with self.write_log(initializer), self.register(initializer), self.write_config(initializer): initializer.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
