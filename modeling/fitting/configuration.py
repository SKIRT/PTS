#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.configuration Contains the FittingConfigurer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.simulation.skifile import LabeledSkiFile
from ...core.tools.logging import log
from ..config.parameters import definition as parameters_definition
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, Configuration, combine_configs
from ..config.parameters import parsing_types_for_parameter_types, unit_parsing_type
from ..config.parameters import default_units, possible_parameter_types_descriptions
from .run import FittingRun

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")

# -----------------------------------------------------------------

class FittingConfigurer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingConfigurer, self).__init__(config)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The default ranges
        self.default_ranges = dict()

        # The ski file template
        self.ski = None

        # The individual configurations
        self.parameters_config = None
        self.descriptions_config = None
        self.types_config = None
        self.units_config = None
        self.ranges_config = None
        self.filters_config = None

        # Additional settings
        self.settings = None

        # The final fitting config
        self.fitting_config = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the necessary input
        self.load_input()

        # 3. Get the fitting parameters
        self.set_parameters()

        # 4. Get parameter descriptions
        self.set_descriptions()

        # 5. Get parameter types
        self.set_types()

        # 6. Get parameter units
        self.set_units()

        # 7. Get the physical parameter ranges
        self.set_ranges()

        # 8. Get the fitting filters
        self.set_filters()

        # 9. Create the fitting configuration
        self.create_config()

        # 10. Adjust the labels of the template ski file
        self.adjust_labels()

        # 11. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

        # Create the fitting run
        self.create_fitting_run()

        # Get the default ranges
        self.default_ranges = kwargs.pop("default_ranges", dict())

        # Get configs as dicts
        if "parameters_config" in kwargs: self.parameters_config = kwargs.pop("parameters_config")
        if "descriptions_config" in kwargs: self.descriptions_config = kwargs.pop("descriptions_config")
        if "types_config" in kwargs: self.types_config = kwargs.pop("types_config")
        if "units_config" in kwargs: self.units_config = kwargs.pop("units_config")
        if "ranges_config" in kwargs: self.ranges_config = kwargs.pop("ranges_config")
        if "filters_config" in kwargs: self.filters_config = kwargs.pop("filters_config")

        # Set settings dict
        if "settings" in kwargs: self.settings = kwargs.pop("settings")

    # -----------------------------------------------------------------

    def create_fitting_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fitting run ...")

        # Create the run
        self.fitting_run = FittingRun(self.config.path, self.config.name, self.config.model_name)

        # Create the run directory
        fs.create_directory(self.fitting_run.path)

    # -----------------------------------------------------------------

    @property
    def run_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.name

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_name

    # -----------------------------------------------------------------

    @property
    def run_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.path

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input ...")

        # Load the template ski file
        self.load_template()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the labeled ski template file
        self.ski = LabeledSkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the free parameters ...")

        # Load or prompt for the parameters
        if self.config.parameters is not None: self.parameters_config = Configuration(free_parameters=self.config.parameters)
        elif isinstance(self.parameters_config, dict): pass
        else: self.prompt_parameters()

    # -----------------------------------------------------------------

    def prompt_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the free parameters ...")

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Free parameters", add_logging=False, add_cwd=False)

        # Create config
        self.parameters_config = setter.run(parameters_definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def set_descriptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter descriptions ...")

        # Load or prompt for the descriptions
        if self.config.descriptions is not None: self.descriptions_config = Configuration(descriptions=self.config.descriptions)
        elif isinstance(self.descriptions_config, dict): pass
        else: self.prompt_descriptions()

    # -----------------------------------------------------------------

    def prompt_descriptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter descriptions ...")

        # Create configuration definition
        definition = ConfigurationDefinition()
        for name in self.parameters_config.free_parameters: definition.add_required(name, "string", "description of the '" + name + "' parameter")

        # Create the configuration setter
        setter = InteractiveConfigurationSetter("Parameter descriptions", add_cwd=False, add_logging=False)

        # Get the config and set the descriptions configuration
        config = setter.run(definition, prompt_optional=False)
        self.descriptions_config = Configuration(descriptions=config)

    # -----------------------------------------------------------------

    def set_types(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter types ...")

        # load or prompt for the tyopes
        if self.config.types is not None: self.types_config = Configuration(types=self.config.types)
        elif isinstance(self.types_config, dict): pass
        else: self.prompt_types()

    # -----------------------------------------------------------------

    def prompt_types(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter types ...")

        # Create definition
        definition = ConfigurationDefinition()
        for name in self.parameters_config.free_parameters: definition.add_required(name, "string", "type of the '" + name + "' parameter", choices=possible_parameter_types_descriptions)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Parameter types", add_cwd=False, add_logging=False)

        # Create the config and set the types configuration
        config = setter.run(definition, prompt_optional=False)
        self.types_config = Configuration(types=config)

    # -----------------------------------------------------------------

    def set_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter units ...")

        # load or prompt for the units
        if self.config.units is not None: self.units_config = Configuration(units=self.config.units)
        elif isinstance(self.units_config, dict): pass
        else: self.prompt_units()

    # -----------------------------------------------------------------

    def prompt_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter units ...")

        # Create definition
        definition = ConfigurationDefinition()
        for name in self.parameters_config.free_parameters:

            # Get the type of quantity for this parameter
            parameter_type = self.types_config.types[name]

            # Don't ask for units for dimensionless quantities
            if parameter_type == "dimless": definition.add_fixed(name, name + " has no unit (dimensionless)", None)
            else: definition.add_optional(name, unit_parsing_type(parameter_type), "unit of the '" + name + "' parameter", default=default_units[parameter_type])

        # Create configuration setter
        setter = InteractiveConfigurationSetter("Parameter units", add_cwd=False, add_logging=False)

        # Create the config and set the units configuration
        config = setter.run(definition, prompt_optional=True)
        self.units_config = Configuration(units=config)

    # -----------------------------------------------------------------

    def set_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Ranges are given
        if len(self.config.ranges) != 0:

            self.ranges_config = dict()
            for label in self.config.ranges:
                self.ranges_config[label + "_range"] = self.config.ranges[label]

        elif isinstance(self.ranges_config, dict): pass
        else: self.prompt_ranges()

    # -----------------------------------------------------------------

    def prompt_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the parameter ranges ...")

        # Create the configuration
        definition = ConfigurationDefinition(write_config=False)

        # Add the options for the ranges
        for label in self.parameters_config.free_parameters:

            # Get the unit
            unit = self.units_config.units[label]
            #in_units_string = " (in " + unit + ")" if unit is not None else " (dimensionless)"
            units_info_string = " (don't forget the units!) " if unit is not None else " (dimensionless)"

            # Get the default range
            default_range = self.default_ranges[label] if label in self.default_ranges else None

            # Get the parsing type
            parsing_type = parsing_types_for_parameter_types[self.types_config.types[label]]

            # Get the description
            parameter_description = self.descriptions_config.descriptions[label]
            description = "range of " + parameter_description + units_info_string

            # Make optional or required depending on whether default is given
            if default_range is not None: definition.add_optional(label + "_range", parsing_type + "_range", description, default=default_range, convert_default=True)
            else: definition.add_required(label + "_range", parsing_type + "_range", description)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("free parameter ranges", add_logging=False, add_cwd=False)

        # Create config, get the range for each chosen free parameter
        self.ranges_config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def set_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the fitting filters ...")

        # Load or prompt for the filters
        if self.config.filters is not None: self.filters_config = Configuration(filters=self.config.filters)
        elif isinstance(self.filters_config, dict): pass
        else: self.prompt_filters()

    # -----------------------------------------------------------------

    def prompt_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the filters ...")

        # Create the configuration
        definition = ConfigurationDefinition(write_config=False)

        # Choose from all the possible filter names
        definition.add_required("filters", "string_list", "the filters for which to use the observed flux as reference for the fitting procedure", choices=self.sed_filter_names)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("filters", add_logging=False, add_cwd=False)

        # Create config, get the filter choices
        self.filters_config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def set_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting additional settings ...")

        # Prompt settings
        if self.settings is None: self.prompt_settings()

    # -----------------------------------------------------------------

    def prompt_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for additional settings ...")

        # ...
        self.settings = dict()

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the fitting run configuration ...")

        # Combine configs
        self.fitting_config = combine_configs(self.parameters_config, self.descriptions_config, self.types_config, self.units_config, self.ranges_config, self.filters_config)

        # Set additional settings
        for label in self.settings: self.fitting_config[label] = self.settings[label]

    # -----------------------------------------------------------------

    def adjust_labels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting labels of the free parameters in the template ski file ...")

        # Loop over the labels currently in the template ski file
        for label in self.ski.labels:

            # If the label is in the list of free parameter labels, skip (don't remove)
            if label in self.parameters_config.free_parameters: continue

            # Otherwise, delabel
            self.ski.delabel(label)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the fitting configuration
        self.write_config()

        # Write the ski file template
        self.write_ski()

        # Write the runs table
        self.write_table()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fitting configuration ...")

        # Write the configuration
        self.fitting_config.saveto(self.fitting_run.fitting_configuration_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the template ski file ...")

        # Save the ski file template
        self.ski.saveto(self.fitting_run.template_ski_path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the runs table ...")

        # Open the runs table, add the new run and save
        table = self.runs_table
        table.add_run(self.fitting_run)
        table.save()

# -----------------------------------------------------------------
