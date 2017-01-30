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
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, Configuration
from ..config.parameters import parsing_types_for_parameter_types, unit_parsing_type
from ..config.parameters import default_units, possible_parameter_types_descriptions

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
        self.get_parameters()

        # 4. Get parameter descriptions
        self.get_descriptions()

        # 5. Get parameter types
        self.get_types()

        # 6. Get parameter units
        self.get_units()

        # 7. Get the physical parameter ranges
        self.get_ranges()

        # 8. Get the fitting filters
        self.get_filters()

        # 9. Adjust the labels of the template ski file
        self.adjust_labels()

        # 10. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

        # Get the default ranges
        self.default_ranges = kwargs.pop("default_ranges", dict())

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

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

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        # The free parameters are specified
        if self.config.parameters is not None: self.parameters_config = Configuration(free_parameters=self.config.parameters)
        else:

            # Create configuration setter
            setter = InteractiveConfigurationSetter("Free parameters", add_logging=False, add_cwd=False)

            # Create config
            self.parameters_config = setter.run(parameters_definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def get_descriptions(self):

        """
        This function ...
        :return:
        """

        if self.config.descriptions is not None: self.descriptions_config = Configuration(descriptions=self.config.descriptions)
        else:

            # Create configuration definition
            definition = ConfigurationDefinition()
            for name in self.parameters_config.free_parameters: definition.add_required(name, "string", "description of the '" + name + "' parameter")

            # Create the configuration setter
            setter = InteractiveConfigurationSetter("Parameter descriptions and types", add_cwd=False, add_logging=False)

            # Get the config and set the descriptions configuration
            config = setter.run(definition, prompt_optional=False)
            self.descriptions_config = Configuration(descriptions=config)

    # -----------------------------------------------------------------

    def get_types(self):

        """
        This function ...
        :return:
        """

        if self.config.types is not None: self.types_config = Configuration(types=self.config.types)
        else:

            # Create definition
            definition = ConfigurationDefinition()
            for name in self.parameters_config.free_parameters: definition.add_required(name, "string", "type of the '" + name + "' parameter", choices=possible_parameter_types_descriptions)

            # Create configuration setter
            setter = InteractiveConfigurationSetter("Parameter types", add_cwd=False, add_logging=False)

            # Create the config and set the types configuration
            config = setter.run(definition, prompt_optional=False)
            self.types_config = Configuration(types=config)

    # -----------------------------------------------------------------

    def get_units(self):

        """
        This function ...
        :return:
        """

        if self.config.units is not None: self.units_config = Configuration(units=self.config.units)
        else:

            # Create definition
            definition = ConfigurationDefinition()
            for name in self.parameters_config.free_parameters:

                # Get the type of quantity for this parameter
                parameter_type = self.types_config.types[name]

                # Don't ask for units for dimensionless quantities
                if parameter_type == "dimless": definition.add_fixed(name, name + " has no unit (dimensionless)", None)
                else:
                    print(parameter_type)
                    print(default_units)
                    print(unit_parsing_type(parameter_type))
                    print(default_units[parameter_type])
                    definition.add_optional(name, unit_parsing_type(parameter_type), "unit of the '" + name + "' parameter", default=default_units[parameter_type])

            # Create configuration setter
            setter = InteractiveConfigurationSetter("Parameter units", add_cwd=False, add_logging=False)

            # Create the config and set the units configuration
            config = setter.run(definition, prompt_optional=True)
            self.units_config = Configuration(units=config)

    # -----------------------------------------------------------------

    def get_ranges(self):

        """
        This function ...
        :return:
        """

        # If ranges are given
        if len(self.config.ranges) != 0: #self.ranges_config = self.config.ranges

            self.ranges_config = dict()
            for label in self.config.ranges:
                self.ranges_config[label + "_range"] = self.config.ranges[label]

        # Ranges are not given
        else:

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

    def get_filters(self):

        """
        This function ...
        :return:
        """

        if self.config.filters is not None: self.filters_config = {"filters": self.config.filters}
        else:

            # Create the configuration
            definition = ConfigurationDefinition(write_config=False)

            # Choose from all the possible filter names
            definition.add_required("filters", "string_list", "the filters for which to use the observed flux as reference for the fitting procedure", choices=self.sed_filter_names)

            # Create configuration setter
            setter = InteractiveConfigurationSetter("filters", add_logging=False, add_cwd=False)

            # Create config, get the filter choices
            self.filters_config = setter.run(definition, prompt_optional=False)

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

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fitting configuration file ...")

        # Combine configs
        self.fitting_config = combine_configs(self.parameters_config, self.descriptions_config, self.types_config, self.units_config, self.ranges_config, self.filters_config)

        # Write the configuration
        self.fitting_config.saveto(self.fitting_configuration_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Save the ski file template
        self.ski.saveto(self.template_ski_path)

# -----------------------------------------------------------------

def combine_configs(*args):

    """
    This function ...
    :param args:
    :return:
    """

    # Initialize a new configuration
    config = Configuration()

    for cfg in args:
        for label in cfg: config[label] = cfg[label]

    # Return the resulting configuration
    return config

# -----------------------------------------------------------------
