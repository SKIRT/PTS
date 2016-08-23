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
from ..config.parameters import descriptions, types_and_ranges
from ..config.parameters import units as parameter_units
from ..config.parameters import definition as parameters_definition
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, Configuration

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

        # The ski file template
        self.ski = None

        # The individual configurations
        self.parameters_config = None
        self.ranges_config = None
        self.filters_config = None

        # The final fitting config
        self.fitting_config = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.load_input()

        # 3. Prompt for the fitting parameters
        self.prompt_parameters()

        # 4. Prompt for the physical parameter ranges
        self.prompt_ranges()

        # 5. Prompt for the fitting filters
        self.prompt_filters()

        # 6. Adjust the labels of the template ski file
        self.adjust_labels()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

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

    def prompt_parameters(self):

        """
        This function ...
        :return:
        """

        # Create configuration setter
        setter = InteractiveConfigurationSetter("free parameters", add_logging=False, add_cwd=False)

        # Create config
        self.parameters_config = setter.run(parameters_definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def prompt_ranges(self):

        """
        This function ...
        :return:
        """

        # Create the configuration
        definition = ConfigurationDefinition(write_config=False)

        # Add the options
        for label in self.parameters_config.free_parameters:
            in_units_string = " (in " + parameter_units[label] + ")" if label in parameter_units else ""
            definition.add_optional(label + "_range", types_and_ranges[label][0] + "_range", "range of the " + descriptions[label] + in_units_string, default=types_and_ranges[label][1], convert_default=True)

        # Create configuration setter
        setter = InteractiveConfigurationSetter("free parameter ranges", add_logging=False, add_cwd=False)

        # Create config, get the range for each chosen free parameter
        self.ranges_config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def prompt_filters(self):

        """
        This function ...
        :return:
        """

        # Create the configuration
        definition = ConfigurationDefinition(write_config=False)

        # Choose from all the possible filter names
        definition.add_required("filters", "string_list", "the filters for which to use the observed flux as reference for the fitting procedure", choices=self.observed_filter_names)

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
        self.fitting_config = combine_configs(self.parameters_config, self.ranges_config, self.filters_config)

        # Write the configuration
        self.fitting_config.save(self.fitting_configuration_path)

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
