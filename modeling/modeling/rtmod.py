#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.rtmod Contains the RTMod class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.configurable import InteractiveConfigurable, InvalidCommandError
from ...core.basics.configuration import ConfigurationDefinition
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import filesystem as fs
from ..core.environment import GalaxyModelingEnvironment
from ...core.tools import strings, types
from ...core.tools.formatting import print_dictionary

# Fitting
from ..fitting.manager import GenerationManager
from ..fitting.refitter import Refitter
from ..fitting.expander import ParameterExpander
from ..fitting.statistics import FittingStatistics

# Analysis
from ..analysis.initialization import AnalysisInitializer
from ..analysis.launcher import AnalysisLauncher
from ..analysis.analysis import Analysis

# -----------------------------------------------------------------

# TODO: MERGE THIS CLASS WITH THE GALAXYMODELER CLASS

# -----------------------------------------------------------------

# Standard commands
_help_command_name = "help"
_history_command_name = "history"
_status_command_name = "status"

# -----------------------------------------------------------------

# Show stuff
_info_command_name = "info"
_properties_command_name = "properties"
_commands_command_name = "commands"

# -----------------------------------------------------------------

# GETTING DATA
_get_properties_command_name = "get_properties"
_get_seds_command_name = "get_seds"
_get_images_command_name = "get_images"
_inspect_data_command_name = "inspect_data"

# PREPARATION
_initialize_preparation_command_name = "initialize_preparation"
_inspect_initialization_command_name = "inspect_initialization"
_prepare_command_name = "prepare"
_inspect_preparation_command_name = "inspect_preparation"

# DECOMPOSITION
_decompose_command_name = "decompose"

# TRUNCATION
_truncate_command_name = "truncate"
_set_levels_command_name = "set_levels"

# PHOTOMETRY
_photometry_command_name = "photometry"

# MAP MAKING
_make_colour_maps_command_name = "make_colour_maps"
_make_ssfr_maps_command_name = "make_ssfr_maps"
_make_tir_maps_command_name = "make_tir_maps"
_make_attenuation_maps_command_name = "make_attenuation_maps"
_make_old_maps_command_name = "make_old_maps"
_make_dust_maps_command_name = "make_dust_maps"
_make_young_maps_command_name = "make_young_maps"
_make_ionizing_maps_command_name = "make_ionizing_maps"
_select_maps_command_name = "select_maps"
_make_component_maps_command_name = "make_component_maps"

# MODEL BUILDING
_build_model_command_name = "build_model"
_generate_representations_command_name = "generate_representations"

# FITTING
_configure_fit_command_name = "configure_fit"
_initialize_fit_command_name = "initialize_fit"
_explore_command_name = "explore"
_fit_sed_command_name = "fit_sed"
_finish_command_name = "finish_fit"
_manage_generation_command_name = "manage_generation" ## Working
_refit_command_name = "refit" ## Working
_expand_command_name = "expand" ## Working
_statistics_command_name = "statistics" ## Working
_remove_generations_command_name = "remove_generations" ## Working

# ANALYSIS
_initialize_analysis_command_name = "initialize_analysis"
_launch_analysis_command_name = "launch_analysis" ## Working
_manage_analysis_command_name = "manage_analysis" ## Working
_analyse_command_name = "analyse" ## Working

_clear_command_name = "clear"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()

# Standard commands
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show analysis status", None)

# Galaxy properties
commands[_info_command_name] = ("show_info", False, "show galaxy info", None)
commands[_properties_command_name] = ("show_properties", False, "show galaxy properties", None)
commands[_commands_command_name] = ("show_commands_command", True, "show history of modeling commands", None)

# DATA
# PREPARATION
# DECOMPOSITION
# TRUNCATION
# PHOTOMETRY
# MAP MAKING
# MODEL BUILDING

# FITTING
commands[_manage_generation_command_name] = ("manage_generation_command", True, "manage simulations of a generation", "fitting_run_generation")
commands[_refit_command_name] = ("refit_command", True, "adjust the fitting mechanism and results", "fitting_run")
commands[_expand_command_name] = ("expand_command", True, "expand the parameter space of a generation", "fitting_run_generation")
commands[_statistics_command_name] = ("statistics_command", True, "view statistics of a fitting run", "fitting_run")
commands[_remove_generations_command_name] = ("remove_generations_command", True, "remove certain generation(s)", "fitting_run_generations")

# ANALYSIS
commands[_initialize_analysis_command_name] = ("initialize_analysis_command", True, "create an analysis run", None)
commands[_launch_analysis_command_name] = ("launch_analysis_command", True, "launch the analysis simulations", "analysis_run")
commands[_manage_analysis_command_name] = ("manage_analysis_command", True, "manage the analysis simulations", "analysis_run")
commands[_analyse_command_name] = ("analyse_command", True, "view analysis results", "analysis_run")

# General
commands[_clear_command_name] = (None, False, "clear the output of one of the RT modeling steps", None)

# -----------------------------------------------------------------

_seds_command_name = "seds"
_images_command_name = "images"
_preparation_command_name = "preparation"
_decomposition_command_name = "decomposition"
_truncation_command_name = "truncation"
_maps_command_name = "maps"
_fitting_command_name = "fitting"
_analysis_command_name = "analysis"
_generation_command_name = "generation"
_simulations_command_name = "simulations"

# -----------------------------------------------------------------

# Clear
clear_commands = OrderedDict()

# Modeling steps
clear_commands[_seds_command_name] = ("clear_seds", False, "clear SEDs", None)
clear_commands[_images_command_name] = ("clear_images", False, "clear images", None)
clear_commands[_preparation_command_name] = ("clear_preparation", False, "clear preparation", None)
clear_commands[_decomposition_command_name] = ("clear_decomposition", False, "clear decomposition", None)
clear_commands[_truncation_command_name] = ("clear_truncation", False, "clear truncation", None)
clear_commands[_photometry_command_name] = ("clear_photometry", False, "clear photometry", None)
clear_commands[_maps_command_name] = ("clear_maps", False, "clear_maps", None)
clear_commands[_fitting_command_name] = ("clear_fitting", False, "clear fitting", None)
clear_commands[_analysis_command_name] = ("clear_analysis", False, "clear analysis", None)

# Fitting output
clear_commands[_generation_command_name] = ("clear_generation_command", True, "clear a generation", "fitting_run_generation")
clear_commands[_simulations_command_name] = ("clear_simulations_command", True, "clear simulations", "fitting_run_generation_simulations")

# -----------------------------------------------------------------

# Set subcommands
subcommands = OrderedDict()
subcommands[_clear_command_name] = clear_commands

# -----------------------------------------------------------------

class RTMod(InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _subcommands = subcommands

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RTMod, self).__init__(*args, **kwargs)

        # The modeling environemnt
        self.environment = None

    # -----------------------------------------------------------------

    @property
    def do_commands(self):

        """
        This function ...
        :return:
        """

        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):

        """
        This function ...
        :return:
        """

        return self.config.interactive

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Run commands
        if self.do_commands: self.run_commands()

        # 3. Interactive
        if self.do_interactive: self.interactive()

        # 4. Show
        self.show()

        # 5. Write the history
        if self.has_commands: self.write_history()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RTMod, self).setup(**kwargs)

        # Load the modeling environment
        self.environment = GalaxyModelingEnvironment(self.config.path)

    # -----------------------------------------------------------------

    @property
    def modeling_configuration(self):

        """
        This fin
        :return:
        """

        return self.environment.modeling_configuration

    # -----------------------------------------------------------------

    @property
    def modeling_history(self):

        """
        This function ...
        :return:
        """

        return self.environment.history

    # -----------------------------------------------------------------

    @property
    def modeling_commands(self):

        """
        This function ...
        :return:
        """

        return self.modeling_history.commands

    # -----------------------------------------------------------------

    @property
    def unique_modeling_commands(self):

        """
        This function ...
        :return:
        """

        return self.modeling_history.unique_commands

    # -----------------------------------------------------------------

    @property
    def finished_modeling_commands(self):

        """
        This function ...
        :return:
        """

        return self.modeling_history.finished_commands

    # -----------------------------------------------------------------

    @property
    def modeling_status(self):

        """
        This function ...
        :return:
        """

        return self.environment.status

    # -----------------------------------------------------------------

    @property
    def galaxy_properties(self):

        """
        This function ...
        :return:
        """

        return self.environment.galaxy_properties

    # -----------------------------------------------------------------

    @property
    def galaxy_info(self):

        """
        This function ...
        :return:
        """

        return self.environment.galaxy_info

    # -----------------------------------------------------------------

    @property
    def maps_collection(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_collection

    # -----------------------------------------------------------------

    @property
    def model_suite(self):

        """
        This function ...
        :return:
        """

        return self.environment.model_suite

    # -----------------------------------------------------------------

    @property
    def fitting_context(self):

        """
        This function ...
        :return:
        """

        return self.environment.fitting_context

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):

        """
        This function ...
        :return:
        """

        return self.environment.fitting_runs

    # -----------------------------------------------------------------

    @property
    def fitting_run_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_runs.names

    # -----------------------------------------------------------------

    @memoize_method
    def get_fitting_run(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.fitting_runs.load(name)

    # -----------------------------------------------------------------

    def get_generation_names(self, fitting_run_name):

        """
        This function ...
        :param fitting_run_name:
        :return:
        """

        return self.get_fitting_run(fitting_run_name).generation_names

    # -----------------------------------------------------------------

    @memoize_method
    def get_generation(self, fitting_run_name, name):

        """
        This function ...
        :param fitting_run_name:
        :param name:
        :return:
        """

        return self.get_fitting_run(fitting_run_name).get_generation(name)

    # -----------------------------------------------------------------

    def get_simulation_names(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        return self.get_generation(fitting_run_name, generation_name).simulation_names

    # -----------------------------------------------------------------

    def get_fitting_run_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("fitting_run", "integer_or_string", "fitting run name", choices=self.fitting_run_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_fitting_run_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]

        # Get the definition
        definition = self.get_fitting_run_generation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get fitting run name
        fitting_run_name = config.pop("fitting_run")

        # Return
        return splitted, fitting_run_name, config

    # -----------------------------------------------------------------

    def get_fitting_run_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, config = self.parse_fitting_run_command(command, name=name, interactive=interactive)

        # Return the fitting run name
        return fitting_run_name

    # -----------------------------------------------------------------

    def get_fitting_run_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, config = self.parse_fitting_run_command(command, command_definition, name=name, interactive=interactive)

        # Return the fitting run name and config
        return fitting_run_name, config

    # -----------------------------------------------------------------

    def get_fitting_run_generation_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("fitting_run", "integer_or_string", "fitting run name", choices=self.fitting_run_names)
        definition.add_required("generation", "integer_or_string", "generation name or index")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_fitting_run_generation_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]

        # Get the definition
        definition = self.get_fitting_run_generation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get fitting run name
        fitting_run_name = config.pop("fitting_run")

        # Get generation name
        if types.is_integer_type(config.generation): generation_name = self.get_generation_names(fitting_run_name)[config.pop("generation")]
        else: generation_name = config.pop("generation")

        # Return
        return splitted, fitting_run_name, generation_name, config

    # -----------------------------------------------------------------

    def get_fitting_run_name_generation_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, generation_name, config = self.parse_fitting_run_generation_command(command, name=name, interactive=interactive)

        # Return the fitting run name and generation name
        return fitting_run_name, generation_name

    # -----------------------------------------------------------------

    def get_fitting_run_name_generation_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, generation_name, config = self.parse_fitting_run_generation_command(command, command_definition, name=name, interactive=interactive)

        # Return the fitting run name, generation name, and config
        return fitting_run_name, generation_name, config

    # -----------------------------------------------------------------

    def get_fitting_run_generations_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("fitting_run", "integer_or_string", "fitting run name", choices=self.fitting_run_names)
        definition.add_required("generations", "integer_and_string_list", "generation names or indices")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_fitting_run_generations_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]

        # Get the definition
        definition = self.get_fitting_run_generations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get fitting run name
        fitting_run_name = config.pop("fitting_run")

        # Get generation names
        generation_names = []
        for index_or_name in config.generations:
            if types.is_string_type(index_or_name):
                generation_name = index_or_name
                if generation_name not in self.get_generation_names(fitting_run_name): raise InvalidCommandError("Invalid generation name '" + generation_name + "'", command)
            elif types.is_integer_type(index_or_name): generation_name = self.get_generation_names(fitting_run_name)[index_or_name]
            else: raise ValueError("Invalid type")
            generation_names.append(generation_name)

        # Return
        return splitted, fitting_run_name, generation_names, config

    # -----------------------------------------------------------------

    def get_fitting_run_name_generation_names_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, generation_names, config = self.parse_fitting_run_generations_command(command, name=name, interactive=interactive)

        # Return the fitting run name and generation names
        return fitting_run_name, generation_names

    # -----------------------------------------------------------------

    def get_fitting_run_name_generation_names_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_name, generation_names, config = self.parse_fitting_run_generations_command(command, command_definition, name=name, interactive=interactive)

        # Return the fitting run name, generation names, and config
        return fitting_run_name, generation_names, config

    # -----------------------------------------------------------------

    def get_fitting_run_generation_simulations_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("fitting_run", "integer_or_string", "fitting run name", choices=self.fitting_run_names)
        definition.add_required("generation", "integer_or_string", "generation name or index")
        definition.add_required("simulations", "integer_and_string_list", "simulation names or indices")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_fitting_run_generation_simulations_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]

        # Get the definition
        definition = self.get_fitting_run_generation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get fitting run name
        fitting_run_name = config.pop("fitting_run")

        # Get generation name
        if types.is_integer_type(config.generation): generation_name = self.get_generation_names(fitting_run_name)[config.pop("generation")]
        else: generation_name = config.pop("generation")

        # Get simulation names
        simulation_names = []
        for index_or_name in config.simulations:
            if types.is_string_type(index_or_name):
                simulation_name = index_or_name
                if simulation_name not in self.get_simulation_names(fitting_run_name, generation_name): raise InvalidCommandError("Invalid simulation name '" + simulation_name + "'", command)
            elif types.is_integer_type(index_or_name): simulation_name = self.get_simulation_names(fitting_run_name, generation_name)[index_or_name]
            else: raise ValueError("Invalid type")
            simulation_names.append(simulation_name)

        # Return
        return splitted, fitting_run_name, generation_name, simulation_names, config

    # -----------------------------------------------------------------

    @property
    def analysis_context(self):

        """
        This function ...
        :return:
        """

        return self.environment.analysis_context

    # -----------------------------------------------------------------

    @property
    def analysis_runs(self):

        """
        This function ...
        :return:
        """

        return self.environment.analysis_runs

    # -----------------------------------------------------------------

    @property
    def analysis_run_names(self):

        """
        This function ...
        :return:
        """

        return self.analysis_runs.names

    # -----------------------------------------------------------------

    @memoize_method
    def get_analysis_run(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.analysis_runs.load(name)

    # -----------------------------------------------------------------

    def get_analysis_run_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("analysis_run", "integer_or_string", "analysis run name", choices=self.analysis_run_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_analysis_run_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only analysis run name

        # Get the definition
        definition = self.get_analysis_run_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get the analysis run name
        analysis_run_name = config.pop("analysis_run")

        # Return
        return splitted, analysis_run_name, config

    # -----------------------------------------------------------------

    def get_analysis_run_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, analysis_run_name, config = self.parse_analysis_run_command(command, name=name, interactive=interactive)

        # Return the names
        return analysis_run_name

    # -----------------------------------------------------------------

    def get_analysis_run_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, analysis_run_name, config = self.parse_analysis_run_command(command, command_definition, name=name, interactive=interactive)

        # Return the names
        return analysis_run_name, config

    # -----------------------------------------------------------------

    def show_info(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the galaxy info ...")

        # Show
        print_dictionary(self.galaxy_info)

    # -----------------------------------------------------------------

    def show_properties(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the galaxy properties ...")

        print("")
        self.galaxy_properties.show_properties()
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_commands_definition(self):

        """
        This function ...
        :return:
        """

        #
        definition = ConfigurationDefinition(write_config=False)

        # Add options
        definition.add_flag("finished", "only show finished commands")
        definition.add_optional("clear", "string_list", "commands to remove from the history")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def show_commands_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.show_commands_definition, **kwargs)

        # Clear?
        if config.clear is not None: self.clear_commands(config.clear)

        # Show
        self.show_commands()

    # -----------------------------------------------------------------

    def clear_commands(self, commands):

        """
        This function ...
        :param commands:
        :return:
        """

        # Loop over the commands
        for command in commands: self.history.remove_entries(command)

        # Save
        self.history.save()

    # -----------------------------------------------------------------

    def show_commands(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the history of modeling commands ...")

        # Show
        for item in self.unique_modeling_commands: print(item)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_status_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def show_status_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.show_status_definition, **kwargs)

        # Show
        self.show_status()

    # -----------------------------------------------------------------

    def show_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing status ...")

        # Show
        self.modeling_status.show(name=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_generation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def manage_generation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name and generation name
        fitting_run_name, generation_name, config = self.get_fitting_run_name_generation_name_and_config_from_command(command, self.manage_generation_definition, **kwargs)

        # Manage
        self.manage_generation(fitting_run_name, generation_name)

    # -----------------------------------------------------------------

    def manage_generation(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        manager = GenerationManager()

        manager.run = fitting_run_name
        manager.generation = generation_name

        manager.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def refit_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def refit_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name
        fitting_run_name, config = self.get_fitting_run_name_and_config_from_command(command, self.refit_definition, **kwargs)

        # Refit
        self.refit(fitting_run_name)

    # -----------------------------------------------------------------

    def refit(self, fitting_run_name):

        """
        This function ...
        :param fitting_run_name:
        :return:
        """

        refitter = Refitter()

        refitter.fitting_run = fitting_run_name

        refitter.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def expand_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def expand_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name and generation name
        fitting_run_name, generation_name, config = self.get_fitting_run_name_generation_name_and_config_from_command(command, self.expand_definition, **kwargs)

        # Expand
        self.expand(fitting_run_name, generation_name)

    # -----------------------------------------------------------------

    def expand(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        expander = ParameterExpander()

        expander.config.run = fitting_run_name
        expander.config.generation = generation_name

        expander.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def statistics_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def statistics_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name
        fitting_run_name = self.get_fitting_run_name_and_config_from_command(command, self.statistics_definition, **kwargs)

        # Statistics
        self.statistics(fitting_run_name)

    # -----------------------------------------------------------------

    def statistics(self, fitting_run_name):

        """
        This function ...
        :param fitting_run_name:
        :return:
        """

        statistics = FittingStatistics()

        statistics.run = fitting_run_name

        statistics.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def remove_generations_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def remove_generations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name and generation names
        fitting_run_name, generation_names, config = self.get_fitting_run_name_generation_names_and_config_from_command(command, self.remove_generations_definition, **kwargs)

        # Remove
        self.remove_generations(fitting_run_name, generation_names)

    # -----------------------------------------------------------------

    def remove_generations(self, fitting_run_name, generation_names):

        """
        This function ...
        :param fitting_run_name:
        :param generation_names:
        :return:
        """

        # Remove
        for generation_name in generation_names: self.remove_generation(fitting_run_name, generation_name)

    # -----------------------------------------------------------------

    def remove_generation(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        # Inform the user
        log.info("Removing generation '" + generation_name + "' ...")

        # Get the fitting run
        fitting_run = self.get_fitting_run(fitting_run_name)

        # Remove this generation from the generations table
        fitting_run.generations_table.remove_entry(generation_name)
        fitting_run.generations_table.save()

        # Get generation path
        generation_path = fitting_run.get_generation_path(generation_name)

        # Remove the generation directory
        fs.remove_directory(generation_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def initialize_analysis_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def initialize_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.initialize_analysis_definition, **kwargs)

        # Initialize
        self.initialize_analysis()

    # -----------------------------------------------------------------

    def initialize_analysis(self):

        """
        This function ...
        :return:
        """

        initializer = AnalysisInitializer()

        initializer.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def launch_analysis_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def launch_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get analysis run name
        analysis_run_name, config = self.get_analysis_run_name_and_config_from_command(command, self.launch_analysis_definition, **kwargs)

        # Launch
        self.launch_analysis(analysis_run_name)

    # -----------------------------------------------------------------

    def launch_analysis(self, analysis_run_name):

        """
        This function ...
        :param analysis_run_name:
        :return:
        """

        launcher = AnalysisLauncher()

        launcher.config.run = analysis_run_name

        launcher.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_analysis_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def manage_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get analysis run name
        analysis_run_name, config = self.get_analysis_run_name_and_config_from_command(command, self.manage_analysis_definition, **kwargs)

        # Manage
        self.manage_analysis(analysis_run_name)

    # -----------------------------------------------------------------

    def manage_analysis(self, analysis_run_name):

        """
        This function ...
        :param analysis_run_name:
        :return:
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def analyse_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the analysis run name
        analysis_run_name, config = self.get_analysis_run_name_and_config_from_command(command, self.analyse_definition, **kwargs)

        # Analyse
        self.analyse(analysis_run_name)

    # -----------------------------------------------------------------

    def analyse(self, analysis_run_name):

        """
        This function ...
        :param analysis_run_name:
        :return:
        """

        analysis = Analysis()

        analysis.config.run = analysis_run_name

        analysis.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add setting
        definition.add_required("step", "string", "the modeling step for which to clear the output")

        # Return
        return definition

    # -----------------------------------------------------------------

    def clear_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the SEDs ...")

        # Clear
        fs.clear_directory(self.environment.data_seds_path)

    # -----------------------------------------------------------------

    def clear_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the images ...")

        # Clear
        fs.clear_directory(self.environment.data_images_path)

    # -----------------------------------------------------------------

    def clear_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the preparation output ...")

        # Clear
        fs.clear_directory(self.environment.prep_path)

    # -----------------------------------------------------------------

    def clear_decomposition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the decomposition output ...")

        # Clear
        fs.clear_directory(self.environment.components_path)

    # -----------------------------------------------------------------

    def clear_truncation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the truncation output ...")

        # Clear
        fs.clear_directory(self.environment.truncation_path)

    # -----------------------------------------------------------------

    def clear_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the photometry output ...")

        # Clear
        fs.clear_directory(self.environment.phot_path)

    # -----------------------------------------------------------------

    def clear_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the maps ...")

        # Clear
        fs.clear_directory(self.environment.maps_path)

    # -----------------------------------------------------------------

    def clear_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the fitting output ...")

        # Clear
        fs.clear_directory(self.environment.fit_path)

    # -----------------------------------------------------------------

    def clear_analysis(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the analysis output ...")

        # Clear
        fs.clear_directory(self.environment.analysis_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_generation_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)

        # Flags
        definition.add_flag("backup", "make backups of non-empty directories", False)
        definition.add_flag("adapt_simulations", "unset retrieved and analysed flags of the corresponding simulations",
                            False)

        # Generations to remove
        definition.add_required("generation", "string", "generation to remove (none means all)")

        # Get configuration
        #config = parse_arguments("clear_generation_output", definition,
        #                         "Clear the output of all simulations of a certain generation")

        return definition

    # -----------------------------------------------------------------

    def clear_generation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def clear_generation(self, generation_name):

        """
        This function ...
        :return:
        """

        # Load the fitting run
        fitting_run = runs.load(config.name)

        # Check
        if not fitting_run.is_generation(config.generation): raise ValueError("Generation doesn't exist")

        # Get generation path
        generation_path = fitting_run.get_generation_path(config.generation)

        # Loop over the simulation paths
        for path, name in fs.directories_in_path(generation_path, returns=["path", "name"]):

            # Determine subpaths
            out_path = fs.join(path, "out")
            extr_path = fs.join(path, "extr")
            plot_path = fs.join(path, "plot")
            misc_path = fs.join(path, "misc")

            # Non-empty output
            if fs.is_directory(out_path) and not fs.is_empty(out_path):

                # Debugging
                log.debug("Clearing output directory of '" + name + "' simulation ...")
                if config.backup: fs.backup_directory(out_path)
                fs.clear_directory(out_path)

            # Non-empty extracted data
            if fs.is_directory(extr_path) and not fs.is_empty(extr_path):

                # Debugging
                log.debug("Clearing extraction directory of '" + name + "' simulation ...")
                if config.backup: fs.backup_directory(extr_path)
                fs.clear_directory(extr_path)

            # Non-empty plotting output
            if fs.is_directory(plot_path) and not fs.is_empty(plot_path):

                # Debugging
                log.debug("Clearing plotting directory of '" + name + "' simulation ...")
                if config.backup: fs.backup_directory(plot_path)
                fs.clear_directory(plot_path)

            # Non-empty misc output
            if fs.is_directory(misc_path) and not fs.is_empty(misc_path):

                # Debugging
                log.debug("Clearing miscellaneous output directory of '" + name + "' simulation ...")
                if config.backup: fs.backup_directory(misc_path)
                fs.clear_directory(misc_path)

        # -----------------------------------------------------------------

        # Clear chi-squared table
        chi_squared_path = fitting_run.chi_squared_table_path_for_generation(config.generation)
        chi_squared = fitting_run.chi_squared_table_for_generation(config.generation)

        # Non-empty
        if len(chi_squared) > 0:

            log.debug("Clearing chi-squared table ...")
            if config.backup: fs.backup_file(chi_squared_path)
            chi_squared.remove_all_rows()
            chi_squared.save()

        # -----------------------------------------------------------------

        # Adapt simulations?
        if config.adapt_simulations:

            # Get the simulations for the generation
            generation = fitting_run.get_generation(config.generation)

            # Check if has assignment table
            if generation.has_assignment_table:

                # Get the simulations
                for simulation in generation.simulations:

                    # Set flag
                    changed = False

                    # Unset retrieved
                    changed = simulation.set_retrieved(False)

                    # Unset analysed
                    # print(simulation.retrieved, simulation.analysed, simulation.analysed_extraction, simulation.analysed_plotting, simulation.analysed_misc, simulation.analysed_batch, simulation.analysed_scaling, simulation.analysed_extra)
                    changed |= simulation.set_analysed(False, all=True)

                    # Has changed
                    if changed:
                        log.debug("Adapting retrieved and analysed flags for simulation '" + simulation.name + "' ...")
                        simulation.save()

            # Give warning
            else: log.warning("Assignment table not found")

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_simulations_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)

        # Generations to remove
        definition.add_positional_optional("generations", "string_list",
                                           "generation(s) for which to remove the simulations")
        definition.add_optional("remote", "string", "remote host for which to remove the simulations (none means all")

        # Remove simulations of generations other than certain generations
        definition.add_optional("other_generations", "string_list",
                                "remove all other generations than these generations")
        definition.add_flag("hard",
                            "remove ALL simulations objects that are not from the 'other_generations' (use with CARE!)")

        # Make backup?
        definition.add_flag("backup", "make backup of simulation files")
        definition.add_optional("backup_path", "directory_path", "backup directory path")

        # Get configuration
        #config = parse_arguments("clear_generation_simulations", definition,
        #                         "Remove the simulation files for a certain generation")

        return definition

    # -----------------------------------------------------------------

    def clear_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

    # -----------------------------------------------------------------

    def clear_simulations(self, generation_name, simulation_names):

        """
        This function ...
        :return:
        """

        # Load the fitting run
        fitting_run = runs.load(config.name)

        # HARD delete
        if config.hard:

            # Initialize a dictionary with the simulation filepaths to keep per host
            keep_paths = defaultdict(list)

            # Loop over the generations of which we CANNOT remove the simulations
            for generation_name in config.other_generations:

                # Check
                if not fitting_run.is_generation(generation_name): raise ValueError(
                    "Generation '" + generation_name + "' doesn't exist")

                # Get generation path
                generation_path = fitting_run.get_generation_path(generation_name)

                # Get the generation
                generation = fitting_run.get_generation(generation_name)

                # Loop over the remote hosts
                for host_id in generation.host_ids:
                    # Get paths of the simulation files
                    filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore",
                                                                         as_dict=False)
                    nfilepaths = len(filepaths)

                    # Add the filepaths
                    keep_paths[host_id].extend(filepaths)

            # Loop over the remote hosts
            if config.remote is not None:
                host_ids = [config.remote]
            else:
                host_ids = keep_paths.keys()

            # Loop over the hosts
            for host_id in host_ids:

                # Get simulation paths
                paths = get_simulation_paths_for_host(host_id)

                # keep = keep_paths[host_id]
                # print(keep)

                # Get list of filepaths to remove
                remove_paths = sequences.get_other(paths, keep_paths[host_id])
                remove_names = [fs.strip_extension(fs.name(filepath)) for filepath in remove_paths]

                # Proceed?
                if not prompt_proceed("proceed removing the simulation files: " + tostr(remove_names)): continue

                # Backup
                if config.backup:
                    log.debug("Creating backup of the simulation files ...")
                    if config.remote is not None:
                        backup_path = config.backup_path
                    else:
                        backup_path = fs.create_directory_in(config.backup_path, host_id)
                    fs.copy_files(remove_paths, backup_path)

                # Remove
                fs.remove_files(remove_paths)

        # Look only at other generations of the same fitting run
        else:

            # Determine generation names
            if config.generations is not None:

                if config.other_generations is not None: raise ValueError("Cannot specify also 'other_generations'")
                generation_names = config.generations

            # Other generations
            elif config.other_generations is not None:

                # Get the other names than those specified
                generation_names = sequences.get_other(fitting_run.generation_names, config.other_generations)

            # Which generations?
            else:
                raise ValueError("Generation names cannot be determined")

            # Show
            log.debug("Removing the simulation files for generations: '" + tostr(generation_names) + "' ...")

            # Loop over the generations
            for generation_name in generation_names:

                # Check
                if not fitting_run.is_generation(generation_name): raise ValueError(
                    "Generation '" + generation_name + "' doesn't exist")

                # Get generation path
                generation_path = fitting_run.get_generation_path(generation_name)

                # Get the generation
                generation = fitting_run.get_generation(generation_name)

                # Loop over the remote hosts
                for host_id in generation.host_ids:

                    # Clear for this host?
                    if config.remote is not None and host_id != config.remote: continue

                    # Get paths of the simulation files
                    filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore")
                    nfilepaths = len(filepaths)
                    if nfilepaths == 0:
                        log.warning("No simulation objects anymore for generation '" + generation_name + "'")
                        continue  # skip generation

                    # Proceed?
                    remove_names = [fs.strip_extension(fs.name(filepath)) for filepath in filepaths.keys()]
                    if not prompt_proceed(
                        "proceed removing the simulation files of generation '" + generation_name + "': " + tostr(
                            remove_names)): continue

                    # Get the filepaths
                    paths = filepaths.values()

                    # Backup
                    if config.backup:
                        log.debug("Creating backup of the simulation files of generation '" + generation_name + "' ...")
                        fs.copy_files(paths, config.backup_path)

                    # Remove
                    fs.remove_files(paths)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "rtmod"

# -----------------------------------------------------------------
