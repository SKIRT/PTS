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
from collections import OrderedDict, defaultdict
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ...core.basics.configurable import InteractiveConfigurable, InvalidCommandError
from ...core.basics.configuration import ConfigurationDefinition, prompt_proceed
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import filesystem as fs
from ..core.environment import GalaxyModelingEnvironment
from ...core.tools import strings, types, sequences, time
from ...core.tools.formatting import print_dictionary
from ...core.tools.stringify import tostr
from ...core.simulation.remote import get_simulation_paths_for_host
from ...core.tools import formatting as fmt
from ...core.tools.serialization import load_dict, write_dict
from ...core.launch.manager import extra_columns
from ...core.remote.host import find_host_ids
from ..analysis.run import AnalysisRunInfo
from ..config.expand import definition as expand_definition
from ..config.manage_generation import definition as manage_generation_definition
from ..config.refit import definition as refit_definition

# Fitting
from ..fitting.manager import GenerationManager
from ..fitting.refitter import Refitter
from ..fitting.expander import ParameterExpander
from ..fitting.statistics import FittingStatistics
from ..fitting.refitter import clone_fitting_run

# Analysis
from ..analysis.initialization import AnalysisInitializer
from ..analysis.launcher import AnalysisLauncher
from ..analysis.analysis import Analysis
from ..analysis.manager import AnalysisManager

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
_component_maps_command_name = "component_maps"

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
_generation_output_command_name = "generation_output" ## Working
_generation_status_command_name = "generation_status" ## Working
_generations_command_name = "generations" ## Working

# ANALYSIS
_initialize_analysis_command_name = "initialize_analysis"
_launch_analysis_command_name = "launch_analysis" ## Working
_manage_analysis_command_name = "manage_analysis" ## Working
_analyse_command_name = "analyse" ## Working

# With subcommands
_clear_command_name = "clear"
_clone_command_name = "clone"

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
clear_commands.description = "clear the output of one of the RT modeling steps"

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
clear_commands[_generation_command_name] = ("clear_generation_command", True, "clear the output of all simulations of a certain generation", "fitting_run_generation")
clear_commands[_simulations_command_name] = ("clear_simulations_command", True, "remove the simulation files for a certain generation", "fitting_run_generation_simulations")

# -----------------------------------------------------------------

# Clone
clone_commands = OrderedDict()
clone_commands.description = "clone a fitting run or analysis run"

# Fitting run and analysis run
clone_commands[_fitting_command_name] = ("clone_fitting_command", True, "clone a fitting run (without the generations)", "fitting_run")
clone_commands[_analysis_command_name] = ("clone_analysis_command", True, "clone an analysis run (without the results)", "analysis_run")

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
commands[_component_maps_command_name] = ("show_component_maps", False, "list the component maps", None)

# MODEL BUILDING

# FITTING
commands[_manage_generation_command_name] = ("manage_generation_command", True, "manage simulations of a generation", "fitting_run_generation")
commands[_refit_command_name] = ("refit_command", True, "adjust the fitting mechanism and results", "fitting_run")
commands[_expand_command_name] = ("expand_command", True, "expand the parameter space of a generation", "fitting_run_generation")
commands[_statistics_command_name] = ("statistics_command", True, "view statistics of a fitting run", "fitting_run")
commands[_remove_generations_command_name] = ("remove_generations_command", True, "remove certain generation(s)", "fitting_run_generations")
commands[_generation_output_command_name] = ("show_generation_output_command", True, "show the output of a generation", "fitting_run_generation")
commands[_generation_status_command_name] = ("show_generation_status_command", True, "show the status of a generation", "fitting_run_generation")
commands[_generations_command_name] = ("show_generations_command", True, "show the generations", "fitting_runs")

# ANALYSIS
commands[_initialize_analysis_command_name] = ("initialize_analysis_command", True, "create an analysis run", None)
commands[_launch_analysis_command_name] = ("launch_analysis_command", True, "launch the analysis simulations", "analysis_run")
commands[_manage_analysis_command_name] = ("manage_analysis_command", True, "manage the analysis simulations", "analysis_run")
commands[_analyse_command_name] = ("analyse_command", True, "view analysis results", "analysis_run")

# General
commands[_clear_command_name] = clear_commands # (None, False, "clear the output of one of the RT modeling steps", None)
commands[_clone_command_name] = clone_commands

# -----------------------------------------------------------------

class RTMod(InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _log_section = "RTMOD"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RTMod, self).__init__(*args, **kwargs)

        # The modeling environment
        self.environment = None

    # -----------------------------------------------------------------

    @property
    def do_commands(self):
        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):
        return self.config.interactive

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Set Matplotlib's backend
        if self.config.mpl_backend is not None: plt.switch_backend(self.config.mpl_backend)

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):
        return self.environment.galaxy_name

    # -----------------------------------------------------------------

    @property
    def modeling_configuration(self):
        return self.environment.modeling_configuration

    # -----------------------------------------------------------------

    @property
    def modeling_history(self):
        return self.environment.history

    # -----------------------------------------------------------------

    @property
    def modeling_commands(self):
        return self.modeling_history.commands

    # -----------------------------------------------------------------

    @property
    def unique_modeling_commands(self):
        return self.modeling_history.unique_commands

    # -----------------------------------------------------------------

    @property
    def finished_modeling_commands(self):
        return self.modeling_history.finished_commands

    # -----------------------------------------------------------------

    @property
    def modeling_status(self):
        return self.environment.status

    # -----------------------------------------------------------------

    @property
    def galaxy_properties(self):
        return self.environment.galaxy_properties

    # -----------------------------------------------------------------

    @property
    def galaxy_info(self):
        return self.environment.galaxy_info

    # -----------------------------------------------------------------

    @property
    def maps_collection(self):
        return self.environment.maps_collection

    # -----------------------------------------------------------------

    @property
    def static_maps_collection(self):
        return self.environment.static_maps_collection

    # -----------------------------------------------------------------

    @property
    def maps_selection(self):
        return self.environment.maps_selection

    # -----------------------------------------------------------------

    @property
    def static_maps_selection(self):
        return self.environment.static_maps_selection

    # -----------------------------------------------------------------

    @property
    def model_suite(self):
        return self.environment.model_suite

    # -----------------------------------------------------------------

    @property
    def fitting_context(self):
        return self.environment.fitting_context

    # -----------------------------------------------------------------

    @property
    def fitting_runs(self):
        return self.environment.fitting_runs

    # -----------------------------------------------------------------

    @property
    def fitting_run_names(self):
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

    def is_generation(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        return generation_name in self.get_generation_names(fitting_run_name)

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

    @memoize_method
    def get_generation_path(self, fitting_run_name, name):

        """
        This function ...
        :param fitting_run_name:
        :param name:
        :return:
        """

        return self.get_fitting_run(fitting_run_name).get_generation_path(name)

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
        definition.add_required("fitting_run", "string", "fitting run name", choices=self.fitting_run_names)

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
        definition = self.get_fitting_run_command_definition(command_definition, required_to_optional=required_to_optional)

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

    def get_fitting_runs_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("fitting_runs", "string_list", "fitting run names", choices=self.fitting_run_names)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_fitting_runs_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

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
        definition = self.get_fitting_runs_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get the configuration
        config = self.get_config_from_definition(name, definition, parse_command, interactive=interactive)

        # Get fitting run name
        fitting_run_names = config.pop("fitting_runs")

        # Return
        return splitted, fitting_run_names, config

    # -----------------------------------------------------------------

    def get_fitting_run_names_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_names, config = self.parse_fitting_runs_command(command, name=name, interactive=interactive)

        # Return the fitting run names
        return fitting_run_names

    # -----------------------------------------------------------------

    def get_fitting_run_names_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, fitting_run_names, config = self.parse_fitting_runs_command(command, command_definition, name=name, interactive=interactive)

        # Return the fitting run names and config
        return fitting_run_names, config

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

        # Check whether generation exists
        if not self.is_generation(fitting_run_name, generation_name): raise InvalidCommandError("Generation '" + generation_name + "' does not exist", command)

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

    def get_fitting_run_name_generation_name_simulation_names_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse
        splitted, fitting_run_name, generation_name, simulation_names, config = self.parse_fitting_run_generation_simulations_command(command, name=name, interactive=interactive)

        # Return
        return fitting_run_name, generation_name, simulation_names

    # -----------------------------------------------------------------

    def get_fitting_run_name_generation_name_simulation_names_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse
        splitted, fitting_run_name, generation_name, simulation_names, config = self.parse_fitting_run_generation_simulations_command(command, command_definition, name=name, interactive=interactive)

        # Return
        return fitting_run_name, generation_name, simulation_names, config

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

    def show_component_maps(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # Old
        self.show_old_component_maps()

        # Young
        self.show_young_component_maps()

        # Ionizing
        self.show_ionizing_component_maps()

        # Dust
        self.show_dust_component_maps()

    # -----------------------------------------------------------------

    def show_old_component_maps(self):

        """
        This function ...
        :return:
        """

        # Old
        print(fmt.blue + "OLD STELLAR DISK" + fmt.reset)
        print("")

        # Loop over the maps
        for name in self.static_maps_selection.old_map_names:

            print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
            print("")

            # Get path
            path = self.static_maps_selection.old_map_paths[name]

            # Get origins
            origins = self.static_maps_collection.get_old_stellar_disk_origins()[name]

            # Get steps path
            steps_path = self.static_maps_selection.get_old_steps_path_for_map(name)

            # Show the info
            show_component_map_info(name, path, origins, steps_path)
            print("")

    # -----------------------------------------------------------------

    def show_young_component_maps(self):

        """
        This function ...
        :return:
        """

        # Young
        print(fmt.blue + "YOUNG STELLAR DISK" + fmt.reset)
        print("")

        # Loop over the maps
        for name in self.static_maps_selection.young_map_names:

            print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
            print("")

            # Get path
            path = self.static_maps_selection.young_map_paths[name]

            # Get origins
            origins = self.static_maps_collection.young_origins_flat[name]

            # Get steps path
            steps_path = self.static_maps_selection.get_young_steps_path_for_map(name)

            # Show the info
            show_component_map_info(name, path, origins, steps_path)
            print("")

    # -----------------------------------------------------------------

    def show_ionizing_component_maps(self):

        """
        This function ...
        :return:
        """

        # Ionizing
        print(fmt.blue + "IONIZING STELLAR DISK" + fmt.reset)
        print("")

        # Loop over the maps
        for name in self.static_maps_selection.ionizing_map_names:

            print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
            print("")

            # Get path
            path = self.static_maps_selection.ionizing_map_paths[name]

            # Get origins
            origins = self.static_maps_collection.ionizing_origins_flat[name]

            # Get steps path
            steps_path = self.static_maps_selection.get_ionizing_steps_path_for_map(name)

            # Show the info
            show_component_map_info(name, path, origins, steps_path)
            print("")

    # -----------------------------------------------------------------

    def show_dust_component_maps(self):

        """
        This function ...
        :return:
        """

        # Dust
        print(fmt.blue + "DUST DISK" + fmt.reset)
        print("")

        # Loop over the maps
        for name in self.static_maps_selection.dust_map_names:

            print(" - " + fmt.green + fmt.underlined + name + fmt.reset)
            print("")

            # Get path
            path = self.static_maps_selection.dust_map_paths[name]

            # Get origins
            origins = self.static_maps_collection.dust_origins_flat[name]

            # Get steps path
            steps_path = self.static_maps_selection.get_dust_steps_path_for_map(name)

            # Show the info
            show_component_map_info(name, path, origins, steps_path)
            print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_generation_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # From manage_generation
        import_names = ["lazy", "find_simulations", "find_remotes", "produce_missing", "check_paths", "correct_paths", "confirm_correction", "check_analysis"]

        # Below: options from manage_simulations that we want the user to be able to control here in this context
        import_names.append("offline")
        import_names.append("dry")
        import_names.append("fix_success")

        import_names.append("backup_simulations")
        import_names.append("backup_assignment")
        import_names.append("backup_path")
        import_names.append("backup_dir_path")
        import_names.append("backup_dirname")

        import_names.append("local")
        import_names.append("warn_local")
        import_names.append("success")

        import_names.append("cache_path")
        import_names.append("cache_root")
        import_names.append("cache_output")
        import_names.append("cache_datacubes")
        import_names.append("cache_misc")
        import_names.append("cache_images")
        import_names.append("cache_after_analysis")

        # Import properties from the manage_generation definition
        definition.import_properties(manage_generation_definition, import_names)

        # Change default
        definition.flags["fix_success"].default = True

        # Write the simulation status to the generation directory
        definition.add_flag("write_status", "write the status table", True)

        # Use previous status
        definition.add_flag("previous_status", "show the previous simulation status", False)

        # Check status
        definition.add_flag("check_status", "check the simulation status", True)

        # Return the definition
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

        # Use previous status?
        if config.previous_status:
            config.pop("previous_status")
            status_filepath = fs.join(self.get_generation_path(fitting_run_name, generation_name), "status.dat")
            if not fs.is_file(status_filepath): status_filepath = None # does not exist
        else: status_filepath = None

        # Manage
        self.manage_generation(fitting_run_name, generation_name, config=config, status_filepath=status_filepath)

    # -----------------------------------------------------------------

    def manage_generation(self, fitting_run_name, generation_name, config=None, status_filepath=None):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :param config:
        :param status_filepath:
        :return:
        """

        # Create the manager
        manager = GenerationManager(config=config, extend_config=True, check_required=False)
        manager.config.path = self.config.path

        # Set fitting run name and generation name
        manager.config.run = fitting_run_name
        manager.config.generation = generation_name

        # Set cache volume name
        manager.config.cache_volume = self.config.cache_volume

        # Set status filepath
        manager.config.status = status_filepath

        # Run the manager
        manager.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def refit_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(refit_definition)
        definition.remove_setting("run")
        #definition.remove_optional("generations")

        # Return the definition
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
        self.refit(fitting_run_name, config=config)

    # -----------------------------------------------------------------

    def refit(self, fitting_run_name, config=None):

        """
        This function ...
        :param fitting_run_name:
        :param config:
        :return:
        """

        # Initialize the refitter
        refitter = Refitter(config=config)

        # Set the modeling path
        refitter.config.path = self.config.path

        # Set the fitting run name
        refitter.config.run = fitting_run_name

        # Run the refitting
        refitter.run()

    # -----------------------------------------------------------------

    @lazyproperty
    def expand_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.import_settings(expand_definition)
        definition.remove_setting("run")
        definition.remove_required("generation")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def expand_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def expand_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get all kwargs
        kwargs.update(self.expand_kwargs)

        # Get fitting run name and generation name
        #fitting_run_name, generation_name, config = self.get_fitting_run_name_generation_name_and_config_from_command(command, self.expand_definition, **kwargs)
        splitted, fitting_run_name, generation_name, config = self.parse_fitting_run_generation_command(command, self.expand_definition, **kwargs)

        # Expand
        self.expand(fitting_run_name, generation_name, config)

    # -----------------------------------------------------------------

    def expand(self, fitting_run_name, generation_name, config):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :param config:
        :return:
        """

        # Create expander
        expander = ParameterExpander(config)
        expander.config.path = self.config.path

        # Set
        expander.config.run = fitting_run_name
        expander.config.generation = generation_name

        # Run
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
        fitting_run_name, config = self.get_fitting_run_name_and_config_from_command(command, self.statistics_definition, **kwargs)

        # Statistics
        self.statistics(fitting_run_name)

    # -----------------------------------------------------------------

    def statistics(self, fitting_run_name):

        """
        This function ...
        :param fitting_run_name:
        :return:
        """

        # Initialize
        statistics = FittingStatistics()
        statistics.config.path = self.config.path

        # Set fitting run name
        statistics.config.run = fitting_run_name

        # Run
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

    def show_generation_output_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation name
        fitting_run_name, generation_name = self.get_fitting_run_name_generation_name_from_command(command, **kwargs)

        # Show
        self.show_generation_output(fitting_run_name, generation_name)

    # -----------------------------------------------------------------

    def show_generation_output(self, fitting_run_name, generation_name):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :return:
        """

        # Inform the user
        log.info("Showing generation output ...")

        # Get the generation
        generation = self.get_generation(fitting_run_name, generation_name)

        print("")
        print(fmt.yellow + "FILES" + fmt.reset)
        print("")

        # Show files
        for filename in fs.files_in_path(generation.path, returns="name", extensions=True): print(" - " + fmt.bold + filename + fmt.reset)
        print("")

        print(fmt.yellow + "SIMULATIONS" + fmt.reset)
        print("")

        # Loop over the simulation names
        for simulation_name in generation.simulation_names:

            # Check whether directory is present
            if not fs.contains_directory(generation.path, simulation_name): print(" - " + fmt.red + simulation_name + ": directory missing")
            else:

                simulation_path = fs.join(generation.path, simulation_name)
                ski_path = fs.join(simulation_path, self.galaxy_name + ".ski")
                out_path = fs.join(simulation_path, "out")
                extr_path = fs.join(simulation_path, "extr")
                plot_path = fs.join(simulation_path, "plot")
                misc_path = fs.join(simulation_path, "misc")

                # Check presence
                has_ski = fs.is_file(ski_path)
                has_out = fs.is_directory(out_path)
                has_extr = fs.is_directory(extr_path)
                has_plot = fs.is_directory(plot_path)
                has_misc = fs.is_directory(misc_path)

                missing = []
                if not has_ski: missing.append("ski file")
                if not has_out: missing.append("output directory")
                if not has_extr: missing.append("extraction directory")
                if not has_plot: missing.append("plotting directory")
                if not has_misc: missing.append("misc directory")
                nmissing = len(missing)

                if nmissing == 0: print(" - " + fmt.green + simulation_name)
                else: print(" - " + fmt.red + simulation_name + ": " + ",".join(missing) + " missing")
        print("")

        # Define necessary output files
        logfiles_name = "logfiles"
        seds_name = "seds"
        datacubes_name = "datacubes"

        print(fmt.yellow + "OUTPUT" + fmt.reset)
        print("")

        # Initialize
        missing_output = OrderedDict()
        missing_output[logfiles_name] = []
        missing_output[seds_name] = []
        if generation.use_images: missing_output[datacubes_name] = []

        # Loop over the simulations
        for simulation_name in generation.simulation_names:

            # Get output
            output = generation.get_simulation_output(simulation_name)

            # Check missing files
            if not output.has_logfiles: missing_output[logfiles_name].append(simulation_name)
            if not output.has_seds: missing_output[seds_name].append(simulation_name)
            if generation.use_images and not output.has_total_images: missing_output[datacubes_name].append(simulation_name)

        # Show
        for output_type in missing_output:

            simulation_names = missing_output[output_type]
            nsimulations = len(simulation_names)

            if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
            else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
        print("")

        # Define necessary output files
        timeline_name = "timeline"
        memory_name = "memory"

        print(fmt.yellow + "EXTRACTION" + fmt.reset)
        print("")

        # Initialize
        missing_extraction = OrderedDict()
        missing_extraction[timeline_name] = []
        missing_extraction[memory_name] = []

        # Loop over the simulations
        for simulation_name in generation.simulation_names:

            # Get extraction output
            extraction = generation.get_extraction_output(simulation_name)

            # Check missing files
            if not extraction.has_timeline: missing_extraction[timeline_name].append(simulation_name)
            if not extraction.has_memory: missing_extraction[memory_name].append(simulation_name)

        # Show
        for output_type in missing_extraction:

            simulation_names = missing_extraction[output_type]
            nsimulations = len(simulation_names)

            if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
            else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
        print("")

        print(fmt.yellow + "PLOTTING" + fmt.reset)
        print("")

        # Initialize
        missing_plotting = OrderedDict()
        missing_plotting[seds_name] = []

        # Loop over the simulations
        for simulation_name in generation.simulation_names:

            # Get plotting output
            plotting = generation.get_plotting_output(simulation_name)

            # Check missing files
            if not plotting.has_seds: missing_plotting[seds_name].append(simulation_name)

        # Show
        for output_type in missing_plotting:

            simulation_names = missing_plotting[output_type]
            nsimulations = len(simulation_names)

            if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
            else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
        print("")

        # Define necessary output files and directories
        differences_name = "differences"
        fluxes_name = "fluxes"
        images_name = "images"

        print(fmt.yellow + "MISC" + fmt.reset)
        print("")

        # Initialize
        missing_misc = OrderedDict()
        missing_misc[differences_name] = []
        missing_misc[fluxes_name] = []
        if generation.use_images: missing_misc[images_name] = []

        # Loop over the simulations
        for simulation_name in generation.simulation_names:

            # Get misc output
            misc = generation.get_misc_output(simulation_name)

            # Check missing files and directories
            if generation.use_images:
                if not misc.has_image_fluxes: missing_misc[fluxes_name].append(simulation_name)
            else:
                if not misc.has_fluxes: missing_misc[fluxes_name].append(simulation_name)
            if generation.use_images and not misc.has_images_for_fluxes: missing_misc[images_name].append(simulation_name)
            if not misc.has_other_filename("differences.dat"): missing_misc[differences_name].append(simulation_name)

        # Show
        for output_type in missing_misc:

            simulation_names = missing_misc[output_type]
            nsimulations = len(simulation_names)

            if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
            else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_generation_status_definition(self):

        """
        This function ...
        :return:
        """

        # Hosts
        all_host_ids = find_host_ids()

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Show parameters
        definition.add_optional("extra", "string_list", "show extra info", choices=extra_columns)

        # Options for the manager
        definition.add_flag("offline", "offline mode")
        definition.add_flag("lazy", "lazy mode")
        definition.add_flag("find_simulations", "find missing simulations by searching on simulation name", True)
        definition.add_optional("find_remotes", "string_list", "find missing simulations in these remote hosts", default=all_host_ids, choices=all_host_ids)
        definition.add_flag("produce_missing", "produce missing simulation files", False)
        definition.add_flag("check_paths", "check simulation paths", False)
        definition.add_flag("correct_paths", "correct simulation paths instead of raising errors", False)
        definition.add_flag("confirm_correction", "confirm before correcting paths", False)
        definition.add_flag("fix_success", "check success flags in assignment table", True)
        definition.add_flag("check_analysis", "check analysis output", False)

        # Write the status table to the generation directory
        definition.add_flag("write_status", "write the status", True)

        # Correct analysis status
        definition.add_flag("correct_status", "correct the status", True)

        # Use previous
        definition.add_flag("previous", "show previous status", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def show_generation_status_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get generation and config
        fitting_run_name, generation_name, config = self.get_fitting_run_name_generation_name_and_config_from_command(command, self.show_generation_status_definition, **kwargs)

        # Use previous status?
        if config.previous:
            status_filepath = fs.join(self.get_generation_path(fitting_run_name, generation_name), "status.dat")
            if not fs.is_file(status_filepath): status_filepath = None  # does not exist
        else: status_filepath = None

        # Show
        self.show_generation_status(fitting_run_name, generation_name, config, status_filepath=status_filepath)

    # -----------------------------------------------------------------

    def show_generation_status(self, fitting_run_name, generation_name, config, status_filepath=None):

        """
        This function ...
        :param fitting_run_name:
        :param generation_name:
        :param config:
        :param status_filepath:
        :return:
        """

        # Create simulation manager
        manager = GenerationManager()

        # Set fitting run and generation name
        manager.config.run = fitting_run_name
        manager.config.generation = generation_name

        # Set options
        manager.config.extra = config.extra
        manager.config.offline = config.offline
        manager.config.lazy = config.lazy
        manager.config.find_simulations = config.find_simulations
        manager.config.find_remotes = config.find_remotes
        manager.config.produce_missing = config.produce_missing
        manager.config.check_paths = config.check_paths
        manager.config.correct_paths = config.correct_paths
        manager.config.confirm_correction = config.confirm_correction
        manager.config.fix_success = config.fix_success
        manager.config.check_analysis = config.check_analysis
        manager.config.write_status = config.write_status
        manager.config.correct_status = config.correct_status

        # Not interactive
        manager.config.interactive = False

        # Set the status table filepath
        manager.config.status = status_filepath

        # Set the status command
        if config.extra is not None: status_command = "status " + ",".join(config.extra)
        else: status_command = "status"

        # Set status command
        manager.config.commands = [status_command]

        # Run the generation manager
        manager.run()

    # -----------------------------------------------------------------

    def show_generations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run names
        fitting_run_names = self.get_fitting_run_names_from_command(command, **kwargs)

        # Show
        self.show_generations(fitting_run_names)

    # -----------------------------------------------------------------

    def show_generations(self, fitting_run_names):

        """
        This function ...
        :param fitting_run_names:
        :return:
        """

        # Inform the user
        log.info("Showing generations ...")

        # Loop over the fitting runs
        for name in fitting_run_names:

            # Get the fitting run
            fitting_run = self.get_fitting_run(name)

            # Show
            print("")
            print(fmt.green + fmt.underlined + name + fmt.reset + ":")
            print("")

            # Loop over the generations
            for generation_name in fitting_run.generation_names: print(" - " + generation_name)

            print("")

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

        # Create the manager
        manager = AnalysisManager()
        manager.config.path = self.config.path

        # Set analysis run name
        manager.config.run = analysis_run_name

        # Run the manager
        manager.run()

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

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Flags
        definition.add_flag("backup", "make backups of non-empty directories", False)
        definition.add_flag("adapt_simulations", "unset retrieved and analysed flags of the corresponding simulations", False)

        # Generations to remove
        definition.add_required("generation", "string", "generation to remove (none means all)")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clear_generation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name and generation name
        fitting_run_name, generation_name, config = self.get_fitting_run_name_generation_name_and_config_from_command(command, self.clear_generation_definition, **kwargs)

        # Clear
        self.clear_generation(fitting_run_name, generation_name, backup=config.backup, adapt_simulations=config.adapt_simulations)

    # -----------------------------------------------------------------

    def clear_generation(self, fitting_run_name, generation_name, backup=False, adapt_simulations=False):

        """
        This function ...
        :param backup:
        :param adapt_simulations:
        :return:
        """

        # Load the fitting run
        fitting_run = self.get_fitting_run(fitting_run_name)

        # Check
        if not fitting_run.is_generation(generation_name): raise ValueError("Generation doesn't exist")

        # Get generation path
        generation_path = fitting_run.get_generation_path(generation_name)

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
                if backup: fs.backup_directory(out_path)
                fs.clear_directory(out_path)

            # Non-empty extracted data
            if fs.is_directory(extr_path) and not fs.is_empty(extr_path):

                # Debugging
                log.debug("Clearing extraction directory of '" + name + "' simulation ...")
                if backup: fs.backup_directory(extr_path)
                fs.clear_directory(extr_path)

            # Non-empty plotting output
            if fs.is_directory(plot_path) and not fs.is_empty(plot_path):

                # Debugging
                log.debug("Clearing plotting directory of '" + name + "' simulation ...")
                if backup: fs.backup_directory(plot_path)
                fs.clear_directory(plot_path)

            # Non-empty misc output
            if fs.is_directory(misc_path) and not fs.is_empty(misc_path):

                # Debugging
                log.debug("Clearing miscellaneous output directory of '" + name + "' simulation ...")
                if backup: fs.backup_directory(misc_path)
                fs.clear_directory(misc_path)

        # Clear chi-squared table
        chi_squared_path = fitting_run.chi_squared_table_path_for_generation(generation_name)
        chi_squared = fitting_run.chi_squared_table_for_generation(generation_name)

        # Non-empty
        if len(chi_squared) > 0:

            log.debug("Clearing chi-squared table ...")
            if backup: fs.backup_file(chi_squared_path)
            chi_squared.remove_all_rows()
            chi_squared.save()

        # Adapt simulations?
        if adapt_simulations:

            # Get the simulations for the generation
            generation = fitting_run.get_generation(generation_name)

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

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Generations to remove
        definition.add_optional("remote", "string", "remote host for which to remove the simulations (none means all")

        # Remove simulations of generations other than certain generations
        definition.add_optional("other_generations", "string_list", "remove all other generations than these generations")
        definition.add_flag("hard", "remove ALL simulations objects that are not from the 'other_generations' (use with CARE!)")

        # Make backup?
        definition.add_flag("backup", "make backup of simulation files")
        definition.add_optional("backup_path", "directory_path", "backup directory path")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clear_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the fitting run name, generation names and config
        fitting_run_name, generation_names, config = self.get_fitting_run_name_generation_names_and_config_from_command(command, self.clear_simulations_definition, **kwargs)

        # Clear
        self.clear_simulations(fitting_run_name, generation_names, hard=config.hard, remote=config.remote,
                               other_generations=config.other_generations, backup=config.backup, backup_path=config.backup_path)

    # -----------------------------------------------------------------

    def clear_simulations(self, fitting_run_name, generation_names, hard=False, remote=None,
                          other_generations=None, backup=False, backup_path=None):

        """
        This function ...
        :param fitting_run_name:
        :param generation_names:
        :param hard:
        :param remote:
        :param other_generations:
        :param backup:
        :param backup_path:
        :return:
        """

        # HARD delete
        if hard:

            if other_generations is None: raise ValueError("'other_generations' has to be specified")
            if generation_names is not None: raise ValueError("Generation names cannot be specified")

            # Clear
            self.clear_all_other_simulations(fitting_run_name, other_generations, remote=remote, backup=backup, backup_path=backup_path)

        # Look only at other generations of the same fitting run
        else: self.clear_simulations_generations(fitting_run_name, generation_names, other_generations=other_generations, remote=remote, backup=backup, backup_path=backup_path)

    # -----------------------------------------------------------------

    def clear_all_other_simulations(self, fitting_run_name, generation_names, remote=None, backup=False, backup_path=None):

        """
        This function ...
        :param fitting_run_name:
        :param generation_names:
        :param remote:
        :param backup:
        :param backup_path:
        :return:
        """

        # Load the fitting run
        fitting_run = self.get_fitting_run(fitting_run_name)

        # Initialize a dictionary with the simulation filepaths to keep per host
        keep_paths = defaultdict(list)

        # Loop over the generations of which we CANNOT remove the simulations
        for generation_name in generation_names:

            # Check
            if not fitting_run.is_generation(generation_name): raise ValueError("Generation '" + generation_name + "' doesn't exist")

            # Get generation path
            #generation_path = fitting_run.get_generation_path(generation_name)

            # Get the generation
            #generation = fitting_run.get_generation(generation_name)
            generation = self.get_generation(fitting_run_name, generation_name)

            # Loop over the remote hosts
            for host_id in generation.host_ids:

                # Get paths of the simulation files
                filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore", as_dict=False)
                nfilepaths = len(filepaths)

                # Add the filepaths
                keep_paths[host_id].extend(filepaths)

        # Loop over the remote hosts
        if remote is not None: host_ids = [remote]
        else: host_ids = keep_paths.keys()

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
            if backup:

                log.debug("Creating backup of the simulation files ...")
                if remote is not None: backup_path = backup_path
                else: backup_path = fs.create_directory_in(backup_path, host_id)
                fs.copy_files(remove_paths, backup_path)

            # Remove
            fs.remove_files(remove_paths)

    # -----------------------------------------------------------------

    def clear_simulations_generations(self, fitting_run_name, generation_names, other_generations=None, remote=None,
                                      backup=False, backup_path=None):

        """
        This function ...
        :param fitting_run_name:
        :param generation_names:
        :param other_generations:
        :param remote:
        :param backup:
        :param backup_path:
        :return:
        """

        # Load the fitting run
        fitting_run = self.get_fitting_run(fitting_run_name)

        # Determine generation names
        #if generation_names is not None:
            #if config.other_generations is not None: raise ValueError("Cannot specify also 'other_generations'")
            #generation_names = config.generations

        # Other generations
        #elif other_generations is not None:

        # Set generation names
        if generation_names is None:
            if other_generations is not None:

                # Get the other names than those specified
                generation_names = sequences.get_other(fitting_run.generation_names, other_generations)

            # Which generations?
            else: raise ValueError("Generation names cannot be determined")

        # Check
        elif other_generations is not None: raise ValueError("Cannot specify also 'other_generations'")

        # Show
        log.debug("Removing the simulation files for generations: '" + tostr(generation_names) + "' ...")

        # Loop over the generations
        for generation_name in generation_names:

            # Check
            if not fitting_run.is_generation(generation_name): raise ValueError("Generation '" + generation_name + "' doesn't exist")

            # Get generation path
            generation_path = fitting_run.get_generation_path(generation_name)

            # Get the generation
            generation = fitting_run.get_generation(generation_name)

            # Loop over the remote hosts
            for host_id in generation.host_ids:

                # Clear for this host?
                if remote is not None and host_id != remote: continue

                # Get paths of the simulation files
                filepaths = generation.get_simulation_paths_for_host(host_id, id_or_name="name", not_exist="ignore")
                nfilepaths = len(filepaths)
                if nfilepaths == 0:
                    log.warning("No simulation objects anymore for generation '" + generation_name + "'")
                    continue  # skip generation

                # Proceed?
                remove_names = [fs.strip_extension(fs.name(filepath)) for filepath in filepaths.keys()]
                if not prompt_proceed("proceed removing the simulation files of generation '" + generation_name + "': " + tostr(remove_names)): continue

                # Get the filepaths
                paths = filepaths.values()

                # Backup
                if backup:
                    log.debug("Creating backup of the simulation files of generation '" + generation_name + "' ...")
                    fs.copy_files(paths, backup_path)

                # Remove
                fs.remove_files(paths)

    # -----------------------------------------------------------------

    @lazyproperty
    def clone_fitting_definition(self):

        """
        This function ...
        :return:
        """

        # Create configuration definition
        definition = ConfigurationDefinition(write_config=False)

        # New name
        definition.add_required("name", "string", "name for the new fitting run", forbidden=self.fitting_run_names)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clone_fitting_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get fitting run name
        fitting_run_name, config = self.get_fitting_run_name_and_config_from_command(command, self.clone_fitting_definition, **kwargs)

        # Clone
        self.clone_fitting_run(fitting_run_name, config.name)

    # -----------------------------------------------------------------

    def clone_fitting_run(self, fitting_run_name, name):

        """
        This function ...
        :param fitting_run_name:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Cloning fitting run ...")

        # Get original fitting run
        fitting_run = self.get_fitting_run(fitting_run_name)

        # Clone the fitting run
        clone_fitting_run(fitting_run, name)

    # -----------------------------------------------------------------

    @lazyproperty
    def clone_analysis_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # New name
        definition.add_optional("name", "string", "name for the new analysis run", forbidden=self.analysis_run_names)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clone_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get analysis run name
        analysis_run_name, config = self.get_analysis_run_name_and_config_from_command(command, self.clone_analysis_definition, **kwargs)

        # Clone
        self.clone_analysis_run(analysis_run_name)

    # -----------------------------------------------------------------

    def clone_analysis_run(self, analysis_run_name, name=None):

        """
        This function ...
        :param analysis_run_name:
        :param name:
        :return:
        """

        # Inform the user
        log.info("Cloning analysis run ...")

        # Get analysis run path
        path = self.analysis_runs.get_path(analysis_run_name)

        # Generate new analysis run name
        if name is not None: new_name = name
        else: new_name = time.unique_name()
        new_path = fs.create_directory_in(self.environment.analysis_path, new_name)

        clear_directories = ["out", "misc", "plot", "residuals", "heating", "extr", "colours", "attenuation"]
        fs.copy_from_directory(path, new_path, exact_not_name=clear_directories, not_extension="sh")
        fs.create_directories_in(new_path, clear_directories)

        # Change the info
        info_path = fs.join(new_path, "info.dat")
        info = AnalysisRunInfo.from_file(info_path)

        # Set properties
        info.name = new_name
        info.path = new_path

        # Save
        info.save()

        # Change the input
        input_path = fs.join(new_path, "input.dat")
        input = load_dict(input_path)

        # Change the dust tree path
        tree_path = fs.join(new_path, "dust grid", "tree.dat")
        input["tree.dat"] = tree_path

        # Change the wavelength grid path
        wavelengths_path = fs.join(new_path, "wavelength_grid.dat")
        input["wavelengths.txt"] = wavelengths_path

        # Write again
        write_dict(input, input_path)

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

info_kwargs = dict()
info_kwargs["path"] = False
info_kwargs["name"] = False
info_kwargs["xsize"] = True
info_kwargs["ysize"] = True
info_kwargs["psf_filter"] = False
info_kwargs["filesize"] = True
info_kwargs["filter"] = False
info_kwargs["wavelength"] = False
info_kwargs["unit"] = False

# -----------------------------------------------------------------

def show_component_map_info(environment, name, filepath, origins=None, steps_path=None):

    """
    This function ...
    :param environment:
    :param name:
    :param filepath:
    :param origins:
    :param steps_path:
    :return:
    """

    from ...magic.tools.info import get_image_info_from_header_file
    from ...core.basics.composite import SimplePropertyComposite

    # Get map info
    info = get_image_info_from_header_file(name, filepath, **info_kwargs)

    # Make property composite
    info = SimplePropertyComposite.from_dict(info)

    # Show the properties
    print(info.to_string(line_prefix="   ", bullet="*"))

    # Show origins
    print("    * " + fmt.bold + "origins" + fmt.reset + ": " + tostr(origins, delimiter=", "))

    # Show origin with highest pixelscale
    highest_pixelscale_filter = get_highest_pixelscale_filter(environment, origins)
    print("    * " + fmt.bold + "pixelgrid reference" + fmt.reset + ": " + str(highest_pixelscale_filter))

    # Show origin with poorest resolution
    highest_fwhm_filter = get_highest_fwhm_filter(environment, origins)
    print("    * " + fmt.bold + "PSF reference" + fmt.reset + ": " + str(highest_fwhm_filter))

    # Load clip info
    clip_info_path = fs.join(steps_path, "clipped.dat")
    clip_info = load_dict(clip_info_path)

    # Show clip image references
    clip_origins = clip_info["origins"]
    print("    * " + fmt.bold + "clipping image references" + fmt.reset + ": " + tostr(clip_origins, delimiter=", "))

    # Show levels
    clipping_path = fs.join(steps_path, "clipping")
    levels = OrderedDict()
    for origin in clip_origins:
        mask_filename = "mask__" + tostr(origin).replace(" ", "_")
        mask_filename = fs.find_file_in_path(clipping_path, extension="fits", startswith=mask_filename, returns="name")
        level = float(mask_filename.split("__")[2])
        levels[origin] = level
    print("    * " + fmt.bold + "signal-to-noise levels" + fmt.reset + ":")
    for origin in levels: print("       " + str(origin) + ": " + str(levels[origin]))

    # Show
    softened = clip_info["soften"]
    print("    * " + fmt.bold + "softened clip mask" + fmt.reset + ": " + str(softened))

# -----------------------------------------------------------------

def get_highest_pixelscale_filter(environment, filters):

    """
    This function ...
    :param environment:
    :param filters:
    :return:
    """

    highest_pixelscale = None
    highest_pixelscale_filter = None

    for fltr in filters:

        image_name = str(fltr)

        # Show origin with highest pixelscale
        # path = environment.get_frame_path_for_filter(fltr)
        path = environment.get_frame_path(image_name)  # FASTER

        info = get_image_info_from_header_file(image_name, path, name=False, filter=False, wavelength=False,
                                               unit=False, fwhm=False, psf_filter=False, xsize=True, ysize=True,
                                               filesize=False, pixelscale=True)
        pixelscale = info["Pixelscale"]

        if highest_pixelscale is None or pixelscale > highest_pixelscale:
            highest_pixelscale = pixelscale
            highest_pixelscale_filter = fltr

    # Return
    return highest_pixelscale_filter

# -----------------------------------------------------------------

def get_highest_fwhm_filter(environment, filters):

    """
    Thisf unction ...
    :param environment:
    :param filters:
    :return:
    """

    highest_fwhm = None
    highest_fwhm_filter = None

    for fltr in filters:

        image_name = str(fltr)

        path = environment.get_frame_path(image_name)

        info = get_image_info_from_header_file(image_name, path, name=False, filter=False, wavelength=False,
                                               unit=False, fwhm=True, psf_filter=True, xsize=False, ysize=False,
                                               filesize=False, pixelscale=False)
        fwhm = info["FWHM"]

        if highest_fwhm is None or fwhm > highest_fwhm:
            highest_fwhm = fwhm
            highest_fwhm_filter = fltr

    # Return
    return highest_fwhm_filter

# -----------------------------------------------------------------
