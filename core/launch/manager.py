#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.manager Contains the SimulationManager class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import traceback
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..basics.configuration import ConfigurationDefinition, parse_arguments, get_usage, get_help, prompt_choice, prompt_settings
from ..basics.log import log
from ..tools.stringify import tostr
from .batchlauncher import SimulationAssignmentTable, SimulationStatusTable
from ..tools import formatting as fmt
from .timing import TimingTable
from .memory import MemoryTable
from ..tools import filesystem as fs
from ..tools import numbers
from ..basics.distribution import Distribution
from ..plot.distribution import plot_distribution
from ..tools.utils import lazyproperty, memoize_method, memoize_method_reset
from ..basics.configuration import prompt_proceed
from .analyser import show_analysis_steps, analyse_simulation, reanalyse_simulation, has_analysed
from ..simulation.remote import get_simulations_for_host
from ..tools import types
from ..basics.containers import DefaultOrderedDict
from ..tools import strings
from ..simulation.remote import SKIRTRemote
from ..remote.host import load_host
from ..basics.containers import create_nested_defaultdict, create_subdict
from ..tools import sequences
from ..basics.log import no_debugging
from ..simulation.remote import is_finished_status, is_running_status, finished_name, is_invalid_or_unknown_status, retrieved_name, analysed_name
from ..basics.configuration import prompt_string, prompt_variable
from ..remote.host import find_host_ids
from ..simulation.shower import show_simulation, show_analysis, compare_simulations, compare_analysis
from ..simulation.adapter import adapt_simulation, adapt_analysis
from ..config.show_simulation_settings import definition as show_simulation_definition
from ..config.adapt_simulation_settings import definition as adapt_simulations_definition
from ..config.show_analysis_options import definition as show_analysis_definition
from ..config.adapt_analysis_options import definition as adapt_analysis_definition
from .batchlauncher import BatchLauncher
from ..basics.table import SmartTable
from ..tools import tables
from ..tools import time
from ..config.analyse_simulation import definition as analyse_simulation_definition
from .analyser import all_steps
from .options import LoggingOptions, SchedulingOptions
from ..remote.load import show_status
from ..remote.mounter import mount_remote
from ..simulation.skifile import show_normalizations
from ..prep.summarize import show_instrument, show_stellar_component, show_dust_component
from ..units.unit import get_common_unit
from ..plot.sed import SEDPlotter
from ..data.sed import load_multiple_seds
from ..misc.fluxes import get_sed_instrument_name
from ..misc.images import get_datacube_instrument_name
from ..simulation.status import show_log_summary
from ..tools import introspection
from ..data.sed import ObservedSED
from ..plot.timeline import plot_timeline
from ..extract.timeline import TimeLineTable, extract_timeline
from ..extract.memory import MemoryUsageTable, extract_memory
from ..units.unit import parse_unit as u
from ..simulation.output import write_cache_path, output_type_choices, extraction_output_type_choices, plotting_output_type_choices, misc_output_type_choices
from ..simulation.output import output_types as output_type_names
from ..simulation.output import misc_output_types as misc_type_names
from ..simulation.jobscript import SKIRTJobScript
from ..simulation.screen import ScreenScript
from ..simulation.logfile import has_valid_timestamp
from ..remote.ensemble import SKIRTRemotesEnsemble

# -----------------------------------------------------------------

# Define types
datacube_types = [output_type_names.total_images, output_type_names.count_images, output_type_names.direct_images, output_type_names.transparent_images, output_type_names.scattered_images, output_type_names.dust_images, output_type_names.dust_scattered_images]
image_types = [misc_type_names.images, misc_type_names.images_for_fluxes]

# -----------------------------------------------------------------

class InvalidCommandError(Exception):

    """
    This class ...
    """

    def __init__(self, message, command):

        """
        Thisf unction ...
        :param message:
        :param command:
        """

        # Call the base class constructor with the parameters it needs
        super(InvalidCommandError, self).__init__(message)

        # The command
        self.command = command

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

plot_x_limits = (0., None)

# -----------------------------------------------------------------

clear_analysis_steps = ["extraction", "plotting", "misc"]

# -----------------------------------------------------------------

_help_command_name = "help"
_history_command_name = "history"
_status_command_name = "status"
_find_command_name = "find"
_hosts_command_name = "hosts"
_parallelizations_command_name = "parallelizations"
_info_command_name = "info"
_open_command_name = "open"
_sed_command_name = "sed"
_datacube_command_name = "datacube"
_fluxes_command_name = "fluxes"
_images_command_name = "images"
_input_command_name = "input"
_output_command_name = "output"
_extraction_command_name = "extraction"
_plotting_command_name = "plotting"
_misc_command_name = "misc"
_instruments_command_name = "instruments"
_stellar_command_name = "stellar"
_dust_command_name = "dust"
_normalizations_command_name = "normalizations"
_show_command_name = "show"
_plot_command_name = "plot"
_move_command_name = "move"
_stop_command_name = "stop"
_cancel_command_name = "cancel"
_remove_command_name = "remove"
_clear_command_name = "clear"
_cache_command_name = "cache"
_unfinish_command_name = "unfinish"
_unretrieve_command_name = "unretrieve"
_unanalyse_command_name = "unanalyse"
_relaunch_command_name = "relaunch"
_log_command_name = "log"
_error_command_name = "error"
_settings_command_name = "settings"
_analysis_command_name = "analysis"
_adapt_command_name = "adapt"
_compare_command_name = "compare"
_retrieve_command_name = "retrieve"
_analyse_command_name = "analyse"
_reanalyse_command_name = "reanalyse"
_assignment_command_name = "assignment"
_runtimes_command_name = "runtimes"
_memory_command_name = "memory"
_timeline_command_name = "timeline"
_base_command_name = "base"
_simulation_command_name = "simulation"

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history", False, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show simulations status", None)
commands[_find_command_name] = ("find_simulations_command", True, "find particular simulation(s) based on info or other parameters", None)
commands[_hosts_command_name] = ("show_hosts_command", True, "show remote hosts of the simulations", "hosts")
commands[_parallelizations_command_name] = ("show_parallelizations_command", True, "show parallelization schemes used per host", "host")
commands[_info_command_name] = ("show_info_command", True, "show info about the simulation (if defined in info tables)", "simulation")
commands[_open_command_name] = (None, None, "open input, output or base simulation directory", "simulation")
commands[_sed_command_name] = ("plot_seds_command", True, "plot SED(s) of a simulation", "simulation")
commands[_datacube_command_name] = ("plot_datacubes_command", True, "plot datacube(s) of a simulation", "simulation")
commands[_fluxes_command_name] = ("plot_fluxes_command", True, "plot mock fluxes calculated for a simulation", "simulation")
commands[_images_command_name] = ("plot_images_command", True, "plot mock images created for a simulation", "simulation")
commands[_input_command_name] = ("show_input_command", True, "show simulation input files", "simulation")
commands[_output_command_name] = ("show_output_command", True, "show simuation output files", "simulation")
commands[_extraction_command_name] = ("show_extraction_command", True, "show simulation extraction output files", "simulation")
commands[_plotting_command_name] = ("show_plotting_command", True, "show simulation plotting output files", "simulation")
commands[_misc_command_name] = ("show_misc_command", True, "show simulation misc output files and directories", "simulation")
commands[_instruments_command_name] = ("show_instruments_command", True, "show simulation instruments", "simulation")
commands[_stellar_command_name] = ("show_stellar_components_command", True, "show stellar components", "simulation")
commands[_dust_command_name] = ("show_dust_components_command", True, "show dust components", "simulation")
commands[_normalizations_command_name] = ("show_normalizations_command", True, "show model normalizations", "simulation")
commands[_show_command_name] = (None, None, "show", None)
commands[_plot_command_name] = (None, None, "plot", None)
commands[_move_command_name] = ("move_simulations_command", True, "move simulations from one remote to another", "simulations")
commands[_stop_command_name] = ("stop_simulations_command", True, "stop running simulations", "simulation")
commands[_cancel_command_name] = ("cancel_simulations_command", True, "cancel queued simulations", "simulations")
commands[_remove_command_name] = ("remove_simulation_command", True, "remove simulation", "simulations")
commands[_clear_command_name] = (None, None, "clear simulation output/input/analysis", "simulation")
commands[_cache_command_name] = (None, None, "cache simulation/analysis output to another directory/filesystem", "simulations")
commands[_unfinish_command_name] = ("unfinish_simulation_command", True, "remove remote output of a simulation and unset finished flag", "simulation")
commands[_unretrieve_command_name] = ("unretrieve_simulation_command", True, "remove local (retrieved) output of a simulation and unset retrieved flag", "simulation")
commands[_unanalyse_command_name] = ("unanalyse_simulation_command", True, "remove analysis output of a simulation and unset analysis flags", "simulation")
commands[_relaunch_command_name] = ("relaunch_simulation_command", True, "relaunch simulations on the original remote host", "simulation")
commands[_log_command_name] = ("show_simulation_log_command", True, "show log output of a simulation", "simulation")
commands[_error_command_name] = ("show_simulation_errors_command", True, "show error output of a simulation", "simulation")
commands[_settings_command_name] = ("show_simulation_settings_command", True, "show simulation settings", "simulation")
commands[_analysis_command_name] = ("show_analysis_options_command", True, "show analysis options", "simulation")
commands[_adapt_command_name] = (None, None, "adapt simulation settings or analysis options", "simulation")
commands[_compare_command_name] = (None, None, "compare simulation settings or analysis options between two simulations", "two_simulations")
commands[_retrieve_command_name] = ("retrieve_simulation_command", True, "retrieve a simulation from the remote host", "simulation")
commands[_analyse_command_name] = ("analyse_simulation_command", True, "analyse a simulation", "simulation")
commands[_reanalyse_command_name] = ("reanalyse_simulation_command", True, "re-analyse a simulation", "simulation")

# -----------------------------------------------------------------

# Define show commands
show_commands = OrderedDict()
show_commands[_assignment_command_name] = ("show_assignment", False, "show the simulation assignment scheme", None)
show_commands[_status_command_name] = ("show_status", False, "show the simulation status", None)
show_commands[_runtimes_command_name] = ("show_runtimes_command", True, "show the simulation runtimes", "host_parallelization")
show_commands[_memory_command_name] = ("show_memory_command", True, "show the simulation memory usages", "host_parallelization")

# -----------------------------------------------------------------

# Define plot commands
plot_commands = OrderedDict()
plot_commands[_runtimes_command_name] = ("plot_runtimes_command", True, "plot simulation runtimes", "host_parallelization")
plot_commands[_memory_command_name] = ("plot_memory_command", True, "plot simulation memory usages", "host_parallelization")
plot_commands[_timeline_command_name] = ("plot_timeline_command", True, "plot simulation timeline", "simulation")

# -----------------------------------------------------------------

# Define open commands
open_commands = OrderedDict()
open_commands[_base_command_name] = ("open_base_command", True, "open simulation base directory", "simulation")
open_commands[_input_command_name] = ("open_input_command", True, "open simulation input directory", "simulation")
open_commands[_output_command_name] = ("open_output_command", True, "open simulation output directory", "simulation")
open_commands[_extraction_command_name] = ("open_extraction_command", True, "open simulation extraction output directory", "simulation")
open_commands[_plotting_command_name] = ("open_plotting_command", True, "open simulation plotting output directory", "simulation")
open_commands[_misc_command_name] = ("open_misc_command", True, "open simulation miscellaneous output directory", "simulation")

# -----------------------------------------------------------------

# Define adapt commands
adapt_commands = OrderedDict()
adapt_commands[_simulation_command_name] = ("adapt_simulation_settings_command", True, "adapt simulation settings", "simulation")
adapt_commands[_analysis_command_name] = ("adapt_analysis_options_command", True, "adapt analysis options", "simulation")

# -----------------------------------------------------------------

# Define compare commands
compare_commands = OrderedDict()
compare_commands[_simulation_command_name] = ("compare_simulation_settings_command", True, "compare simulation settings", "two_simulations")
compare_commands[_analysis_command_name] = ("compare_analysis_options_command", True, "compare analysis options", "two_simulations")

# -----------------------------------------------------------------

# Define clear commands
clear_commands = OrderedDict()
clear_commands[_input_command_name] = ("clear_simulation_input_command", True, "clear simulation input", "simulation")
clear_commands[_output_command_name] = ("clear_simulation_output_command", True, "clear simulation output", "simulation")
clear_commands[_analysis_command_name] = ("clear_simulation_analysis_command", True, "clear simulation analysis output", "simulation")
clear_commands[_extraction_command_name] = ("clear_extraction_command", True, "clear simulation extraction output", "simulation")
clear_commands[_plotting_command_name] = ("clear_plotting_command", True, "clear simulation plotting output", "simulation")
clear_commands[_misc_command_name] = ("clear_misc_command", True, "clear simulation misc output", "simulation")

# -----------------------------------------------------------------

# Define cache commands
cache_commands = OrderedDict()
cache_commands[_output_command_name] = ("cache_simulations_output_command", True, "cache simulation output", "simulations")
cache_commands[_analysis_command_name] = ("cache_simulations_analysis_command", True, "cache simulation analysis output", "simulations")
cache_commands[_extraction_command_name] = ("cache_simulations_extraction_command", True, "cache simulation extraction output", "simulations")
cache_commands[_plotting_command_name] = ("cache_simulations_plotting_command", True, "cache simulation plotting output", "simulations")
cache_commands[_misc_command_name] = ("cache_simulations_misc_command", True, "cache simulation misc output", "simulations")
cache_commands[_datacube_command_name] = ("cache_simulations_datacubes_command", True, "cache datacubes", "simulations")
cache_commands[_images_command_name] = ("cache_simulations_images_command", True, "cache images", "simulations")

# -----------------------------------------------------------------

_simulation_subject_name = "simulation"
_simulations_subject_name = "simulations"
_two_simulations_subject_name = "two_simulations"
_host_subject_name = "host"
_hosts_subject_name = "hosts"
_host_parallelization_subject_name = "host_parallelization"

# -----------------------------------------------------------------

class MovedSimulationsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Original host ID"] = (str, None, "original host ID")
    _column_info["Original clustername"] = (str, None, "original cluster name")
    _column_info["Original ID"] = (int, None, "original simulation ID")
    _column_info["New host ID"] = (str, None, "new host ID")
    _column_info["New clustername"] = (str, None, "new cluster name")
    _column_info["New ID"] = (int, None, "new simulation ID")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MovedSimulationsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Find index of the simulation
        return tables.find_index(self, simulation_name)

    # -----------------------------------------------------------------

    def add_simulation(self, simulation, new_host, new_id=None):

        """
        This function ...
        :param simulation:
        :param new_host:
        :param new_id:
        :return:
        """

        # Set values
        values = [simulation.name, simulation.host_id, simulation.cluster_name, simulation.id, new_host.id, new_host.cluster_name, new_id]

        # Add the row
        self.add_row(values)

    # -----------------------------------------------------------------

    def set_new_id(self, simulation_name, new_id):

        """
        This function ...
        :param simulation_name:
        :param new_id:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        self.set_value("New ID", index, new_id)

# -----------------------------------------------------------------

class RelaunchedSimulationsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Original ID"] = (int, None, "original simulation ID")
    _column_info["New ID"] = (int, None, "new simulation ID")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RelaunchedSimulationsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Find index of the simulation
        return tables.find_index(self, simulation_name)

    # -----------------------------------------------------------------

    def add_simulation(self, simulation, new_id=None):

        """
        This function ...
        :param simulation:
        :param new_id:
        :return:
        """

        # Set values
        values = [simulation.name, simulation.id, new_id]

        # Add the row
        self.add_row(values)

    # -----------------------------------------------------------------

    def set_new_id(self, simulation_name, new_id):

        """
        This function ...
        :param simulation_name:
        :param new_id:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        self.set_value("New ID", index, new_id)

# -----------------------------------------------------------------

_screen_extra_name = "screen"
_job_extra_name = "job"
_disk_extra_name = "disk"
_runtime_extra_name = "runtime"
_memory_extra_name = "memory"

# -----------------------------------------------------------------

# Define extra columns
extra_columns = OrderedDict()
extra_columns[_screen_extra_name] = "screen session name"
extra_columns[_job_extra_name] = "job ID"
extra_columns[_disk_extra_name] = "simulation output disk size"
extra_columns[_runtime_extra_name] = "total simulation runtime"
extra_columns[_memory_extra_name] = "peak simulation memory usage"

# Define extra column names
extra_column_names = dict()
extra_column_names[_screen_extra_name] = "Screen name"
extra_column_names[_job_extra_name] = "Job ID"
extra_column_names[_disk_extra_name] = "Disk size"
extra_column_names[_runtime_extra_name] = "Runtime"
extra_column_names[_memory_extra_name] = "Peak memory"

# Define extra column units
extra_column_units = dict()
extra_column_units[_disk_extra_name] = "GB"
extra_column_units[_runtime_extra_name] = "h"
extra_column_units[_memory_extra_name] = "GB"

# -----------------------------------------------------------------

class SimulationManager(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SimulationManager, self).__init__(*args, **kwargs)

        # The simulation assignment scheme
        self.assignment = None

        # Flags: has the assignment scheme been adapted or created?
        self._adapted_assignment = False
        self._new_assignment = False

        # The simulations
        self.simulations = DefaultOrderedDict(OrderedDict)

        # Timing and memory table
        self.timing = None
        self.memory = None

        # The moved simulations
        self.moved = MovedSimulationsTable()

        # The relaunched simulations
        self.relaunched = RelaunchedSimulationsTable()

        # The commands that have been executed
        self.commands = []

        # The batch simulation launcher
        self.launcher = BatchLauncher()

        # Additional info tables
        self.info = OrderedDict()

        # The ensemble of remotes
        self.remotes = SKIRTRemotesEnsemble()

        # Mount paths
        self.mount_paths = dict()

    # -----------------------------------------------------------------

    @property
    def ninfo(self):

        """
        This function ...
        :return:
        """

        return len(self.info)

    # -----------------------------------------------------------------

    @property
    def has_info(self):

        """
        This function ...
        :return:
        """

        return self.ninfo > 0

    # -----------------------------------------------------------------

    @memoize_method
    def get_info_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Initialize dictionary for the info
        info = OrderedDict()

        # Loop over the tables
        for name in self.info:

            # Get the table
            table = self.info[name]

            # Check
            if "Simulation name" not in table.column_names:
                log.warning("Table '" + name + "' does not contain a column with the simulation names: skipping ...")
                continue

            # Find the index of the simulation
            index = table.find_index(simulation_name, column_name="Simulation name")
            if index is None: continue  # skip this table

            # Loop over the column names
            for column_name in table.column_names:

                # Check column name
                if column_name == "Simulation name": continue
                if column_name in info:
                    log.warning("Already obtained the '" + column_name + "' parameter from another table: skipping ...")
                    continue

                # Get the value
                value = table.get_value(column_name, index)

                # Set the value
                info[column_name] = value

        # Return the info for the simulation
        return info

    # -----------------------------------------------------------------

    def ninfo_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return len(self.get_info_for_simulation(simulation_name))

    # -----------------------------------------------------------------

    def has_info_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.ninfo_for_simulation(simulation_name) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def info_names(self):

        """
        This function ...
        :return:
        """

        names = DefaultOrderedDict(int)
        for simulation_name in self.simulation_names:
            if not self.has_info_for_simulation(simulation_name): continue
            info = self.get_info_for_simulation(simulation_name)
            for name in info: names[name] += 1
        return names.keys()

    # -----------------------------------------------------------------

    @memoize_method
    def get_info_values(self, name):

        """
        This function returns the info value for a particular info property name, for each simulation
        :param name:
        :return:
        """

        values = dict()

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Get the value (or None if not defined)
            if self.has_info_for_simulation(simulation_name):
                info = self.get_info_for_simulation(simulation_name)
                if name in info: value = info[name]
                else: value = None
            else: value = None

            # Set the value
            values[simulation_name] = value

        # Return the values
        return values

    # -----------------------------------------------------------------

    @memoize_method
    def get_info_values_scalar(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the unit
        unit = self.info_units[name]

        # Get the values
        values = self.get_info_values(name)

        # Initialize dictionary for the scalar values
        scalar = dict()

        # Loop over the simulation
        for simulation_name in values:
            value = values[simulation_name]
            if value is not None and unit is not None: value = value.to(unit).value
            scalar[simulation_name] = value

        # Return the scalar values
        return scalar

    # -----------------------------------------------------------------

    @lazyproperty
    def info_units(self):

        """
        This function ...
        :return:
        """

        units = dict()

        # Loop over the info properties
        for name in self.info_names:

            # Get the values for the different simulations
            values = self.get_info_values(name)

            # Get the common unit
            unit = get_common_unit(values.values())

            # Set the unit
            units[name] = unit

        # Return the units
        return units

    # -----------------------------------------------------------------

    @property
    def info_unit_strings(self):

        """
        This function ...
        :return:
        """

        strings = []
        for name in self.info_names:
            unit = self.info_units[name]
            if unit is None: string = ""
            else: string = tostr(unit)
            strings.append(string)
        return strings

    # -----------------------------------------------------------------

    @property
    def has_moving(self):

        """
        This function ...
        :return:
        """

        return self.config.move

    # -----------------------------------------------------------------

    @property
    def has_showing(self):

        """
        This function ...
        :return:
        """

        return self.config.show and (self.config.show_assignment or self.config.show_status or self.config.show_runtimes or self.config.show_memory)

    # -----------------------------------------------------------------

    @property
    def has_plotting(self):

        """
        This function ...
        :return:
        """

        return self.config.plot and (self.config.plot_runtimes or self.config.plot_memory)

    # -----------------------------------------------------------------

    @property
    def has_writing(self):

        """
        This function ...
        :return:
        """

        return self.config.write and (self.config.write_assignment or self.config.write_status or self.config.write_moved)

    # -----------------------------------------------------------------

    @property
    def has_analysis(self):

        """
        This function ...
        :return:
        """

        return self.config.analyse

    # -----------------------------------------------------------------

    @property
    def has_reanalysis(self):

        """
        This function ...
        :return:
        """

        return self.config.reanalysis is not None

    # -----------------------------------------------------------------

    @property
    def has_any(self):

        """
        This function ...
        :return:
        """

        return self.has_moving or self.has_showing or self.has_plotting or self.has_writing or self.has_analysis or self.has_reanalysis

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

        if self.config.interactive is None: return not self.has_any and not self.do_commands
        else: return self.config.interactive

    # -----------------------------------------------------------------

    @property
    def do_moving(self):

        """
        This function ...
        :return:
        """

        return self.config.move

    # -----------------------------------------------------------------

    @property
    def do_launch(self):

        """
        This function ...
        :return:
        """

        return self.launcher.has_queued

    # -----------------------------------------------------------------

    @property
    def do_showing(self):

        """
        This function ...
        :return:
        """

        return self.config.show

    # -----------------------------------------------------------------

    @property
    def do_writing(self):

        """
        This function ...
        :return:
        """

        return self.config.write

    # -----------------------------------------------------------------

    @property
    def do_plotting(self):

        """
        This function ...
        :return:
        """

        return self.config.plot

    # -----------------------------------------------------------------

    @property
    def do_reanalysis(self):

        """
        This function ...
        :return:
        """

        return self.config.reanalyse is not None

    # -----------------------------------------------------------------

    @property
    def do_analysis(self):

        """
        This function ...
        :return:
        """

        return self.config.analyse

    # -----------------------------------------------------------------

    @property
    def do_caching(self):

        """
        This function ...
        :return:
        """

        return self.config.cache_output or self.config.cache_datacubes or self.config.cache_misc or self.config.cache_images

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

        # 4. Move simulations
        if self.do_moving: self.move()

        # 5. Launch
        if self.do_launch: self.launch()

        # 6. Show
        if self.do_showing: self.show()

        # 7. Write
        if self.do_writing: self.write()

        # 8. Plotting
        if self.do_plotting: self.plot()

        # 9. Re-analyse
        if self.do_reanalysis: self.reanalyse()

        # 10. Analyse
        if self.do_analysis: self.analyse()

        # 11. Cache
        if self.do_caching: self.cache()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationManager, self).setup(**kwargs)

        # Check options
        if self.config.prompt_simulations_reanalysis is None: self.config.prompt_simulations_reanalysis = self.config.reanalyse_simulations is None
        if self.config.prompt_simulations_move is None: self.config.prompt_simulations_move = self.config.move_simulations is None

        # Check caching options
        if self.config.cache_output and self.config.cache_datacubes: self.config.cache_datacubes = False # datacubes will be cached, as all output
        if self.config.cache_misc and self.config.cache_images: self.config.cache_images = False # images will be cached, as all misc output

        # Get the status
        if "status" in kwargs: self.status = kwargs.pop("status")
        elif self.config.status is not None: self.status = SimulationStatusTable.from_file(self.config.status)

        # Initialize simulations and assignment scheme
        self.initialize(**kwargs)

        # Get timing table
        self.get_timing_table(**kwargs)

        # Get memory table
        self.get_memory_table(**kwargs)

        # Get the info tables
        if "info_tables" in kwargs:
            tables = kwargs.pop("info_tables")
            if types.is_dictionary(tables): self.info.update(tables)
            elif types.is_sequence(tables):
                for table in tables: self.info[table.filename] = table
            else: raise ValueError("Invalid type for 'info_tables'")
        elif self.config.info_tables is not None:
            for path in self.config.info_tables:
                filename = fs.strip_extension(fs.name(path))
                self.info[filename] = SmartTable.from_file(path)

        # Get remotes
        if "remotes" in kwargs:
            remotes = kwargs.pop("remotes")
            if types.is_dictionary(remotes):
                for name in remotes: self.set_remote(remotes[name], name)
            elif types.is_sequence(remotes):
                for remote in remotes: self.set_remote(remote)
            else: raise ValueError("Invalid type for 'remotes'")

        # Set launcher options
        self.set_launcher_options()

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting options for the batch simulation launcher ...")

        # Set the working directory for the batch launcher
        self.launcher.config.path = self.config.path

        # Set options
        self.launcher.config.dry = self.config.dry
        self.launcher.config.retrieve = False
        self.launcher.config.analyse = False
        self.launcher.config.shared_input = self.config.shared_input

        # Showing
        self.launcher.config.show = True
        self.launcher.config.show_info = True
        self.launcher.config.show_status = False
        self.launcher.config.show_parallelizations = True

        # Writing
        self.launcher.config.write = True
        self.launcher.config.write_assignment = False
        self.launcher.config.write_queues = True

        # Clear existing (for moving and relaunching)
        self.launcher.config.clear_existing = True

    # -----------------------------------------------------------------

    def set_remote_input_path_for_host(self, host_id, path):

        """
        This function ...
        :param host_id:
        :param path:
        :return:
        """

        self.launcher.set_remote_input_path_for_host(host_id, path)

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_path(self):

        """
        This function ...
        :return:
        """

        if self.config.backup_path is not None:
            if self.config.backup_dir_path is not None: raise ValueError("Cannot specify both 'backup_path' and 'backup_dir_path'")
            return self.config.backup_path
        elif self.config.backup_dir_path is not None:
            if self.config.backup_dirname is not None: return fs.create_directory_in(self.config.backup_dir_path, self.config.backup_dirname)
            else: return fs.create_directory_in(self.config.backup_dir_path, time.unique_name("manage"))
        else:
            if self.config.backup_dirname is not None: return self.output_path_directory(self.config.backup_dirname, create=True)
            else: return self.output_path_directory(time.unique_name("manage"), create=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_simulations_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.backup_path, "simulations")

    # -----------------------------------------------------------------

    @property
    def backup_assignment_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.backup_path, "assignment.dat")

    # -----------------------------------------------------------------

    @memoize_method
    def get_backup_simulations_path_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return fs.create_directory_in(self.backup_simulations_path, host_id)

    # -----------------------------------------------------------------

    @property
    def has_timing(self):

        """
        This function ...
        :return:
        """

        return self.timing is not None

    # -----------------------------------------------------------------

    @property
    def has_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory is not None

    # -----------------------------------------------------------------

    @property
    def all_simulations(self):

        """
        This function ...
        :return:
        """

        for host_id in self.simulations:
            for simulation_name in self.simulations[host_id]:
                yield self.simulations[host_id][simulation_name]

    # -----------------------------------------------------------------

    @property
    def all_retrieved_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if not simulation.retrieved: continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_not_retrieved_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if simulation.retrieved: continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_analysed_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if not simulation.analysed: continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_not_analysed_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if simulation.analysed: continue
            yield simulation

    # -----------------------------------------------------------------

    def get_timing_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "timing" in kwargs: self.timing = kwargs.pop("timing")
        elif self.config.timing is not None: self.timing = TimingTable.from_file(self.config.timing)
        else:

            table_paths = []
            for simulation in self.all_simulations: table_paths.append(simulation.analysis.timing_table_path)
            timing_table_path = sequences.get_all_equal_value(table_paths, ignore_none=True, return_none=True)
            if timing_table_path is not None: self.timing = TimingTable.from_file(timing_table_path)

    # -----------------------------------------------------------------

    def get_memory_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "memory" in kwargs: self.memory = kwargs.pop("memory")
        elif self.config.memory is not None: self.memory = MemoryTable.from_file(self.config.memory)
        else:

            table_paths = []
            for simulation in self.all_simulations: table_paths.append(simulation.analysis.memory_table_path)
            memory_table_path = sequences.get_all_equal_value(table_paths, ignore_none=True, return_none=True)
            if memory_table_path is not None: self.memory = MemoryTable.from_file(memory_table_path)

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Initializing simulations ...")

        # Load the simulation assignment scheme
        if "assignment" in kwargs: self.initialize_from_assignment(kwargs.pop("assignment"))

        # From simulations
        elif "simulations" in kwargs: self.initialize_from_simulations(kwargs.pop("simulations"))

        # From simulation names
        elif "simulation_names" in kwargs: self.initialize_from_simulation_names(kwargs.pop("simulation_names"))

        # From assignment
        elif self.config.assignment is not None:

            assignment = SimulationAssignmentTable.from_file(self.config.assignment)
            self.initialize_from_assignment(assignment)

        # From simulation names
        elif self.config.simulation_names is not None: self.initialize_from_simulation_names(self.config.simulation_names)

        # From directories
        elif self.config.from_directories:
            names = fs.directories_in_path(self.config.path, returns="name")
            self.initialize_from_simulation_names(names)

        # Not enough input
        else: raise ValueError("Not enough input to initialize simulations")

    # -----------------------------------------------------------------

    @lazyproperty
    def simulations_for_hosts(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        simulations = dict()

        # Get simulations for each remote host
        for host_id in self.config.remotes:

            # Get the simulations
            simulations_host = get_simulations_for_host(host_id, as_dict=True)

            # Set the simulations
            simulations[host_id] = simulations_host

        # Return
        return simulations

    # -----------------------------------------------------------------

    def find_host_id_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        the_host_id = None

        # Loop over the remotes
        for host_id in self.simulations_for_hosts:

            # Check whether the simulation name is in the
            if simulation_name in self.simulations_for_hosts[host_id]:
                the_host_id = host_id
                break

        # Return the host ID
        return the_host_id

    # -----------------------------------------------------------------

    @memoize_method
    def host_id_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.find_host_id_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def get_simulation(self, simulation_name, host_id=None):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :return:
        """

        # Get the host ID if necessary
        if host_id is None: host_id = self.host_id_for_simulation(simulation_name)

        # Checks
        if host_id not in self.simulations: raise ValueError("No simulations for host '" + host_id + "'")
        if simulation_name not in self.simulations[host_id]: raise ValueError("Simulation '" + simulation_name + "' not in simulations for host '" + host_id + "'")

        # Get the simulation
        return self.simulations[host_id][simulation_name]

    # -----------------------------------------------------------------

    def get_simulation_prefix(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).prefix()

    # -----------------------------------------------------------------

    def get_simulation_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).path

    # -----------------------------------------------------------------

    def get_log_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).log_file_path

    # -----------------------------------------------------------------

    def has_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).has_logfile

    # -----------------------------------------------------------------

    @memoize_method
    def get_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).log_file

    # -----------------------------------------------------------------

    @memoize_method
    def get_logfiles(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).logfiles()

    # -----------------------------------------------------------------

    def get_runtime(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_logfile(simulation_name): return self.get_logfile(simulation_name).total_runtime * u("s")
        else: return None

    # -----------------------------------------------------------------

    def get_peak_memory(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_logfile(simulation_name): return self.get_logfile(simulation_name).peak_memory * u("GB")
        else: return None

    # -----------------------------------------------------------------

    def get_disk_size(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        output = self.get_output_disk_size(simulation_name)
        extraction = self.get_extraction_disk_size(simulation_name)
        plotting = self.get_plotting_disk_size(simulation_name)
        misc = self.get_misc_disk_size(simulation_name)
        return output + extraction + plotting + misc

    # -----------------------------------------------------------------

    def get_output_disk_size(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_output(simulation_name).disk_size

    # -----------------------------------------------------------------

    def get_extraction_disk_size(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_extraction_output(simulation_name).disk_size

    # -----------------------------------------------------------------

    def get_plotting_disk_size(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_plotting_output(simulation_name).disk_size

    # -----------------------------------------------------------------

    def get_misc_disk_size(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).disk_size

    # -----------------------------------------------------------------

    def get_simulation_host_id(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).host_id

    # -----------------------------------------------------------------

    def get_execution_handle(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).handle

    # -----------------------------------------------------------------

    def is_screen_execution(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_execution_handle(simulation_name).is_screen

    # -----------------------------------------------------------------

    def is_job_execution(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_execution_handle(simulation_name).is_job

    # -----------------------------------------------------------------

    def get_screen_name(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Not screen
        if not self.is_screen_execution(simulation_name): return None

        # Return the screen name
        return self.get_execution_handle(simulation_name).value

    # -----------------------------------------------------------------

    def get_remote_screen_output_path(self, simulation_name):

        """
        This function ...
        :return:
        """

        # Not screen
        if not self.is_screen_execution(simulation_name): return None

        # Return the output path
        return self.get_execution_handle(simulation_name).remote_screen_output_path

    # -----------------------------------------------------------------

    def get_remote_screen_script_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Not screen
        if not self.is_screen_execution(simulation_name): return None

        # Return the output path
        return self.get_execution_handle(simulation_name).remote_screen_script_path

    # -----------------------------------------------------------------

    def has_remote_screen_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_remote_screen_script_path(simulation_name) is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def local_screen_scripts(self):

        """
        This property ...
        :return:
        """

        # Initialize dictionary
        scripts = dict()
        if self.config.screen_scripts is None: return scripts

        # Loop over the paths
        for filepath in self.config.screen_scripts:

            # Load script
            script = ScreenScript.from_file(filepath)

            # Add under the name
            scripts[script.name] = script

        # Return
        return scripts

    # -----------------------------------------------------------------

    @memoize_method
    def get_local_screen_script_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the screen name
        screen_name = self.get_screen_name(simulation_name)

        # Check whether there is a local screen script with this name
        if screen_name not in self.local_screen_scripts: return None

        # Extra check: has simulation?
        script = self.local_screen_scripts[screen_name]
        if not script.has_simulation(simulation_name): raise RuntimeError("The screen script at '" + script.path + "' does not contain the simulation '" + simulation_name + "'")

        # Return the screen script path
        return script.path

    # -----------------------------------------------------------------

    def has_local_screen_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_local_screen_script_path(simulation_name) is not None

    # -----------------------------------------------------------------

    @memoize_method
    def get_screen_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Local screen script found
        if self.has_local_screen_script(simulation_name):

            # Return the screen script
            screen_name = self.get_screen_name(simulation_name)
            return self.local_screen_scripts[screen_name]

        # Remote screen script found
        elif self.has_remote_screen_script(simulation_name):

            # Open the remote screen script
            host_id = self.get_host_id_for_simulation(simulation_name)
            screen_script_path = self.get_remote_screen_script_path(simulation_name)
            return ScreenScript.from_remote_file(screen_script_path, remote=self.get_remote(host_id))

        # No screen script found
        else: return None

    # -----------------------------------------------------------------

    def has_screen_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_screen_script(simulation_name) is not None

    # -----------------------------------------------------------------

    def get_screen_log_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Try to get remote output path
        output_path = self.get_remote_screen_output_path(simulation_name)
        if output_path is None and self.has_screen_script(simulation_name): output_path = self.get_screen_script(simulation_name).output_path
        if output_path is None: return None

        # Get the remote
        remote = self.get_remote_for_simulation(simulation_name)

        # Try to find screenlog file
        filename = "screenlog.0"
        filepath = fs.join(output_path, filename)
        if not remote.is_file(filepath): return None

        # Return the screen log path
        return filepath

    # -----------------------------------------------------------------

    def get_job_id(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Not job
        if not self.is_job_execution(simulation_name): return None

        # Return the job ID
        return self.get_execution_handle(simulation_name).value

    # -----------------------------------------------------------------

    def get_job_id_string(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        job_id = self.get_job_id(simulation_name)
        if job_id is None: return None
        else: return str(job_id)

    # -----------------------------------------------------------------

    def get_remote_job_script_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Not job
        if not self.is_job_execution(simulation_name): return None

        # Return the script path
        return self.get_execution_handle(simulation_name).remote_job_script_path

    # -----------------------------------------------------------------

    def has_remote_job_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_remote_job_script_path(simulation_name) is not None

    # -----------------------------------------------------------------

    @memoize_method
    def get_local_job_script_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Look for sh files in simulation base directory
        base_path = self.get_base_path(simulation_name)

        # Look for sh files in the simulation path
        filepaths = fs.files_in_path(base_path, extension="sh")
        nfiles = len(filepaths)

        # Return the single filepath
        if nfiles == 0: return None
        elif nfiles > 1:
            log.warning("Multiple potential job script files found in the simulation directory of simulation '" + simulation_name + "'")
            return None
        else: return filepaths[0]

    # -----------------------------------------------------------------

    def has_local_job_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_local_job_script_path(simulation_name) is not None

    # -----------------------------------------------------------------

    @memoize_method
    def get_job_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Local job script is found
        if self.has_local_job_script(simulation_name):

            filepath = self.get_local_job_script_path(simulation_name)
            return SKIRTJobScript.from_file(filepath)

        # Remote job script is found
        elif self.has_remote_job_script(simulation_name):

            host_id = self.get_host_id_for_simulation(simulation_name)
            filepath = self.get_remote_job_script_path(simulation_name)
            return SKIRTJobScript.from_remote_file(filepath, remote=self.get_remote(host_id))

        # No job script is found
        else: return None

    # -----------------------------------------------------------------

    def has_job_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_job_script(simulation_name) is not None

    # -----------------------------------------------------------------

    def get_job_output_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_job_script(simulation_name): return self.get_job_script(simulation_name).output_path
        else:
            remote_output_path = self.get_simulation(simulation_name).remote_output_path
            filepaths = self.get_remote_for_simulation(simulation_name).files_in_path(remote_output_path, extension="txt", contains="out")
            nfiles = len(filepaths)
            if nfiles == 0: return None
            elif nfiles > 1: return None
            else: return filepaths[0]

    # -----------------------------------------------------------------

    @memoize_method
    def get_job_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        filepath = self.get_job_output_path(simulation_name)
        if filepath is None: return []
        return self.get_remote_for_simulation(simulation_name).get_lines(filepath)

    # -----------------------------------------------------------------

    def get_job_error_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_job_script(simulation_name): return self.get_job_script(simulation_name).error_path
        else:
            remote_output_path = self.get_simulation(simulation_name).remote_output_path
            filepaths = self.get_remote_for_simulation(simulation_name).files_in_path(remote_output_path, extension="txt", contains="err")
            nfiles = len(filepaths)
            if nfiles == 0: return None
            elif nfiles > 1: return None
            else: return filepaths[0]

    # -----------------------------------------------------------------

    @memoize_method
    def get_job_error(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        filepath = self.get_job_error_path(simulation_name)
        if filepath is None: return []
        return self.get_remote_for_simulation(simulation_name).get_lines(filepath)

    # -----------------------------------------------------------------

    def get_simulation_extraction_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).extraction_path

    # -----------------------------------------------------------------

    def get_simulation_plotting_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).plotting_path

    # -----------------------------------------------------------------

    def get_simulation_misc_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).misc_path

    # -----------------------------------------------------------------

    @memoize_method
    def get_input(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).input

    # -----------------------------------------------------------------

    @memoize_method
    def get_remote_input(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        simulation = self.get_simulation(simulation_name)
        return simulation.get_remote_input(remote=self.get_remote(simulation.host_id))

    # -----------------------------------------------------------------

    def get_base_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).base_path

    # -----------------------------------------------------------------

    def get_remote_base_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).remote_simulation_path

    # -----------------------------------------------------------------

    def get_output_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).output_path

    # -----------------------------------------------------------------

    def get_remote_output_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).remote_output_path

    # -----------------------------------------------------------------

    def get_extraction_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).extraction_path

    # -----------------------------------------------------------------

    def get_plotting_path(self, simulation_name):

        """
        THis function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).plotting_path

    # -----------------------------------------------------------------

    def get_misc_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).misc_path

    # -----------------------------------------------------------------

    @memoize_method
    def get_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).output

    # -----------------------------------------------------------------

    @memoize_method
    def get_noutput(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return len(self.get_output(simulation_name))

    # -----------------------------------------------------------------

    def has_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_noutput(simulation_name) > 0

    # -----------------------------------------------------------------

    def get_noutput_files(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_output(simulation_name).get_nfiles_types(otypes)

    # -----------------------------------------------------------------

    def get_noutput_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_output(simulation_name).get_ndirectories_types(otypes)

    # -----------------------------------------------------------------

    def get_noutput_files_and_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_noutput_files(simulation_name, otypes) + self.get_noutput_directories(simulation_name, otypes)

    # -----------------------------------------------------------------

    def has_datacubes(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_noutput_files(simulation_name, datacube_types) > 0

    # -----------------------------------------------------------------

    def is_cached_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_output(simulation_name).has_cached

    # -----------------------------------------------------------------

    @memoize_method
    def get_extraction_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).extraction_output

    # -----------------------------------------------------------------

    @memoize_method
    def get_nextraction_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return len(self.get_extraction_output(simulation_name))

    # -----------------------------------------------------------------

    def has_extraction_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_nextraction_output(simulation_name) > 0

    # -----------------------------------------------------------------

    def get_nextraction_files(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_extraction_output(simulation_name).get_nfiles_types(otypes)

    # -----------------------------------------------------------------

    def get_nextraction_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_extraction_output(simulation_name).get_ndirectories_types(otypes)

    # -----------------------------------------------------------------

    def get_nextraction_files_and_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_nextraction_files(simulation_name, otypes) + self.get_nextraction_directories(simulation_name, otypes)

    # -----------------------------------------------------------------

    def is_cached_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_extraction_output(simulation_name).has_cached

    # -----------------------------------------------------------------

    def has_analysis_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.has_extraction_output(simulation_name) or self.has_plotting_output(simulation_name) or self.has_misc_output(simulation_name)

    # -----------------------------------------------------------------

    def has_timeline_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        extraction = self.get_extraction_output(simulation_name)
        return extraction.has_single_timeline

    # -----------------------------------------------------------------

    @memoize_method
    def get_timeline_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Timeline has already been created
        if self.has_timeline_for_simulation(simulation_name):

            extraction = self.get_extraction_output(simulation_name)
            return TimeLineTable.from_file(extraction.single_timeline)

        # Create timeline
        else: return extract_timeline(self.get_simulation(simulation_name))

    # -----------------------------------------------------------------

    def has_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        extraction = self.get_extraction_output(simulation_name)
        return extraction.has_single_memory

    # -----------------------------------------------------------------

    @memoize_method
    def get_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Memory table has already been created
        if self.has_memory_for_simulation(simulation_name):

            extraction = self.get_extraction_output(simulation_name)
            return MemoryUsageTable.from_file(extraction.single_memory)

        # Create memory table
        else: return extract_memory(self.get_simulation(simulation_name))

    # -----------------------------------------------------------------

    @memoize_method
    def get_plotting_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).plotting_output

    # -----------------------------------------------------------------

    @memoize_method
    def get_nplotting_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return len(self.get_plotting_output(simulation_name))

    # -----------------------------------------------------------------

    def has_plotting_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_nplotting_output(simulation_name) > 0

    # -----------------------------------------------------------------

    def get_nplotting_files(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_plotting_output(simulation_name).get_nfiles_types(otypes)

    # -----------------------------------------------------------------

    def get_nplotting_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_plotting_output(simulation_name).get_ndirectories_types(otypes)

    # -----------------------------------------------------------------

    def get_nplotting_files_and_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_nplotting_files(simulation_name, otypes) + self.get_nplotting_directories(simulation_name, otypes)

    # -----------------------------------------------------------------

    def is_cached_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_plotting_output(simulation_name).has_cached

    # -----------------------------------------------------------------

    def has_sed_plots_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        plotting = self.get_plotting_output(simulation_name)
        return plotting.has_seds

    # -----------------------------------------------------------------

    def get_sed_plot_paths_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        plotting = self.get_plotting_output(simulation_name)
        return plotting.seds

    # -----------------------------------------------------------------

    @memoize_method
    def get_misc_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).misc_output

    # -----------------------------------------------------------------

    @memoize_method
    def get_nmisc_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return len(self.get_misc_output(simulation_name))

    # -----------------------------------------------------------------

    def has_misc_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_nmisc_output(simulation_name) > 0

    # -----------------------------------------------------------------

    def get_nmisc_files(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_misc_output(simulation_name).get_nfiles_types(otypes)

    # -----------------------------------------------------------------

    def get_nmisc_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_misc_output(simulation_name).get_ndirectories_types(otypes)

    # -----------------------------------------------------------------

    def get_nmisc_files_and_directories(self, simulation_name, otypes):

        """
        This function ...
        :param simulation_name:
        :param otypes:
        :return:
        """

        return self.get_nmisc_files(simulation_name, otypes) + self.get_nmisc_directories(simulation_name, otypes)

    # -----------------------------------------------------------------

    def has_images(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_nmisc_files_and_directories(simulation_name, image_types) > 0

    # -----------------------------------------------------------------

    def is_cached_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_misc_output(simulation_name).has_cached

    # -----------------------------------------------------------------

    def has_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_fluxes

    # -----------------------------------------------------------------

    @memoize_method
    def get_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return ObservedSED.from_file(misc.single_fluxes)

    # -----------------------------------------------------------------

    def has_image_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_image_fluxes

    # -----------------------------------------------------------------

    @memoize_method
    def get_image_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return ObservedSED.from_file(misc.single_image_fluxes)

    # -----------------------------------------------------------------

    def has_images_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_images

    # -----------------------------------------------------------------

    def get_images_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        from ...magic.core.frame import Frame

        # Initialize a dictionary for the images
        images = OrderedDict()

        # Get misc output
        misc = self.get_misc_output(simulation_name)

        # Loop over the files
        for filepath in misc.images:
            name = fs.strip_extension(fs.name(filepath))
            frame = Frame.from_file(filepath)
            images[name] = frame

        # Return the dictionary of images
        return images

    # -----------------------------------------------------------------

    def has_images_for_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_images_for_fluxes

    # -----------------------------------------------------------------

    def get_images_for_fluxes_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        from ...magic.core.frame import Frame

        # Initialize a dictionary for the images
        images = OrderedDict()

        # Get misc output
        misc = self.get_misc_output(simulation_name)

        # Loop over the files
        for filepath in misc.images_for_fluxes:
            name = fs.strip_extension(fs.name(filepath))
            frame = Frame.from_file(filepath)
            images[name] = frame

        # Return the dictionary of images
        return images

    # -----------------------------------------------------------------

    def has_fluxes_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_fluxes_plot

    # -----------------------------------------------------------------

    def get_fluxes_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.single_fluxes_plot

    # -----------------------------------------------------------------

    def has_image_fluxes_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_images_fluxes_plot

    # -----------------------------------------------------------------

    def get_image_fluxes_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.single_images_fluxes_plot

    # -----------------------------------------------------------------

    def has_images_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_images_plot

    # -----------------------------------------------------------------

    def get_images_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.single_images_plot

    # -----------------------------------------------------------------

    def has_images_for_fluxes_plot_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.has_single_images_for_fluxes_plot

    # -----------------------------------------------------------------

    def get_images_for_fluxes_plot_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        misc = self.get_misc_output(simulation_name)
        return misc.single_images_for_fluxes_plot

    # -----------------------------------------------------------------

    @memoize_method
    def get_remote_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        simulation = self.get_simulation(simulation_name)
        return simulation.get_remote_output(remote=self.get_remote(simulation.host_id))

    # -----------------------------------------------------------------

    @memoize_method
    def get_wavelength_grid(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        ski = self.get_skifile(simulation_name)
        input_paths = self.get_input(simulation_name)
        return ski.get_wavelengths(input_paths, as_grid=True)

    # -----------------------------------------------------------------

    @memoize_method
    def get_skifile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).ski_file

    # -----------------------------------------------------------------

    def get_skifile_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).ski_path

    # -----------------------------------------------------------------

    def get_remote_skifile_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).remote_ski_path

    # -----------------------------------------------------------------

    @memoize_method
    def get_instrument_names(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get skifile
        ski = self.get_skifile(simulation_name)

        # Return the names
        return ski.get_instrument_names()

    # -----------------------------------------------------------------

    @memoize_method
    def get_instruments(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Initialize instruments dictionary
        instruments = OrderedDict()

        # Get skifile
        ski = self.get_skifile(simulation_name)

        # Loop over the instrument names
        for name in self.get_instrument_names(simulation_name):

            # Get the instrument
            instrument = ski.get_instrument_object(name)

            # Add the instrument
            instruments[name] = instrument

        # Return the instruments dictionary
        return instruments

    # -----------------------------------------------------------------

    @memoize_method
    def get_stellar_component_ids(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get ski file
        ski = self.get_skifile(simulation_name)

        # Return names
        return ski.get_stellar_component_ids()

    # -----------------------------------------------------------------

    @memoize_method
    def get_dust_component_ids(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get ski file
        ski = self.get_skifile(simulation_name)

        # Return names
        return ski.get_dust_component_ids()

    # -----------------------------------------------------------------

    def get_host_id_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).host_id

    # -----------------------------------------------------------------

    def get_cluster_name_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).cluster_name

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_names_for_host_id(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        names = []

        # Loop over the simulation names
        for simulation_name in self.simulation_names:
            if self.get_host_id_for_simulation(simulation_name) != host_id: continue
            names.append(simulation_name)

        # Return the names
        return names

    # -----------------------------------------------------------------

    def get_remote_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_remote(self.get_host_id_for_simulation(simulation_name))

    # -----------------------------------------------------------------

    @memoize_method
    def get_host_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        simulation = self.get_simulation(simulation_name)
        return load_host(simulation.host_id, simulation.cluster_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_parallelization_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the simulation
        simulation = self.get_simulation(simulation_name) #, host_id=host_id)
        return simulation.parallelization

    # -----------------------------------------------------------------

    @memoize_method
    def get_parallelizations_for_host(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        parallelizations = set()

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Get the host for the simulation
            host_simulation = self.get_host_for_simulation(simulation_name)
            if host_simulation != host: continue

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)
            parallelizations.add(parallelization)

        # Return
        return list(parallelizations)

    # -----------------------------------------------------------------

    def has_simulation_for_host_id(self, host_id, simulation_name):

        """
        This function ...
        :param host_id:
        :param simulation_name:
        :return:
        """

        if host_id not in self.simulations_for_hosts: raise ValueError("Simulations for host '" + host_id + "' are not loaded")
        return simulation_name in self.simulations_for_hosts[host_id]

    # -----------------------------------------------------------------

    def initialize_from_assignment(self, assignment):

        """
        This function ...
        :param assignment:
        :return:
        """

        # Set the assignment scheme
        self.assignment = assignment

        # Loop over the simulations
        for index in range(self.assignment.nsimulations):

            # Get simulation name, host ID and cluster name
            simulation_name = self.assignment.get_simulation_name_for_index(index)
            host_id = self.assignment.get_host_id_for_index(index)
            #print(simulation_name, host_id)

            # Does this simulation exist?
            if self.has_simulation_for_host_id(host_id, simulation_name):

                # Get the simulation
                simulation = self.simulations_for_hosts[host_id][simulation_name]

                # Add the simulation
                self.simulations[host_id][simulation_name] = simulation

            # Simulation is not found for the host specified in the assignment scheme
            else:

                # Find actual host ID
                actual_host_id = self.find_host_id_for_simulation(simulation_name)

                # No host ID found
                if actual_host_id is None:

                    # Assume local?
                    if self.config.local:

                        # Give warning
                        if self.config.warn_local: log.warning("Simulation '" + simulation_name + "' has not been found for any remote host: assuming the simulation has run locally ...")

                        # Set local
                        self.assignment.set_local_for_simulation(simulation_name)
                        self._adapted_assignment = True

                    # Cannot assume local
                    else: raise ValueError("Cannot determine the host for simulation '" + simulation_name + "'")

                # Host ID found
                else:

                    # Give warning
                    log.warning("Host ID for simulation '" + simulation_name + "' in assignment scheme (" + host_id + ") is not correct: simulation found assigned to host '" + actual_host_id + "'")

                    # Get the simulation for the actual host ID
                    simulation = self.simulations_for_hosts[actual_host_id][simulation_name]

                    # Get simulation ID
                    actual_id = simulation.id

                    # Get cluster name
                    cluster_name = simulation.cluster_name

                    # Change the host and simulation ID for this simulation
                    self.assignment.set_id_and_host_for_simulation(simulation_name, actual_id, actual_host_id, cluster_name=cluster_name)
                    self._adapted_assignment = True

                    # Add the simulation
                    self.simulations[actual_host_id][simulation_name] = simulation

    # -----------------------------------------------------------------

    def initialize_from_simulations(self, simulations):

        """
        This function ...
        :param simulations:
        :return:
        """

        # Initialize assignment table
        self.assignment = SimulationAssignmentTable()
        self._new_assignment = True

        # Get list of simulations
        if types.is_dictionary(simulations): simulations = simulations.values()
        elif types.is_sequence(simulations): pass
        else: raise ValueError("Invalid type for 'simulations'")

        # Loop over the simulations
        for simulation in simulations:

            # Get simulation properties
            simulation_name = simulation.name
            simulation_id = simulation.id
            host_id = simulation.host_id
            cluster_name = simulation.cluster_name

            # Add to assignment
            self.assignment.add_remote_simulation(simulation_name, host_id, cluster_name=cluster_name,
                                                  simulation_id=simulation_id, success=self.config.success)

            # Add the simulation
            self.simulations[host_id][simulation_name] = simulation

    # -----------------------------------------------------------------

    def initialize_from_simulation_names(self, simulation_names):

        """
        This function ...
        :param simulation_names:
        :return:
        """

        # Initialize assignment table
        self.assignment = SimulationAssignmentTable()
        self._new_assignment = True

        # Loop over the simulation names and look for matches
        for simulation_name in simulation_names:

            # Find the remote host ID
            host_id = self.find_host_id_for_simulation(simulation_name)

            # No remote host for this simulation
            if host_id is None:

                # Assume local?
                if self.config.local:

                    # Give warning
                    if self.config.warn_local: log.warning("Simulation '" + simulation_name + "' has not been found for any remote host: assuming the simulation has run locally ...")

                    # Set local
                    self.assignment.add_local_simulation(simulation_name, success=self.config.success)

                # Cannot assume local
                else: raise ValueError("Cannot determine the host for simulation '" + simulation_name + "'")

            # Remote host was found
            else:

                # Get the simulation
                simulation = self.simulations_for_hosts[host_id][simulation_name]

                # Get simulation properties
                simulation_id = simulation.id
                cluster_name = simulation.cluster_name

                # Add to assignment
                self.assignment.add_remote_simulation(simulation_name, host_id, cluster_name=cluster_name, simulation_id=simulation_id, success=self.config.success)

                # Add the simulation
                self.simulations[host_id][simulation_name] = simulation

    # -----------------------------------------------------------------

    @lazyproperty
    def host_ids(self):

        """
        This function ...
        :return:
        """

        return self.assignment.unique_host_ids

    # -----------------------------------------------------------------

    @lazyproperty
    def hosts(self):

        """
        This function ...
        :return:
        """

        return self.assignment.unique_hosts

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.assignment.names

    # -----------------------------------------------------------------

    @lazyproperty
    def status(self):

        """
        This function ...
        :return:
        """

        # Initialize lists
        status_list = []

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Get the simulation
            simulation = self.get_simulation(simulation_name)

            # Analysed
            if simulation.analysed: simulation_status = "analysed"

            # Partly analysed
            elif simulation.analysed_any:

                analysed = []
                if simulation.analysed_all_extraction: analysed.append("extraction")
                if simulation.analysed_all_plotting: analysed.append("plotting")
                if simulation.analysed_all_misc: analysed.append("misc")
                if simulation.analysed_batch: analysed.append("batch")
                if simulation.analysed_scaling: analysed.append("scaling")
                if simulation.analysed_all_extra: analysed.append("extra")

                if len(analysed) > 0: simulation_status = "analysed: " + ", ".join(analysed)
                else: simulation_status = "analysed: started"

            # Retrieved
            elif simulation.retrieved: simulation_status = "retrieved"

            # Finished
            elif simulation.finished: simulation_status = "finished"

            # Not yet retrieved
            else:

                host_id = simulation.host_id
                screen_states = self.screens[host_id]
                jobs_status = self.jobs[host_id]
                with no_debugging(): simulation_status = self.get_remote(host_id).get_simulation_status(simulation, screen_states=screen_states, jobs_status=jobs_status)

            # Check success flag in assignment
            if self.config.fix_success and not self.assignment.is_launched(simulation.name) and not is_invalid_or_unknown_status(simulation_status):
                log.warning("Settting the launch of simulation '" + simulation.name + "' as succesful in the assignment table as this was not yet done")
                self.assignment.set_success_for_simulation(simulation.name)

            # Retrieve finished simulations?
            if simulation_status == finished_name and self.config.retrieve:
                self.retrieve_simulation(simulation_name)
                simulation_status = "retrieved"

            # Add the status
            status_list.append(simulation_status)

        # Create the table and return
        return SimulationStatusTable.from_columns(self.simulation_names, status_list)

    # -----------------------------------------------------------------

    def retrieve_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Retrieve
        self.retrieve_simulation(simulation_name)

    # -----------------------------------------------------------------

    def retrieve_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Retrieving simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Get the remote and retrieve the simulation
        remote = self.get_remote(simulation.host_id)
        remote.retrieve_simulation(simulation)

    # -----------------------------------------------------------------

    def reset_status(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Resetting the simulation status ...")

        # Reset the screen states and job status
        self.remotes.reset_jobs()
        self.remotes.reset_screens()

        # Reset the simulation status table
        del self.status

    # -----------------------------------------------------------------

    @property
    def has_status(self):

        """
        This function ...
        :return:
        """

        return "status" in self.__dict__

    # -----------------------------------------------------------------

    @property
    def screens(self):

        """
        This function ...
        :return:
        """

        return self.remotes.screens

    # -----------------------------------------------------------------

    @property
    def jobs(self):

        """
        This function ...
        :return:
        """

        return self.remotes.jobs

    # -----------------------------------------------------------------

    def get_status(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.get_status(simulation_name)

    # -----------------------------------------------------------------

    def is_queued(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_queued(simulation_name)

    # -----------------------------------------------------------------

    def is_running(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_running(simulation_name)

    # -----------------------------------------------------------------

    def is_queued_or_running(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_queued_or_running(simulation_name)

    # -----------------------------------------------------------------

    def is_finished(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_finished(simulation_name)

    # -----------------------------------------------------------------

    def is_running_or_finished(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_running_or_finished(simulation_name)

    # -----------------------------------------------------------------

    def is_retrieved(self, simulation_name):

        """
        This function ....
        :param simulation_name:
        :return:
        """

        return self.status.is_retrieved(simulation_name)

    # -----------------------------------------------------------------

    def is_analysed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_analysed(simulation_name)

    # -----------------------------------------------------------------

    def is_cancelled(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_cancelled(simulation_name)

    # -----------------------------------------------------------------

    def is_aborted(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_aborted(simulation_name)

    # -----------------------------------------------------------------

    def is_crashed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_crashed(simulation_name)

    # -----------------------------------------------------------------

    def is_running_or_finished_or_aborted_or_crashed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_running_or_finished_or_aborted_or_crashed(simulation_name)

    # -----------------------------------------------------------------

    def is_failed(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.status.is_failed(simulation_name)

    # -----------------------------------------------------------------

    @property
    def nfinished(self):

        """
        This function ...
        :return:
        """

        return self.status.nfinished

    # -----------------------------------------------------------------

    @property
    def relative_nfinished(self):

        """
        This function ...
        :return:
        """

        return self.status.relative_nfinished

    # -----------------------------------------------------------------

    @property
    def percentage_nfinished(self):

        """
        This function ...
        :return:
        """

        return self.status.percentage_nfinished

    # -----------------------------------------------------------------

    @property
    def nretrieved(self):

        """
        This function ...
        :return:
        """

        return self.status.nretrieved

    # -----------------------------------------------------------------

    @property
    def relative_nretrieved(self):

        """
        This function ...
        :return:
        """

        return self.status.relative_nretrieved

    # -----------------------------------------------------------------

    @property
    def percentage_nretrieved(self):

        """
        This function ...
        :return:
        """

        return self.status.percentage_nretrieved

    # -----------------------------------------------------------------

    @property
    def nanalysed(self):

        """
        This function ...
        :return:
        """

        return self.status.nanalysed

    # -----------------------------------------------------------------

    @property
    def relative_nanalysed(self):

        """
        This function ...
        :return:
        """

        return self.status.relative_nanalysed

    # -----------------------------------------------------------------

    @property
    def percentage_nanalysed(self):

        """
        This function ...
        :return:
        """

        return self.status.percentage_nanalysed

    # -----------------------------------------------------------------

    def run_commands(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running commands ...")

        # Loop over the commands
        for command in self.config.commands:

            # Debugging
            log.debug("Running '" + command + "' ...")

            # Process command, give error if fails
            self.process_command(command)

            # Add command
            self.commands.append(command)

    # -----------------------------------------------------------------

    def interactive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Entering interactive mode ...")

        # Enter loop
        while True:

            # Get next command, break if no command is given
            command = prompt_string("command", "command to be executed")
            if not command: break

            # DEVELOPER COMMAND
            if command.startswith("$"):
                self._run_developer_command(command)
                continue

            # Process command
            success = True
            try: self.process_command(command)
            except InvalidCommandError as e:
                log.warning("Invalid command: '" + e.command + "'")
                success = False
            except Exception as e:
                traceback.print_exc()
                log.error(str(e))
                success = False

            # Add command, if succesful
            if success: self.commands.append(command)

    # -----------------------------------------------------------------

    @property
    def ncommands(self):

        """
        This function ...
        :return:
        """

        return len(self.commands)

    # -----------------------------------------------------------------

    @property
    def has_commands(self):

        """
        This function ...
        :return:
        """

        return self.ncommands > 0

    # -----------------------------------------------------------------

    def has_subcommands(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        return command in [_show_command_name, _plot_command_name, _open_command_name, _adapt_command_name, _compare_command_name, _clear_command_name, _cache_command_name]

    # -----------------------------------------------------------------

    def get_subcommands(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        if command == _show_command_name: return show_commands
        elif command == _plot_command_name: return plot_commands
        elif command == _open_command_name: return open_commands
        elif command == _adapt_command_name: return adapt_commands
        elif command == _compare_command_name: return compare_commands
        elif command == _clear_command_name: return clear_commands
        elif command == _cache_command_name: return cache_commands
        else: raise InvalidCommandError("Invalid command: '" + command + "'", command)

    # -----------------------------------------------------------------

    def _run_developer_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Clean
        command = command[1:].strip()

        # Create evaluation command
        eval_command = "self." + command

        # Evaluate
        eval(eval_command)

    # -----------------------------------------------------------------

    def _run_command(self, command, cmds):

        """
        This function ...
        :param command:
        :param cmds:
        :return:
        """

        # Get first word
        first = command.split(" ")[0]

        # Check whether interactive
        if first.startswith("*"):
            interactive = True
            command = command[1:]
            first = first[1:]
        else: interactive = False

        # Find key
        if first not in cmds: raise InvalidCommandError("Invalid command: '" + first + "'", command)
        key = first

        # Has subcommands
        if self.has_subcommands(key): self._run_subcommand(command)

        # Regular command
        else: self._run_command_impl(command, cmds, interactive=interactive)

    # -----------------------------------------------------------------

    def _run_subcommand(self, command, interactive=False):

        """
        This function ...
        :param command:
        :param interactive:
        :return:
        """

        # Get first word == key
        key = command.split(" ")[0]

        # Get the possible subcommands
        subcommands = self.get_subcommands(key)

        # Get command without main command
        subcommand = strings.split_at_first(command, key)[1].strip()

        # No subcommand?
        if not strings.startswith_any(subcommand, subcommands):

            # Show help for the main command
            if "-h" in subcommand or "--help" in subcommand:

                # Set subcommands with descriptions
                subcommands_descriptions = OrderedDict()
                for subkey in subcommands:
                    function_name, pass_command, description, subject = subcommands[subkey]
                    subcommands_descriptions[subkey] = description

                # Create definition for the subcommands
                definition = ConfigurationDefinition(write_config=False)
                definition.add_required("subcommand", "string", "subcommand", choices=subcommands_descriptions)
                help = get_help(key, definition, add_logging=False, add_cwd=False)

                # Show help
                for line in help: print(fmt.red + line + fmt.reset)

            # Nothing more than the main command is given as input
            elif subcommand == "":

                # Interactive mode: prompt between the different subcommands
                if interactive:

                    # Set subcommands with descriptions
                    subcommands_descriptions = OrderedDict()
                    for subkey in subcommands:
                        function_name, pass_command, description, subject = subcommands[subkey]
                        subcommands_descriptions[subkey] = description

                    # Prompt for the subcommand
                    subcommand = prompt_string("subcommand", "subcommand", choices=subcommands_descriptions)

                    # Run the subcommand in interactive mode
                    self._run_command_impl(subcommand, subcommands, main_command=key, interactive=True)

                # Not enough input
                else: raise InvalidCommandError("Not enough input for '" + key + "' command", command)

            # Invalid command
            else: raise InvalidCommandError("Invalid command: '" + subcommand + "'", command)

        # Run the command
        else: self._run_command_impl(subcommand, subcommands, main_command=key, interactive=interactive)

    # -----------------------------------------------------------------

    def _run_command_impl(self, command, cmds, main_command=None, interactive=False):

        """
        This function ...
        :param command:
        :param cmds:
        :param main_command:
        :param interactive:
        :return:
        """

        # Get first word == key
        key = command.split(" ")[0]

        # Get function name and description
        function_name, pass_command, description, subject = cmds[key]

        # Show help for the command
        if "-h" in command or "--help" in command:

            # Get the help info
            help = self.get_help_for_key(key, cmds, main_command=main_command)

            # Show help
            if help is None: print(fmt.red + "no input required" + fmt.reset)
            else:
                for line in help: print(fmt.red + line + fmt.reset)

        # Actually run
        else:

            # Get the function
            function = getattr(self, function_name)

            # Call the function
            if pass_command: function(command, interactive=interactive)
            else: function(interactive=interactive)

    # -----------------------------------------------------------------

    def process_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Run the command
        self._run_command(command, commands)

    # -----------------------------------------------------------------

    def get_definition_for_function_name(self, function_name):

        """
        This function ...
        :param function_name:
        :return:
        """

        # Get definition
        definition_property_name = function_name.split("_command")[0] + "_definition"
        definition = getattr(self, definition_property_name, None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def get_kwargs_for_function_name(self, function_name):

        """
        This function ...
        :param function_name:
        :return:
        """

        # Get kwargs
        kwargs_property_name = function_name.split("_command")[0] + "_kwargs"
        kwargs = getattr(self, kwargs_property_name, {})

        # Return
        return kwargs

    # -----------------------------------------------------------------

    def get_usage_for_key(self, key, cmds, main_command=None):

        """
        This function ...
        :param key:
        :param cmds:
        :param main_command:
        :return:
        """

        # Get properties
        function_name, pass_command, description, subject = cmds[key]
        #if subject is None: return None # no usage for this command (simple command that doesn't need input)

        # Get the definition
        definition = self.get_definition_for_function_name(function_name)

        # Get the kwargs
        kwargs = self.get_kwargs_for_function_name(function_name)

        # Set command name
        if main_command is None: name = key
        else: name = main_command + " " + key

        # Get usage lines
        if subject is None: usage = self.get_usage_command(definition, name=name)
        elif subject == _simulation_subject_name: usage = self.get_usage_simulation_command(definition, name=name, **kwargs)
        elif subject == _simulations_subject_name: usage = self.get_usage_simulations_command(definition, name=name, **kwargs)
        elif subject == _two_simulations_subject_name: usage = self.get_usage_two_simulations_command(definition, name=name, **kwargs)
        elif subject == _host_subject_name: usage = self.get_usage_host_command(definition, name=name, **kwargs)
        elif subject == _hosts_subject_name: usage = self.get_usage_hosts_command(definition, name=name, **kwargs)
        elif subject == _host_parallelization_subject_name: usage = self.get_usage_host_and_parallelization_command(definition, name=name, **kwargs)
        else: raise ValueError("Invalid subject '" + subject + "'")

        # Return the usage info
        return usage

    # -----------------------------------------------------------------

    def get_help_for_key(self, key, cmds, main_command=None):

        """
        This function ...
        :param key:
        :param cmds:
        :param main_command:
        :return:
        """

        # Get properties
        function_name, pass_command, description, subject = cmds[key]
        if subject is None: return None # no help for this command (simple command that doesn't need input)

        # Get the definition
        definition = self.get_definition_for_function_name(function_name)

        # Get the kwargs
        kwargs = self.get_kwargs_for_function_name(function_name)

        # Set command name
        if main_command is None: name = key
        else: name = main_command + " " + key

        # Get help lines
        if subject is None: help = self.get_help_command(definition, name=name)
        elif subject == _simulation_subject_name: help = self.get_help_simulation_command(definition, name=name, **kwargs)
        elif subject == _simulations_subject_name: help = self.get_help_simulations_command(definition, name=name, **kwargs)
        elif subject == _two_simulations_subject_name: help = self.get_help_two_simulations_command(definition, name=name, **kwargs)
        elif subject == _host_subject_name: help = self.get_help_host_command(definition, name=name, **kwargs)
        elif subject == _hosts_subject_name: help = self.get_help_hosts_command(definition, name=name, **kwargs)
        elif subject == _host_parallelization_subject_name: help = self.get_help_host_and_parallelization_command(definition, name=name, **kwargs)
        else: raise ValueError("Invalid subject '" + subject + "'")

        # Return the help info
        return help

    # -----------------------------------------------------------------

    def show_help(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing help ...")

        print("")

        # Loop over the commands
        for key in commands:

            # Get properties
            function_name, pass_command, description, subject = commands[key]

            # Show
            print(" - " + fmt.bold + key + fmt.reset + ": " + description)

            if function_name is None and not self.has_subcommands(key): raise RuntimeError("Something went wrong")

            # Show usage
            if pass_command:

                # Get usage
                usage = self.get_usage_for_key(key, commands)

                # Show
                for line in usage: print("    " + fmt.blue + line + fmt.reset)

            # Command with subcommands
            elif self.has_subcommands(key):

                # Set subcommands with descriptions
                subcommands = self.get_subcommands(key)
                subcommands_descriptions = OrderedDict()
                for subkey in subcommands:
                    function_name, pass_command, description, subject = subcommands[subkey]
                    subcommands_descriptions[subkey] = description

                # Create definition for the subcommands
                definition = ConfigurationDefinition(write_config=False)
                definition.add_required("subcommand", "string", "subcommand", choices=subcommands_descriptions)
                usage = get_usage(key, definition, add_logging=False, add_cwd=False)

                # Show usage
                for line in usage: print("    " + fmt.blue + line + fmt.reset)

                # Show help for subcommands
                self._show_help_subcommands(subcommands, key)

            # No input needed
            else: print("    " + fmt.blue + "no input" + fmt.reset)

        print("")

    # -----------------------------------------------------------------

    def _show_help_subcommands(self, subcommands, main_command):

        """
        This function ...
        :param subcommands:
        :param main_command:
        :return:
        """

        # Loop over the commands
        for key in subcommands:

            # Get description
            function_name, pass_command, description, subject = subcommands[key]

            # Show
            print("    * " + fmt.bold + key + fmt.reset + ": " + description)

            # Show usage
            # if not pass_command or subject is None: continue
            if pass_command:

                # Get usage
                usage = self.get_usage_for_key(key, subcommands, main_command=main_command)

                # Show
                for line in usage: print("      " + fmt.blue + line + fmt.reset)

            # No input expected
            else: print("        " + fmt.blue + "no input" + fmt.reset)

    # -----------------------------------------------------------------

    def show_history(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing history of commands ...")

        print("")
        for command in self.commands:
            print(command)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_hosts_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("properties", "show host properties")
        definition.add_flag("status", "show host status")

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def show_hosts_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required"] = False
        kwargs["choices"] = self.hosts
        return kwargs

    # -----------------------------------------------------------------

    def show_hosts_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing hosts ...")

        # Get all kwargs
        kwargs.update(self.show_hosts_kwargs)

        # Get the host
        splitted, hosts, config = self.parse_hosts_command(command, command_definition=self.show_hosts_definition, **kwargs)

        print("")

        # Loop over the hosts
        for host in hosts: self.show_host(host, properties=config.properties, status=config.status)

        print("")

    # -----------------------------------------------------------------

    def show_host(self, host, properties=False, status=False):

        """
        This function ...
        :param host:
        :param properties:
        :param status:
        :return:
        """

        # Show name
        print(fmt.blue + fmt.underlined + tostr(host) + fmt.reset)

        # Show properties?
        if properties:

            # Show whether scheduler
            print(" - " + fmt.bold + "uses scheduling system: " + fmt.reset + str(host.scheduler))

            # Show cluster
            if host.has_clusters:

                info = host.cluster
                numa_domains_per_node = info.numa_domains_per_node
                cores_per_socket = info.cores_per_socket
                sockets_per_node = info.sockets_per_node
                multi_node_communication = info.multi_node_communication
                memory = info.memory
                nodes = info.nodes
                threads_per_core = info.threads_per_core

                # Show
                print("    - " + fmt.bold + "number of NUMA domains per node: " + fmt.reset + str(numa_domains_per_node))
                print("    - " + fmt.bold + "number of nodes: " + fmt.reset + str(nodes))
                print("    - " + fmt.bold + "number of sockets per node: " + fmt.reset + str(sockets_per_node))
                print("    - " + fmt.bold + "number of cores per socket: " + fmt.reset + str(cores_per_socket))
                print("    - " + fmt.bold + "number of threads per core: " + fmt.reset + str(threads_per_core))
                print("    - " + fmt.bold + "suitable for multinode communication: " + fmt.reset + str(multi_node_communication))
                print("    - " + fmt.bold + "virtual memory per node: " + fmt.reset + str(memory))

        # Show status?
        if status: show_status(self.get_remote(host))

        print("")

    # -----------------------------------------------------------------

    def show_parallelizations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the host
        host = self.get_host_from_command(command, **kwargs)

        # Show
        self.show_parallelizations(host)

    # -----------------------------------------------------------------

    def show_parallelizations(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        # Get parallelizations
        parallelizations = self.get_parallelizations_for_host(host)

        print("")
        for parallelization in parallelizations:
            print(" - " + tostr(parallelization))
        print("")

    # -----------------------------------------------------------------

    def show_info_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show info
        self.show_simulation_info(simulation_name)

    # -----------------------------------------------------------------

    def show_simulation_info(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing info for simulation '" + simulation_name + "' ...")

        # Check whether info is defined
        if not self.has_info_for_simulation(simulation_name): print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ": no info")

        # Info is defined
        else:

            # Show simulation name
            print(fmt.blue + fmt.underlined + simulation_name + fmt.reset)
            print("")

            # Get the info
            info = self.get_info_for_simulation(simulation_name)

            # Show the info
            for name in info: print(" - " + fmt.bold + name + fmt.reset + ": " + tostr(info[name]))
            print("")

    # -----------------------------------------------------------------

    def mount_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Mount, if necessary
        if host_id not in self.mount_paths: self.mount_paths[host_id] = mount_remote(host_id)

        # Get paths
        remote_home_path = self.get_home_path(host_id)
        mount_path = self.mount_paths[host_id]

        # Return
        return mount_path, remote_home_path

    # -----------------------------------------------------------------

    def mount_remote_path(self, host_id, path):

        """
        This function ...
        :param host_id:
        :param path:
        :return:
        """

        # Get mount path
        mount_path, remote_home_path = self.mount_remote(host_id)

        # Determine path relative to the home directory
        relative_path = fs.relative_to(path, remote_home_path)

        # Set the mounted path
        mounted_path = fs.join(mount_path, relative_path)

        # Return the mounted path
        return mounted_path

    # -----------------------------------------------------------------

    def mount_and_open_path(self, host_id, path):

        """
        This function ...
        :param host_id:
        :param path:
        :return:
        """

        # Get mount path
        output_path = self.mount_remote_path(host_id, path)

        # Open
        fs.open_directory(output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def open_base_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "open on remote filesystem (determined automatically based on the simulation status by default)", None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def open_base_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _open_command_name + " " + _base_command_name

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.open_base_definition, name=name, **kwargs)

        # Open
        self.open_base(simulation_name, remote=config.remote)

    # -----------------------------------------------------------------

    def open_base(self, simulation_name, remote=None):

        """
        This function ...
        :param simulation_name:
        :param remote:
        :return:
        """

        # Debugging
        log.debug("Opening the base directory of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set remote flag
        if remote is None: remote = False

        # On the remote
        if remote: self.mount_and_open_path(simulation.host_id, self.get_remote_base_path(simulation_name))

        # Locally
        else: fs.open_directory(self.get_base_path(simulation_name))

    # -----------------------------------------------------------------

    @lazyproperty
    def open_input_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "open on remote filesystem (determined automatically based on the simulation status by default)", None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def open_input_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set name
        name = _open_command_name + " " + _input_command_name

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.open_input_definition, name=name, **kwargs)

        # Open
        self.open_input(simulation_name, remote=config.remote)

    # -----------------------------------------------------------------

    def open_input(self, simulation_name, remote=None):

        """
        This function ...
        :param simulation_name:
        :param remote:
        :return:
        """

        # Debugging
        log.debug("Opening the input directory of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Check whether the simulation has input
        if not simulation.has_input:
            log.warning("The simulation '" + simulation_name + "' has no input: skipping ...")
            return

        # Set remote flag
        if remote is None: remote = False

        # On the remote
        if remote: self.mount_and_open_path(simulation.host_id, simulation.remote_input_path)

        # Local, single output directory
        elif simulation.has_input_directory: fs.open_directory(simulation.input_path)

        # Local, input files in the same directory
        elif self.get_input(simulation_name).has_single_directory: fs.open_directory(self.get_input(simulation_name).single_directory_path)

        # Local, input files in different directories (scattered)
        else:
            log.warning("The input files for simulation '" + simulation_name + "' are in different directories: opening all directories ...")
            #directory_paths = self.get_input(simulation_name).directory_paths
            #for name in directory_paths:
            #    path = directory_paths[name]
            paths = self.get_input(simulation_name).paths
            for name in paths:
                fs.show_file_in_directory(paths[name])
                time.wait(1) # wait one second to avoid errors

    # -----------------------------------------------------------------

    @lazyproperty
    def open_output_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "open on remote filesystem (determined automatically based on the simulation status by default)", None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def open_output_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set name
        name = _open_command_name + " " + _output_command_name

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.open_output_definition, name=name, **kwargs)

        # Open
        self.open_output(simulation_name, remote=config.remote)

    # -----------------------------------------------------------------

    def open_output(self, simulation_name, remote=None):

        """
        This function ...
        :param simulation_name:
        :param remote:
        :return:
        """

        # Debugging
        log.debug("Opening the output directory of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set remote flag
        if remote is None:
            if self.is_retrieved(simulation_name): remote = False
            elif self.is_running_or_finished(simulation_name): remote = True
            elif self.is_crashed(simulation_name) or self.is_aborted(simulation_name):
                log.warning("The simulation '" + simulation_name + "' is aborted or has crashed, the output may be partial or missing")
                remote = True
            else:
                log.warning("The simulation '" + simulation_name + "' has no output yet, opening local output directory ...")
                remote = False

        # On the remote
        if remote: self.mount_and_open_path(simulation.host_id, simulation.remote_output_path)

        # Locally
        else: fs.open_directory(self.get_output_path(simulation_name))

    # -----------------------------------------------------------------

    def open_extraction_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _open_command_name + " " + _extraction_command_name

        # Parse the command
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Open
        self.open_extraction(simulation_name)

    # -----------------------------------------------------------------

    def open_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Opening the extraction output directory of simulation '" + simulation_name + "' ...")

        # Get the path
        path = self.get_simulation_extraction_path(simulation_name)

        # Open
        fs.open_directory(path)

    # -----------------------------------------------------------------

    def open_plotting_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _open_command_name + " " + _plotting_command_name

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Open
        self.open_plotting(simulation_name)

    # -----------------------------------------------------------------

    def open_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Opening the plotting output directory of simulation '" + simulation_name + "' ...")

        # Get the path
        path = self.get_simulation_plotting_path(simulation_name)

        # Open
        fs.open_directory(path)

    # -----------------------------------------------------------------

    def open_misc_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _open_command_name + " " + _misc_command_name

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Open
        self.open_misc(simulation_name)

    # -----------------------------------------------------------------

    def open_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Opening the misc output directory of simulation '" + simulation_name + "' ...")

        # Get the path
        path = self.get_simulation_misc_path(simulation_name)

        # Open
        fs.open_directory(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_seds(self):

        """
        This function ...
        :return:
        """

        # Initialize
        seds = OrderedDict()
        if self.config.reference_seds is None: return seds

        # Loop over the SEDs
        for name in self.config.reference_seds:

            # Get the path
            path = self.config.reference_seds[name]

            # Load the sed
            sed = ObservedSED.from_file(path)

            # Add
            seds[name] = sed

        # Return
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_seds_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("from_file", "use SED plot from file", None)
        definition.add_optional("path", "string", "destination filepath")
        return definition

    # -----------------------------------------------------------------

    def plot_seds_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.plot_seds_definition, **kwargs)

        # Plot
        self.plot_simulation_seds(simulation_name, path=config.path, from_file=config.from_file)

    # -----------------------------------------------------------------

    def plot_simulation_seds(self, simulation_name, path=None, from_file=None):

        """
        This function ...
        :param simulation_name:
        :param path:
        :param from_file:
        :return:
        """

        # Debugging
        log.debug("Plotting simulated SEDs for simulation '" + simulation_name + "' ...")

        # Check whether plot file is present
        if from_file is None: from_file = self.has_sed_plots_for_simulation(simulation_name)
        elif from_file and not self.has_sed_plots_for_simulation(simulation_name): raise IOError("Plot file not present")

        # Use file plot
        if from_file: self.get_sed_plots_for_simulation(simulation_name, path=path)

        # Make plots
        else: self.make_sed_plots_for_simulation(simulation_name, path=path)

    # -----------------------------------------------------------------

    def get_sed_plots_for_simulation(self, simulation_name, path=None):

        """
        This function ...
        :param simulation_name:
        :param path:
        :return:
        """

        # Get the paths
        filepaths = self.get_sed_plot_paths_for_simulation(simulation_name)
        filepath = strings.get_shortest(filepaths)

        # Copy plot file to destination?
        if path is not None: fs.copy_file(filepath, path)

        # No destination path given, show the plot
        else: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def make_sed_plots_for_simulation(self, simulation_name, path=None):

        """
        This function ...
        :param simulation_name:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Making plot ...")

        # Get simulation prefix
        prefix = self.get_simulation_prefix(simulation_name)

        # Get simulation output
        output = self.get_output(simulation_name)

        # Create SED plotter
        plotter = SEDPlotter()

        # Add reference SEDs
        for name in self.reference_seds:
            sed = self.reference_seds[name]
            plotter.add_sed(sed, label=name)

        # Plot the SEDs
        for path in output.seds:

            # Get instrument name
            instr_name = get_sed_instrument_name(path, prefix)

            # Load SEDs
            seds = load_multiple_seds(path, as_dict=True)

            # Loop over the SEDs
            for contribution in seds:

                # Determine label
                label = instr_name + " " + contribution

                # Get the SED
                sed = seds[contribution]

                # Add the SED to the plotter
                plotter.add_sed(sed, label=label)

        # BECAUSE FOR SOME REASON INTERACTIVE PLOTTING IS NOT WORKING

        # Set filepath
        if path is not None: filepath = path
        else: filepath = fs.join(introspection.pts_temp_dir, "seds.pdf")

        # Run the plotter
        plotter.run(output=filepath)

        # Open the file, if path was not specified
        if path is None: fs.open_file(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_datacubes_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        contributions = ["total", "direct", "scattered", "dust", "dustscattered", "transparent"]
        default_contributions = ["total"]
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("contributions", "string_list", "contributions for which to plot the datacubes", default_contributions, choices=contributions)
        definition.add_optional("instruments", "string_list", "instruments for which to plot the datacubes")

        # Return
        return definition

    # -----------------------------------------------------------------

    def plot_datacubes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.plot_datacubes_definition, **kwargs)

        # Plot
        self.plot_simulation_datacubes(simulation_name, contributions=config.contributions, instruments=config.instruments)

    # -----------------------------------------------------------------

    def plot_simulation_datacubes(self, simulation_name, contributions=None, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param contributions:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting simulated datacubes for simulation '" + simulation_name + "' ...")

        # Total
        if contributions is None or "total" in contributions: self.plot_total_datacubes(simulation_name, instruments=instruments)

        # Direct
        if contributions is None or "direct" in contributions: self.plot_direct_datacubes(simulation_name, instruments=instruments)

        # Transparent
        if contributions is None or "transparent" in contributions: self.plot_transparent_datacubes(simulation_name, instruments=instruments)

        # Scattered
        if contributions is None or "scattered" in contributions: self.plot_scattered_datacubes(simulation_name, instruments=instruments)

        # Dust
        if contributions is None or "dust" in contributions: self.plot_dust_datacubes(simulation_name, instruments=instruments)

        # Dust-scattered
        if contributions is None or "dustscattered" in contributions: self.plot_dustscattered_datacubes(simulation_name, instruments=instruments)

    # -----------------------------------------------------------------

    def _plot_datacubes(self, simulation_name, paths, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param paths:
        :param instruments:
        :return:
        """

        from ...magic.tools import plotting
        from ...magic.core.datacube import DataCube

        # Get simulation prefix
        prefix = self.get_simulation_prefix(simulation_name)

        # Get the wavelength grid
        wavelength_grid = self.get_wavelength_grid(simulation_name)

        # Show the datacubes
        for path in paths:

            # Get instrument name
            instr_name = get_datacube_instrument_name(path, prefix)
            if instruments is not None and instr_name not in instruments: continue

            # Load the datacube
            datacube = DataCube.from_file(path, wavelength_grid)

            # Plot
            plotting.plot_datacube(datacube, title=instr_name, share_normalization=False, show_axes=False)

    # -----------------------------------------------------------------

    def plot_total_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Showing total datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.total_images, instruments=instruments)

    # -----------------------------------------------------------------

    def plot_direct_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting direct datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.direct_images, instruments=instruments)

    # -----------------------------------------------------------------

    def plot_transparent_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting transparent datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.transparent_images, instruments=instruments)

    # -----------------------------------------------------------------

    def plot_scattered_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting scattered datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.scattered_images, instruments=instruments)

    # -----------------------------------------------------------------

    def plot_dust_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting dust datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.dust_images, instruments=instruments)

    # -----------------------------------------------------------------

    def plot_dustscattered_datacubes(self, simulation_name, instruments=None):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :return:
        """

        # Debugging
        log.debug("Plotting dust-scattered datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Show
        self._plot_datacubes(simulation_name, output.dust_scattered_images, instruments=instruments)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_fluxes_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("from_images", "plot the fluxes calculated from images (by default, determined automatically)", None)
        definition.add_flag("from_file", "use an already existing plot file (by default, determined automatically)", None)
        definition.add_optional("path", "string", "destination filepath")
        return definition

    # -----------------------------------------------------------------

    def plot_fluxes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.plot_fluxes_definition, **kwargs)

        # Determine flag
        if config.from_images is None:

            config.pop("from_images")

            # Plot already present?
            if self.has_fluxes_plot_for_simulation(simulation_name): from_images = False
            elif self.has_image_fluxes_plot_for_simulation(simulation_name): from_images = True

            # Plot not yet present, but fluxes are present?
            elif self.has_fluxes_for_simulation(simulation_name): from_images = False
            elif self.has_image_fluxes_for_simulation(simulation_name): from_images = True

            # Neither is present
            else: raise RuntimeError("No fluxes for simulation '" + simulation_name + "'")

        else: from_images = config.pop("from_images")

        # From images
        if from_images: self.plot_simulation_fluxes_from_images(simulation_name, **config)

        # Regular
        else: self.plot_simulation_fluxes(simulation_name, **config)

    # -----------------------------------------------------------------

    def plot_simulation_fluxes(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Plotting mock fluxes for simulation '" + simulation_name + "' ...")

        # Check whether plot is already present
        from_plot_file = kwargs.pop("from_file", None)
        if from_plot_file is None: from_plot_file = self.has_fluxes_plot_for_simulation(simulation_name)
        elif from_plot_file and not self.has_fluxes_plot_for_simulation(simulation_name): raise IOError("Plot file not present")

        # Use file plot
        if from_plot_file: self.get_fluxes_plot_for_simulation(simulation_name, **kwargs)

        # Make the plot
        else: self.make_fluxes_plot_for_simulation(simulation_name, **kwargs)

    # -----------------------------------------------------------------

    def get_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the plot path
        filepath = self.get_fluxes_plot_path_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Copy plot file to destination?
        if path is not None: fs.copy_file(filepath, path)

        # No destination path given, show the plot
        else: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def make_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Making plot ...")

        # Get the fluxes
        fluxes = self.get_fluxes_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Create SED plotter
        plotter = SEDPlotter()

        # Add reference SEDs
        for name in self.reference_seds:
            sed = self.reference_seds[name]
            plotter.add_sed(sed, label=name)

        # Add the mock SED to the plotter
        plotter.add_sed(fluxes, label="Mock observation")

        # BECAUSE FOR SOME REASON INTERACTIVE PLOTTING IS NOT WORKING

        # Set filepath
        if path is not None: filepath = path
        else: filepath = fs.join(introspection.pts_temp_dir, "fluxes.pdf")

        # Run the plotter
        plotter.run(output=filepath)

        # Open the file, if path was not specified
        if path is None: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def plot_simulation_fluxes_from_images(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Plotting mock fluxes from images for simulation '" + simulation_name + "' ...")

        # Check whether plot is already present
        from_plot_file = kwargs.pop("from_file", None)
        if from_plot_file is None: from_plot_file = self.has_image_fluxes_plot_for_simulation(simulation_name)
        elif from_plot_file and not self.has_image_fluxes_plot_for_simulation(simulation_name): raise IOError("Plot file not present")

        # Use file plot
        if from_plot_file: self.get_image_fluxes_plot_for_simulation(simulation_name, **kwargs)

        # Make the plot
        else: self.make_image_fluxes_plot_for_simulation(simulation_name, **kwargs)

    # -----------------------------------------------------------------

    def get_image_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the plot path
        filepath = self.get_image_fluxes_plot_path_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Copy plot file to destination?
        if path is not None: fs.copy_file(filepath, path)

        # No destination path given, show the plot
        else: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def make_image_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Making plot ...")

        # Get the fluxes
        fluxes = self.get_image_fluxes_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Create SED plotter
        plotter = SEDPlotter()

        # Add reference SEDs
        for name in self.reference_seds:
            sed = self.reference_seds[name]
            plotter.add_sed(sed, label=name)

        # Add the mock SED to the plotter
        plotter.add_sed(fluxes, label="Mock observation")

        # BECAUSE FOR SOME REASON INTERACTIVE PLOTTING IS NOT WORKING

        # Set filepath
        if path is not None: filepath = path
        else: filepath = fs.join(introspection.pts_temp_dir, "fluxes.pdf")

        # Run the plotter
        plotter.run(output=filepath)

        # Open the file, if path was not specified
        if path is None: fs.open_file(filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_images_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("for_fluxes", "plot the images created for calculating fluxes (by default, determined automatically)", None)
        definition.add_flag("from_file", "use an already existing plot file (by default, determined automatically)", None)
        definition.add_optional("path", "string", "destination filepath")
        return definition

    # -----------------------------------------------------------------

    def plot_images_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.plot_images_definition, **kwargs)

        # Determine flag
        if config.for_fluxes is None:

            config.pop("for_fluxes")

            # Check if plots are already created
            if self.has_images_plot_for_simulation(simulation_name): for_fluxes = False
            elif self.has_images_for_fluxes_plot_for_simulation(simulation_name): for_fluxes = True

            # Check if images are present
            elif self.has_images_for_simulation(simulation_name): for_fluxes = False
            elif self.has_images_for_fluxes_for_simulation(simulation_name): for_fluxes = True

            # Neither is present
            else: raise RuntimeError("No images for simulation '" + simulation_name + "'")

        # Flag is passed
        else: for_fluxes = config.pop("for_fluxes")

        # For fluxes
        if for_fluxes: self.plot_simulation_images_for_fluxes(simulation_name, **config)

        # Regular
        else: self.plot_simulation_images(simulation_name, **config)

    # -----------------------------------------------------------------

    def plot_simulation_images(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Plotting mock images for simulation '" + simulation_name + "' ...")

        # Check whether the plot is already present
        from_plot_file = kwargs.pop("from_file", None)
        if from_plot_file is None: from_plot_file = self.has_images_plot_for_simulation(simulation_name)
        elif from_plot_file and not self.has_images_plot_for_simulation(simulation_name): raise IOError("Plot file not present")

        # Use file plot
        if from_plot_file: self.get_images_plot_for_simulation(simulation_name, **kwargs)

        # Make plot
        else: self.make_images_plot_for_simulation(simulation_name, **kwargs)

    # -----------------------------------------------------------------

    def get_images_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the plot path
        filepath = self.get_images_plot_path_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Copy plot file to destination?
        if path is not None: fs.copy_file(filepath, path)

        # No destination path given, show the plot
        else: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def make_images_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the images
        images = self.get_images_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    def plot_simulation_images_for_fluxes(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Plotting mock images (for fluxes) for simulation '" + simulation_name + "' ...")

        # Check whether the plot is already present
        from_plot_file = kwargs.pop("from_file", None)
        if from_plot_file is None: from_plot_file = self.has_images_for_fluxes_plot_for_simulation(simulation_name)
        elif from_plot_file and not self.has_images_for_fluxes_plot_for_simulation(simulation_name): raise IOError("Plot file not present")

        # Use file plot
        if from_plot_file: self.get_images_for_fluxes_plot_for_simulation(simulation_name, **kwargs)

        # Make plot
        else: self.make_images_for_fluxes_plot_for_simulation(simulation_name, **kwargs)

    # -----------------------------------------------------------------

    def get_images_for_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the plot path
        filepath = self.get_images_for_fluxes_plot_path_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

        # Copy plot file to destination?
        if path is not None: fs.copy_file(filepath, path)

        # No destination path given, show the plot
        else: fs.open_file(filepath)

    # -----------------------------------------------------------------

    def make_images_for_fluxes_plot_for_simulation(self, simulation_name, **kwargs):

        """
        This function ...
        :param simulation_name:
        :param kwargs:
        :return:
        """

        # Get the images
        images = self.get_images_for_fluxes_for_simulation(simulation_name)

        # Get destination path
        path = kwargs.pop("path", None)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_input_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "show remote input", False)
        return definition

    # -----------------------------------------------------------------

    def show_input_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_input_definition, **kwargs)

        # Show
        if config.remote: self.show_input_remote(simulation_name)
        else: self.show_input_local(simulation_name)

    # -----------------------------------------------------------------

    def show_input_local(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing local input files of simulation '" + simulation_name + "' ...")

        # Get the input
        input = self.get_input(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the input
        input.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_input_remote(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing remote input files of simulation '" + simulation_name + "' ...")

        # Get the input
        input = self.get_input(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the input
        input.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_output_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "show remote output", None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def show_output_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_output_definition, **kwargs)

        # Set remote flag
        if config.remote is None:
            if self.is_retrieved(simulation_name): remote = False
            elif self.is_running_or_finished(simulation_name): remote = True
            elif self.is_crashed(simulation_name) or self.is_aborted(simulation_name):
                log.warning("The simulation '" + simulation_name + "' is aborted or has crashed, the output may be partial or missing")
                remote = True
            else:
                log.warning("The simulation '" + simulation_name + "' has no output yet, skipping ...")
                return
        else: remote = config.remote

        # Show
        if remote: self.show_output_remote(simulation_name)
        else: self.show_output_local(simulation_name)

    # -----------------------------------------------------------------

    def show_output_local(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing local output files of simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the output
        output.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_output_remote(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing remote output files of simulation '" + simulation_name + "' ...")

        # Get remote simulation output
        output = self.get_remote_output(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the output
        output.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_extraction_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_extraction(simulation_name)

    # -----------------------------------------------------------------

    def show_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing extraction output of simulation '" + simulation_name + "' ...")

        # Get extraction output
        extraction = self.get_extraction_output(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the output
        extraction.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_plotting_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_plotting(simulation_name)

    # -----------------------------------------------------------------

    def show_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing plotting output of simulation '" + simulation_name + "' ...")

        # Get plotting output
        plotting = self.get_plotting_output(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the output
        plotting.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_misc_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_misc(simulation_name)

    # -----------------------------------------------------------------

    def show_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing miscellaneous output of simulation '" + simulation_name + "' ...")

        # Get misc output
        misc = self.get_misc_output(simulation_name)

        # Print the simulation name
        print(fmt.blue + fmt.underlined + simulation_name + fmt.reset + ":")
        print("")

        # Show the output
        misc.show(line_prefix="  ")
        print("")

    # -----------------------------------------------------------------

    def show_instruments_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show instruments
        self.show_instruments(simulation_name)

    # -----------------------------------------------------------------

    def show_instruments(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing instruments of simulation '" + simulation_name + "' ...")

        print("")
        # Loop over the instruments
        ski = self.get_skifile(simulation_name)
        for name in self.get_instrument_names(simulation_name): show_instrument(ski, name)
        print("")

    # -----------------------------------------------------------------

    def show_stellar_components_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_stellar_components(simulation_name)

    # -----------------------------------------------------------------

    def show_stellar_components(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing stellar components of simulation '" + simulation_name + "' ...")

        # Loop over the stellar components
        ski = self.get_skifile(simulation_name)
        for id in self.get_stellar_component_ids(simulation_name): show_stellar_component(ski, id)
        print("")

    # -----------------------------------------------------------------

    def show_dust_components_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_dust_components(simulation_name)

    # -----------------------------------------------------------------

    def show_dust_components(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Showing dust components of simulation '" + simulation_name + "' ...")

        # Loop over the dust components
        ski = self.get_skifile(simulation_name)
        for id in self.get_dust_component_ids(simulation_name): show_dust_component(ski, id)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_normalizations_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("flux", "photometric_unit", "also show flux in a particular unit")

        # Return
        return definition

    # -----------------------------------------------------------------

    def show_normalizations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_normalizations_definition, **kwargs)

        # Show
        self.show_normalizations(simulation_name, flux_unit=config.flux)

    # -----------------------------------------------------------------

    def show_normalizations(self, simulation_name, flux_unit=None):

        """
        This function ...
        :param simulation_name:
        :param flux_unit:
        :return:
        """

        # Debugging
        log.debug("Showing normalizations of simulation '" + simulation_name + "' ...")

        # Get the ski file
        ski = self.get_skifile(simulation_name)

        # Show
        show_normalizations(ski, flux_unit=flux_unit)

    # -----------------------------------------------------------------

    def get_config_from_command(self, command, definition, name=None, index=1, interactive=False):

        """
        This function ...
        :param command:
        :param definition:
        :param name:
        :param index:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        parse_command = splitted[index:]

        # Interactively get the settings
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False,
                                 initialize=False, add_logging=False, add_cwd=False)

        # Return the configuration
        return config

    # -----------------------------------------------------------------

    def get_usage_command(self, definition, name):

        """
        This function ...
        :param definition:
        :param name:
        :return:
        """

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_command(self, definition, name):

        """
        This function ...
        :param definition:
        :param name:
        :return:
        """

        # Return the usage
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_host_command_definition(self, command_definition=None, required=True, choices=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required:
        :param choices:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Add host setting
        if required: definition.add_required("host", "host", "remote host", choices=choices)
        else: definition.add_positional_optional("host", "host", "remote host", choices=choices)

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_host_command(self, command, command_definition=None, name=None, index=1, required=True, choices=None,
                           required_to_optional=True, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required:
        :param choices:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        parse_command = splitted[index:]

        # Get the definition
        definition = self.get_host_command_definition(command_definition, required=required, choices=choices, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get the host
        host = config.pop("host")

        # Return
        return splitted, host, config

    # -----------------------------------------------------------------

    def get_usage_host_command(self, command_definition=None, name=None, required=True, choices=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required:
        :param choices:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_host_command_definition(command_definition, required=required, choices=choices, required_to_optional=required_to_optional)

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_host_command(self, command_definition=None, name=None, required=True, choices=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required:
        :param choices:
        :param required_to_optional:
        :return:
        """

        # Get the efinition
        definition = self.get_host_command_definition(command_definition, required=required, choices=choices, required_to_optional=required_to_optional)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_host_from_command(self, command, name=None, index=1, required=True, choices=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :param required:
        :param choices:
        :param interactive:
        :return:
        """

        # Parse
        splitted, host, config = self.parse_host_command(command, name=name, index=index, required=required, choices=choices, interactive=interactive)

        # Return the host
        return host

    # -----------------------------------------------------------------

    def get_hosts_command_definition(self, command_definition=None, required=False, choices=None):

        """
        This function ...
        :param command_definition:
        :param required:
        :param choices:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Hosts setting
        if required: definition.add_required("hosts", "host_list", "remote hosts", choices=choices)
        else: definition.add_positional_optional("hosts", "host_list", "remote hosts", choices=choices)

        # Add definition settings
        if command_definition is not None: definition.import_settings(command_definition, required_to="optional")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_hosts_command(self, command, command_definition=None, name=None, index=1, required=False, choices=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required:
        :param choices:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only host list

        # Get the definition
        definition = self.get_hosts_command_definition(command_definition, required=required, choices=choices)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Return
        hosts = config.pop("hosts")
        return splitted, hosts, config

    # -----------------------------------------------------------------

    def get_usage_hosts_command(self, command_definition=None, name=None, required=False, choices=None):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required:
        :param choices:
        :return:
        """

        # Get the definition
        definition = self.get_hosts_command_definition(command_definition, required=required, choices=choices)

        # Return usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_hosts_command(self, command_definition=None, name=None, required=False, choices=None):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required:
        :param choices:
        :return:
        """

        # Get the definition
        definition = self.get_hosts_command_definition(command_definition, required=required, choices=choices)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_hosts_from_command(self, command, name=None, index=1, required=False, choices=None):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :param required:
        :param choices:
        :return:
        """

        # Parse command
        splitted, hosts, config = self.parse_hosts_command(command, name=name, index=index, required=required, choices=choices)

        # Return
        return hosts

    # -----------------------------------------------------------------

    def get_host_and_parallelization_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Host and parallelization
        definition.add_positional_optional("host", "host", "remote host")
        definition.add_positional_optional("parallelization", "parallelization", "parallelization scheme")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_host_and_parallelization_command(self, command, command_definition=None, name=None, index=1,
                                               required_to_optional=True, interactive=False):

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
        parse_command = splitted[index:]

        # Get the definition
        definition = self.get_host_and_parallelization_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get host and parallelization
        host = config.pop("host")
        parallelization = config.pop("parallelization")

        # Return
        return splitted, host, parallelization, config

    # -----------------------------------------------------------------

    def get_usage_host_and_parallelization_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_host_and_parallelization_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_host_and_parallelization_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_host_and_parallelization_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_host_and_parallelization_from_command(self, command, name=None, index=1, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :param interactive:
        :return:
        """

        # Parse
        splitted, host, parallelization, config = self.parse_host_and_parallelization_command(command, name=name, index=index, interactive=interactive)

        # Return
        return host, parallelization

    # -----------------------------------------------------------------

    def get_simulation_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("simulation", "integer_or_string", "simulation index or name")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_simulation_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True,
                                 interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :param command_definition:
        :param required_to_optional:
        :param interactive:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        if command_definition is not None: parse_command = splitted[index:]
        else: parse_command = splitted[index:index + 1]  # only simulation name

        # Get the definition
        definition = self.get_simulation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation name
        if types.is_integer_type(config.simulation): simulation_name = self.simulation_names[config.pop("simulation")]
        else: simulation_name = config.pop("simulation")

        # Return
        return splitted, simulation_name, config

    # -----------------------------------------------------------------

    def get_usage_simulation_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_simulation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_simulation_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_simulation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_simulations_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("simulations", "integer_list", "simulation indices")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_simulations_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

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
        else: parse_command = splitted[index:index + 1]  # only simulation indices

        # Get definition
        definition = self.get_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation names
        simulation_names = []
        for index in config.simulations: simulation_names.append(self.simulation_names[index])

        # Return
        return splitted, simulation_names, config

    # -----------------------------------------------------------------

    def get_usage_simulations_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_simulations_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_two_simulations_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("simulation_a", "integer_or_string", "simulation index or name")
        definition.add_required("simulation_b", "integer_or_string", "simulation index or name")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_two_simulations_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):

        """
        This function ....
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

        # Get the definition
        definition = self.get_two_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=splitted[index:], error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation_a name
        if types.is_integer_type(config.simulation_a): simulation_a_name = self.simulation_names[config.pop("simulation_a")]
        else: simulation_a_name = config.pop("simulation_a")

        # Get simulation_b name
        if types.is_integer_type(config.simulation_b): simulation_b_name = self.simulation_names[config.pop("simulation_b")]
        else: simulation_b_name = config.pop("simulation_b")

        # Return
        return splitted, simulation_a_name, simulation_b_name, config

    # -----------------------------------------------------------------

    def get_usage_two_simulations_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_two_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_two_simulations_command(self, command_definition=None, name=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param name:
        :param required_to_optional:
        :return:
        """

        # Get the definition
        definition = self.get_two_simulations_command_definition(command_definition, required_to_optional=required_to_optional)

        # Return the help info
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_simulation_names_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, simulation_names, config = self.parse_simulations_command(command, name=name, interactive=interactive)

        # Return the simulation names
        return simulation_names

    # -----------------------------------------------------------------

    def get_simulation_name_from_command(self, command, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, name=name, interactive=interactive)

        # Return the simulation name
        return simulation_name

    # -----------------------------------------------------------------

    def get_simulation_name_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, command_definition=command_definition, name=name, interactive=interactive)

        # Return simulation name and config
        return simulation_name, config

    # -----------------------------------------------------------------

    @lazyproperty
    def show_simulation_log_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("summarize", "show summarized log output")
        definition.add_flag("short", "show short output")
        definition.add_flag("screen", "show screen log output instead of SKIRT log file")
        definition.add_flag("job", "show job log output instead of SKIRT log file")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def show_simulation_log_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_simulation_log_definition, **kwargs)

        # Check
        if not self.is_running_or_finished_or_aborted_or_crashed(simulation_name): raise ValueError("Simulation '" + simulation_name + "' cannot have log output (yet)")
        if self.is_running(simulation_name): log.warning("Simulation '" + simulation_name + "' is still running")
        elif self.is_failed(simulation_name): log.warning("Simulation '" + simulation_name + "' has not finished succesfully")

        # Show screen log
        if config.screen:

            # Check
            if not self.is_screen_execution(simulation_name): raise ValueError("Simulation '" + simulation_name + "' is/was not executed in a screen session")

            # Show output
            self.show_simulation_screen_output(simulation_name, summarize=config.summarize, short=config.short)

        # Show job log
        elif config.job:

            # Check
            if not self.is_job_execution(simulation_name): raise ValueError("Simulation '" + simulation_name + "' is/was not executed as a job")

            # Show output
            self.show_simulation_job_output(simulation_name, summarize=config.summarize, short=config.short)

        # Show simulation log
        else: self.show_simulation_log(simulation_name, summarize=config.summarize, short=config.short)

    # -----------------------------------------------------------------

    def show_simulation_log(self, simulation_name, summarize=False, short=False):

        """
        This function ...
        :param simulation_name:
        :param summarize:
        :param short:
        :return:
        """

        # Debugging
        log.debug("Showing log output of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Simulation is retrieved
        if simulation.retrieved:

            # Determine the path to the simulation log file
            local_log_file_path = simulation.log_file_path

            # Read the log file
            lines = fs.read_lines(local_log_file_path)

        # Not yet retrieved
        else:

            # Get the remote
            remote = self.get_remote(simulation.host_id)

            # The path to the simulation log file
            remote_log_file_path = simulation.remote_log_file_path

            # Check whether the log file exists
            if not remote.is_file(remote_log_file_path): raise RuntimeError("The log file does not (yet) exist remotely")

            # Read the log file
            lines = remote.read_lines(remote_log_file_path)

        # Print the lines of the log file
        print("")
        if summarize: show_log_summary(lines, debug_output=True)
        elif short: show_log_summary(lines, debug_output=False)
        else:
            for line in lines: print(line)
        print("")

    # -----------------------------------------------------------------

    def show_simulation_screen_output(self, simulation_name, summarize=False, short=False):

        """
        This function ...
        :param simulation_name:
        :param summarize:
        :param short:
        :return:
        """

        # Get the screen name
        screen_name = self.get_screen_name(simulation_name)

        # Debugging
        log.debug("Showing output of screen session '" + screen_name + "' for simulation '" + simulation_name + "' ...")

        # Get the screenlog path
        filepath = self.get_screen_log_path(simulation_name)
        #print(filepath)

        remote_ski_path = self.get_remote_skifile_path(simulation_name)

        # Get the lines
        lines = self.get_remote_for_simulation(simulation_name).get_lines_between(filepath, remote_ski_path, "Finished simulation")
        
        # Show
        print("")
        if summarize: show_log_summary(lines, debug_output=True)
        elif short: show_log_summary(lines, debug_output=False)
        else:
            for line in lines:

                # Show the line
                if not has_valid_timestamp(line): print(fmt.red + line + fmt.reset)
                else: print(line)

        print("")

    # -----------------------------------------------------------------

    def show_simulation_job_output(self, simulation_name, summarize=False, short=False):

        """
        This function ...
        :param simulation_name:
        :param summarize:
        :param short:
        :return:
        """

        # Get the job ID
        job_id = self.get_job_id(simulation_name)

        # Debugging
        log.debug("Showing output for job '" + str(job_id) + "' ...")

        # Get the output lines
        lines = self.get_job_output(simulation_name)
        nlines = len(lines)
        if nlines == 0:
            log.warning("No output for job '" + str(job_id) + "'")
            return

        # Show
        print("")
        if summarize: show_log_summary(lines, debug_output=True)
        elif short: show_log_summary(lines, debug_output=False)
        else:
            for line in lines:

                # Show the line
                if not has_valid_timestamp(line): print(fmt.red + line + fmt.reset)
                else: print(line)

        print("")

    # -----------------------------------------------------------------

    def show_simulation_errors_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_simulation_errors(simulation_name)

    # -----------------------------------------------------------------

    def show_simulation_errors(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Show error output of simulation '" + simulation_name + "'...")

        # Executed in screen
        if self.is_screen_execution(simulation_name): self.show_simulation_errors_screen(simulation_name)

        # Executed in job
        elif self.is_job_execution(simulation_name): self.show_simulation_errors_job(simulation_name)

        # Not supported
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    def show_simulation_errors_screen(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the screen name
        screen_name = self.get_screen_name(simulation_name)

        # Debugging
        log.debug("Showing error output of screen session '" + screen_name + "' for simulation '" + simulation_name + "' ...")

        # Get the screenlog path
        filepath = self.get_screen_log_path(simulation_name)
        print(filepath)

    # -----------------------------------------------------------------

    def show_simulation_errors_job(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the job ID
        job_id = self.get_job_id(simulation_name)

        # Debugging
        log.debug("Showing error output for job '" + str(job_id) + "' ...")

        # Get the error lines
        lines = self.get_job_error(simulation_name)
        nlines = len(lines)
        if nlines == 0:
            print(fmt.red + "no error output" + fmt.reset)
            return

        # Show
        print("")
        break_after = False
        for line in lines:

            # Skip regular SKIRT output messages
            if has_valid_timestamp(line): continue

            # Show the line
            print(fmt.red + line + fmt.reset)
            if break_after: break

            # Skip after certain MPI message
            if "YOU CAN IGNORE THE BELOW CLEANUP MESSAGES" in line: break_after = True

        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def show_simulation_settings_definition(self):

        """
        This function ...
        :return:
        """

        return show_simulation_definition

    # -----------------------------------------------------------------

    def show_simulation_settings_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_simulation_settings_definition, **kwargs)

        # Show
        self.show_simulation_settings(simulation_name, config=config)

    # -----------------------------------------------------------------

    def show_simulation_settings(self, simulation_name, config=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Showing settings of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Show
        show_simulation(simulation, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_analysis_options_definition(self):

        """
        This function ...
        :return:
        """

        return show_analysis_definition

    # -----------------------------------------------------------------

    def show_analysis_options_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.show_analysis_options_definition, **kwargs)

        # Show
        self.show_analysis_options(simulation_name, config=config)

    # -----------------------------------------------------------------

    def show_analysis_options(self, simulation_name, config=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Showing analysis options for simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Show
        show_analysis(simulation, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def adapt_simulation_settings_definition(self):

        """
        This function ...
        :return:
        """

        return adapt_simulations_definition

    # -----------------------------------------------------------------

    def adapt_simulation_settings_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _adapt_command_name + " " + _simulation_command_name

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.adapt_simulation_settings_definition, name=name, **kwargs)

        # Adapt
        self.adapt_simulation_settings(simulation_name, config=config)

    # -----------------------------------------------------------------

    def adapt_simulation_settings(self, simulation_name, config=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Adapting simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Adapt
        adapt_simulation(simulation, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def adapt_analysis_options_definition(self):

        """
        This function ...
        :return:
        """

        return adapt_analysis_definition

    # -----------------------------------------------------------------

    def adapt_analysis_options_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _adapt_command_name + " " + _analysis_command_name

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.adapt_analysis_options_definition, name=name, **kwargs)

        # Adapt
        self.adapt_analysis_options(simulation_name, config=config)

    # -----------------------------------------------------------------

    def adapt_analysis_options(self, simulation_name, config=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Adapting analysis options of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Adapt
        adapt_analysis(simulation, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def compare_simulation_settings_definition(self):

        """
        This function ...
        :return:
        """

        return show_simulation_definition

    # -----------------------------------------------------------------

    def compare_simulation_settings_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _compare_command_name + " " + _simulation_command_name

        # Parse
        splitted, simulation_a_name, simulation_b_name, config = self.parse_two_simulations_command(command, self.compare_simulation_settings_definition, name=name, **kwargs)

        # Compare
        self.compare_simulation_settings(simulation_a_name, simulation_b_name, config=config)

    # -----------------------------------------------------------------

    def compare_simulation_settings(self, simulation_a_name, simulation_b_name, config=None):

        """
        This function ...
        :param simulation_a_name:
        :param simulation_b_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Comparing simulations '" + simulation_a_name + "' and '" + simulation_b_name + "' ...")

        # Get the simulations
        simulation_a = self.get_simulation(simulation_a_name)
        simulation_b = self.get_simulation(simulation_b_name)

        # Compare
        compare_simulations(simulation_a, simulation_b, config=config)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def compare_analysis_options_definition(self):

        """
        This function ...
        :return:
        """

        return show_analysis_definition

    # -----------------------------------------------------------------

    def compare_analysis_options_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _compare_command_name + " " + _simulation_command_name

        # Parse
        splitted, simulation_a_name, simulation_b_name, config = self.parse_two_simulations_command(command, self.compare_analysis_options_definition, name=name, **kwargs)

        # Compare
        self.compare_analysis_options(simulation_a_name, simulation_b_name, config=config)

    # -----------------------------------------------------------------

    def compare_analysis_options(self, simulation_a_name, simulation_b_name, config=None):

        """
        This function ...
        :param simulation_a_name:
        :param simulation_b_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Comparing analysis options between simulations '" + simulation_a_name + "' and '" + simulation_b_name + "' ...")

        # Get the simulations
        simulation_a = self.get_simulation(simulation_a_name)
        simulation_b = self.get_simulation(simulation_b_name)

        # Compare
        compare_analysis(simulation_a, simulation_b, config=config)
        print("")

    # -----------------------------------------------------------------

    def move(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Moving simulations ...")

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Get the simulation
            simulation = self.get_simulation(simulation_name)

            # Check remote
            if self.config.move_remotes is not None and simulation.host_id not in self.config.move_remotes: continue

            # Check simulation name
            if self.config.move_simulations is not None and simulation_name not in self.config.move_simulations: continue

            # Get the status
            status = self.get_status(simulation_name)

            # Skip finished and running
            if is_finished_status(status): continue
            if not self.config.move_running and is_running_status(status): continue

            # Get display name
            display_name = self.get_display_name(simulation, id_size=3)

            # Move?
            print(" - " + display_name + ": " + status)
            if self.config.prompt_simulations_move and not prompt_proceed("move simulation " + display_name + "?"): continue

            # To which host?
            if self.config.move_to_host is not None: host = self.config.move_to_host
            else: host = prompt_variable("host", "host", "remote host to move simulation '" + simulation_name + "' to", choices=all_host_ids)
            if host == simulation.host: raise ValueError("Simulation '" + simulation_name + "' is already queued/running on host '" + tostr(host) + "'")

            # Move simulation
            self.move_simulation(simulation_name, host)

    # -----------------------------------------------------------------

    @lazyproperty
    def move_simulations_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("host", "host", "remote host to move the simulation to")
        definition.add_positional_optional("parallelization", "parallelization", "parallelization scheme for the simulation")
        definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def move_simulations_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def move_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set all kwargs
        kwargs.update(self.move_simulations_kwargs)

        # Get the simulation names
        splitted, simulation_names, config = self.parse_simulations_command(command, self.move_simulations_definition, **kwargs)

        # Loop over the simulation names
        for simulation_name in simulation_names:

            # Check status
            if self.is_finished(simulation_name): raise ValueError("Simulation is already finished")
            if self.is_running(simulation_name): log.warning("Simulation is already running")

            # Get host and check
            if config.host == self.get_simulation(simulation_name).host: raise ValueError("Simulation '" + simulation_name + "' is already queued/running on host '" + tostr(config.host) + "'")

            # Move simulation
            self.move_simulation(simulation_name, config.host, config.parallelization, config.logging, config.scheduling)

    # -----------------------------------------------------------------

    def move_simulation(self, simulation_name, host, parallelization=None, logging_options=None, scheduling_options=None):

        """
        This function ...
        :param simulation_name:
        :param host:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :return:
        """

        # Debugging
        log.debug("Moving simulation '" + simulation_name + "' to host '" + tostr(host) + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Cancel or stop the simulation
        if self.is_running(simulation_name): self.stop_simulation(simulation_name)
        elif self.is_queued(simulation_name): self.cancel_simulation(simulation_name)

        # TODO: support different clusters
        host_id = host.id

        # Add the simulation to the queue
        self.launcher.add_to_queue(simulation.definition, simulation_name, host_id=host_id, parallelization=parallelization,
                                   logging_options=logging_options, scheduling_options=scheduling_options,
                                   analysis_options=simulation.analysis)

        # Add entry to moved table
        self.moved.add_simulation(simulation, host)

        # Remove original simulation?
        if not self.config.dry: self.remove_simulation(simulation_name) # backup is done here

        # Don't remove
        else: log.warning("[DRY] Not removing the simulation from '" + simulation.path + "' ...")

    # -----------------------------------------------------------------

    @property
    def nmoved(self):

        """
        This function ...
        :return:
        """

        return len(self.moved)

    # -----------------------------------------------------------------

    @property
    def has_moved(self):

        """
        This function ...
        :return:
        """

        return self.nmoved > 0

    # -----------------------------------------------------------------

    @property
    def moved_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.moved.simulation_names

    # -----------------------------------------------------------------

    def is_moved(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.moved_simulation_names

    # -----------------------------------------------------------------

    def cancel_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation names
        simulation_names = self.get_simulation_names_from_command(command, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether queued
            if not self.is_queued(simulation_name): raise ValueError("The simulation '" + simulation_name + "' is not queued")

            # Cancel the simulation
            self.cancel_simulation(simulation_name)

    # -----------------------------------------------------------------

    def cancel_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Cancelling simulation '" + simulation_name + "' ...")

        # Screen
        if self.is_screen_execution(simulation_name): self.cancel_simulation_screen(simulation_name)

        # Job
        elif self.is_job_execution(simulation_name): self.cancel_simulation_job(simulation_name)

        # Not supported
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    def cancel_simulation_screen(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get screen name
        screen_name = self.get_screen_name(simulation_name)

        # Debugging
        log.debug("Cancelling simulation '" + simulation_name + "' queued in screen session '" + screen_name + "' ...")

        # Get the remote
        remote = self.get_remote_for_simulation(simulation_name)

        # Get the remote ski file path
        ski_path = self.get_remote_skifile_path(simulation_name)

        if remote.is_file(ski_path):

            # Debugging
            log.debug("Removing the ski file '" + ski_path + "' from remote host '" + remote.host_id + "' ...")

            # Remove the remote ski file
            remote.remove_file(ski_path)

        else: log.warning("The ski file '" + ski_path + "' is already removed from remote host '" + remote.host_id + "'")

    # -----------------------------------------------------------------

    def cancel_simulation_job(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get job ID
        job_id = self.get_job_id(simulation_name)

        # Debugging
        log.debug("Cancelling simulation '" + simulation_name + "' by deleting job '" + str(job_id) + "' ...")

        # Get the cluster name
        cluster_name = self.get_cluster_name_for_simulation(simulation_name)

        # Get the remote
        remote = self.get_remote_for_simulation(simulation_name)

        # Stop the job
        remote.kill_job(job_id, cluster_name=cluster_name, show_output=True)

        # Get the remote ski file path
        ski_path = self.get_remote_skifile_path(simulation_name)

        if remote.is_file(ski_path):

            # Debugging
            log.debug("Removing the ski file '" + ski_path + "' from remote host '" + remote.host_id + "' ...")

            # Remove the remote ski file
            remote.remove_file(ski_path)

        else: log.warning("The ski file '" + ski_path + "' is already removed from remote host '" + remote.host_id + "'")

    # -----------------------------------------------------------------

    def stop_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_names = self.get_simulation_names_from_command(command, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether running
            if not self.is_running(simulation_name): raise ValueError("The simulation '" + simulation_name + "' is not running")

            # Stop the simulation
            self.stop_simulation(simulation_name)

    # -----------------------------------------------------------------

    def stop_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Stopping simulation '" + simulation_name + "' ...")

        # Screen
        if self.is_screen_execution(simulation_name): self.stop_simulation_screen(simulation_name)

        # Job
        elif self.is_job_execution(simulation_name): self.stop_simulation_job(simulation_name)

        # Not supported
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    def stop_simulation_screen(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Give warning
        log.warning("Stopping a running simulation within a screen session is not supported")

    # -----------------------------------------------------------------

    def stop_simulation_job(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get job ID
        job_id = self.get_job_id(simulation_name)

        # Debugging
        log.debug("Stopping simulation by aborting job '" + str(job_id) + "'...")

        # Get the remote
        remote = self.get_remote_for_simulation(simulation_name)

        # Get the cluster name
        cluster_name = self.get_cluster_name_for_simulation(simulation_name)

        # Stop the job
        remote.kill_job(job_id, cluster_name=cluster_name, show_output=True)

    # -----------------------------------------------------------------

    def remove_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Remove
        self.remove_simulation(simulation_name)

    # -----------------------------------------------------------------

    def remove_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the simulation path
        path = self.get_simulation_path(simulation_name)

        # Make backup of the simulation file?
        if self.config.backup_simulations:
            host_id = self.get_simulation_host_id(simulation_name)
            log.debug("Making backup of the simulation '" + simulation_name + "' ...")
            fs.copy_file(path, self.get_backup_simulations_path_for_host(host_id))

        # Remove the simulation file
        fs.remove_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_simulation_input_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "clear on the remote host", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clear_simulation_input_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _clear_command_name + " " + _input_command_name

        # Parse command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.clear_simulation_input_definition, name=name, **kwargs)

        # Execute the command
        if config.analysis_steps is not None: raise ValueError("Cannot specify 'analysis_steps'")
        if config.remote: self.clear_simulation_input_remote(simulation_name)
        else: self.clear_simulation_input_local(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_simulation_output_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("remote", "clear on the remote host", False)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clear_simulation_output_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set the command name
        name = _clear_command_name + " " + _output_command_name

        # Parse command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.clear_simulation_output_definition, name=name, **kwargs)

        # Execute the command
        if config.analysis_steps is not None: raise ValueError("Cannot specify 'analysis_steps'")
        if config.remote: self.clear_simulation_output_remote(simulation_name)
        else: self.clear_simulation_output_local(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def clear_simulation_analysis_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("analysis_steps", "string_list", "analysis steps for which to clear the output", choices=clear_analysis_steps)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def clear_simulation_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _clear_command_name + " " + _analysis_command_name

        # Parse command
        splitted, simulation_name, config = self.parse_simulation_command(command, self.clear_simulation_analysis_definition, name=name, **kwargs)

        # Execute
        self.clear_simulation_analysis(simulation_name, steps=config.analysis_steps)

    # -----------------------------------------------------------------

    def clear_simulation_input(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing input of simulation '" + simulation_name + "' ...")

        # Local
        self.clear_simulation_input_local(simulation_name)

        # Remote
        self.clear_simulation_input_remote(simulation_name)

    # -----------------------------------------------------------------

    def clear_simulation_input_local(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing local input of simulation '" + simulation_name + "' ...")

        # Get simulation and check
        simulation = self.get_simulation(simulation_name)
        if not simulation.has_input:
            log.warning("The simulation '" + simulation_name + "' has no input")
            return

        # Local, single output directory
        if simulation.has_input_directory: fs.clear_directory(simulation.input_path)

        # Local, input file paths
        else:
            paths = self.get_input(simulation_name).paths
            for name in paths:
                path = paths[name]
                fs.remove_file(path)

    # -----------------------------------------------------------------

    def clear_simulation_input_remote(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing remote input of simulation '" + simulation_name + "' ...")

        # Get simulation
        simulation = self.get_simulation(simulation_name)
        if not simulation.has_input:
            log.warning("The simulation '" + simulation_name + "' has no input")
            return

        # Get the remote
        remote = self.get_remote(simulation.host_id)

        # Clear the input directory
        remote.clear_directory(simulation.remote_input_path)

    # -----------------------------------------------------------------

    def unfinish_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Unfinish
        self.unfinish_simulation(simulation_name)

    # -----------------------------------------------------------------

    def unfinish_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Unfinishing simulation '" + simulation_name + "' ...")

        # Clear remote output
        self.clear_simulation_output_remote(simulation_name)

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set not finished
        simulation.set_finished(False)

        # Save the simulation
        simulation.save()

    # -----------------------------------------------------------------

    def unretrieve_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Unretrieve
        self.unretrieve_simulation(simulation_name)

    # -----------------------------------------------------------------

    def unretrieve_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Unretrieving simulation '" + simulation_name + "' ...")

        # Clear local output
        self.clear_simulation_output_local(simulation_name)

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set not retrieved
        simulation.set_retrieved(False)

        # Save the simulation
        simulation.save()

    # -----------------------------------------------------------------

    def clear_simulation_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing output of simulation '" + simulation_name + "' ...")

        # Local
        self.clear_simulation_output_local(simulation_name)

        # Remote
        self.clear_simulation_output_remote(simulation_name)

    # -----------------------------------------------------------------

    def clear_simulation_output_local(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing local output of simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Remove
        output.remove_all()

    # -----------------------------------------------------------------

    def clear_simulation_output_remote(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing remote output of simulation '" + simulation_name + "' ...")

        # Get the remote output path
        output_path = self.get_remote_output_path(simulation_name)

        # Remove all output files
        self.get_remote_for_simulation(simulation_name).remove_files_in_path(output_path, startswith=self.get_simulation_prefix(simulation_name))

    # -----------------------------------------------------------------

    @lazyproperty
    def unanalyse_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("steps", "string_list", "analysis steps to revert", choices=clear_analysis_steps)

        # Return
        return definition

    # -----------------------------------------------------------------

    def unanalyse_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Parse
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.unanalyse_simulation_definition, **kwargs)

        # Unanalyse
        self.unanalyse_simulation(simulation_name, steps=config.steps)

    # -----------------------------------------------------------------

    def unanalyse_simulation(self, simulation_name, steps=None):

        """
        This function ...
        :param simulation_name:
        :param steps:
        :return:
        """

        # Debugging
        log.debug("Unanalysing simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set flags
        extraction = simulation.analysed_any_extraction and (steps is None or "extraction" in steps)
        plotting = simulation.analysed_any_plotting and (steps is None or "plotting" in steps)
        misc = simulation.analysed_any_misc and (steps is None or "misc" in steps)

        # Extraction
        if extraction: self.unanalyse_extraction(simulation_name)

        # Plotting
        if plotting: self.unanalyse_plotting(simulation_name)

        # Misc
        if misc: self.unanalyse_misc(simulation_name)

    # -----------------------------------------------------------------

    def unanalyse_extraction(self, simulation_name):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Unanalysing extraction of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear the extraction
        self.clear_extraction(simulation_name)

        # Unset flag
        simulation.unset_analysed_extraction()

        # Save simulation
        simulation.save()

    # -----------------------------------------------------------------

    def unanalyse_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Unanalysing plotting of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear the plotting
        self.clear_plotting(simulation_name)

        # Unset flag
        simulation.unset_analysed_plotting()

        # Save simulation
        simulation.save()

    # -----------------------------------------------------------------

    def unanalyse_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Unanalysing misc analysis of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear the misc analysis
        self.clear_misc(simulation_name)

        # Unset flag
        simulation.unset_analysed_misc()

        # Save simulation
        simulation.save()

    # -----------------------------------------------------------------

    def clear_simulation_analysis(self, simulation_name, steps=None):

        """
        This function ...
        :param simulation_name:
        :param steps:
        :return:
        """

        # Debugging
        log.debug("Clearing analysis of simulation '" + simulation_name + "' ...")

        # Extraction
        if steps is None or "extraction" in steps: self.clear_extraction(simulation_name)

        # Plotting
        if steps is None or "plotting" in steps: self.clear_plotting(simulation_name)

        # Misc
        if steps is None or "misc" in steps: self.clear_misc(simulation_name)

    # -----------------------------------------------------------------

    def clear_extraction_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _clear_command_name + " " + _extraction_command_name

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Clear
        self.clear_extraction(simulation_name)

    # -----------------------------------------------------------------

    def clear_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing extraction analysis output for simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear
        if simulation.has_extraction_output: fs.clear_directory(simulation.analysis.extraction.path)

    # -----------------------------------------------------------------

    def clear_plotting_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _clear_command_name + " " + _plotting_command_name

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Clear
        self.clear_plotting(simulation_name)

    # -----------------------------------------------------------------

    def clear_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing plotting analysis output for simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear
        if simulation.has_plotting_output: fs.clear_directory(simulation.analysis.plotting.path)

    # -----------------------------------------------------------------

    def clear_misc_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _clear_command_name + " " + _misc_command_name

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, name=name, **kwargs)

        # Clear
        self.clear_misc(simulation_name)

    # -----------------------------------------------------------------

    def clear_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Clearing misc analysis output for simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Clear
        if simulation.has_misc_output: fs.clear_directory(simulation.analysis.misc.path)

    # -----------------------------------------------------------------

    def _cache_directories(self, output, directory_path, cache_path, output_types=None):

        """
        This function ...
        :param output:
        :param directory_path:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching directories ...")

        # Get the directory paths
        if output_types is not None: dirpaths = output.get_all_directory_paths_not_in_directory_for_types(output_types)
        else: dirpaths = output.all_directory_paths_not_in_directory

        # Loop over the directories
        for dirpath in dirpaths:

            # Determine relative directory path
            relpath = fs.relative_to(dirpath, directory_path)

            # Determine new path
            new_path = fs.join(cache_path, relpath)

            # Get directory name and containing directory
            dirname = fs.name(dirpath)
            cache_directory_path = fs.directory_of(dirpath)
            if not fs.is_directory(cache_directory_path): fs.create_directory(cache_directory_path, recursive=True)

            # Debugging
            log.debug("Caching '" + dirname + "' directory from '" + fs.directory_of(relpath) + "' to '" + cache_directory_path + "' ...")

            # Copy the directory
            fs.copy_directory(dirpath, cache_directory_path)

            # Check whether the directory is present
            if not fs.is_directory(new_path): raise RuntimeError("Something went wrong")

            # Remove the original directory
            fs.remove_directory(dirpath)

    # -----------------------------------------------------------------

    def _cache_files(self, output, directory_path, cache_path, output_types=None):

        """
        This function ...
        :param output:
        :param directory_path:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching files ...")

        # Get the file paths
        if output_types is not None: filepaths = output.get_all_file_paths_not_in_directory_for_types(output_types)
        else: filepaths = output.all_file_paths_not_in_directory

        # Loop over the files
        for filepath in filepaths:

            # Determine relative file path
            relpath = fs.relative_to(filepath, directory_path)

            # Determine the new path
            new_path = fs.join(cache_path, relpath)

            # Get filename and containing directory
            filename = fs.name(filepath)
            cache_directory_path = fs.directory_of(filepath)
            if not fs.is_directory(cache_directory_path): fs.create_directory(cache_directory_path, recursive=True)

            # Debugging
            log.debug("Caching '" + filename + "' file from '" + fs.directory_of(relpath) + "' to '" + cache_directory_path + "' ...")

            # Copy the file
            fs.copy_file(filepath, cache_directory_path)

            # Check whether the file is present
            if not fs.is_file(new_path): raise RuntimeError("Something went wrong")

            # Remove the original file
            fs.remove_file(filepath)

    # -----------------------------------------------------------------

    def _cache_output(self, output, directory_path, cache_path, output_types=None):

        """
        This function ...
        :param output:
        :param directory_path: local output directory
        :param cache_path: cache directory path
        :param output_types:
        :return:
        """

        # Set cache path (write cached.dat file)
        write_cache_path(directory_path, cache_path)

        # Directories
        self._cache_directories(output, directory_path, cache_path, output_types=output_types)

        # Files
        self._cache_files(output, directory_path, cache_path, output_types=output_types)

        # Remove original output files and directories: NO, depends on output_types
        #output.remove_all()

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_output_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("output_types", "string_list", "output types to cache", choices=output_type_choices)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_output_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _output_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_output_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether there is output
            if not self.has_output(simulation_name):
                log.warning("No output for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_output(simulation_name):
                log.warning("Output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_output(simulation_name, cache_path=config.path, output_types=config.output_types)

    # -----------------------------------------------------------------

    @memoize_method
    def get_cache_output_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Root path is given
        if self.config.cache_root is not None:

            # Get the relative simulation output path
            simulation_output_path = self.get_output_path(simulation_name)
            relative_output_path = fs.relative_to(simulation_output_path, self.config.cache_root)

            # Determine the output cache path
            output_cache_path = fs.create_directory_in(self.config.cache_path, relative_output_path, recursive=True)

        # Root path is not given: create directory for the simulation in the given cache path
        else:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(self.config.cache_path, simulation_name)

            # Create output cache path
            output_cache_path = fs.create_directory_in(simulation_cache_path, "out")

        # Return the path
        return output_cache_path

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_datacubes_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_datacubes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _datacube_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_datacubes_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether the simulation has datacubes
            if not self.has_datacubes(simulation_name):
                log.warning("No datacubes for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_output(simulation_name):
                log.warning("Output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_datacubes(simulation_name, cache_path=config.path)

    # -----------------------------------------------------------------

    def cache_simulation_datacubes(self, simulation_name, cache_path=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :return:
        """

        # Cache
        self.cache_simulation_output(simulation_name, cache_path=cache_path, output_types=datacube_types)

    # -----------------------------------------------------------------

    def cache_simulation_output(self, simulation_name, cache_path=None, output_types=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching output of simulation '" + simulation_name + "' ...")

        # Get the simulation output path
        output_path = self.get_output_path(simulation_name)

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Create the cache output path
        if cache_path is not None:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(cache_path, simulation_name)

            # Create output cache path
            cache_output_path = fs.create_directory_in(simulation_cache_path, "out")

        # Get the cache output path
        else: cache_output_path = self.get_cache_output_path_for_simulation(simulation_name)

        # Cache
        self._cache_output(output, output_path, cache_output_path, output_types=output_types)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_analysis_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _analysis_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_analysis_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether there is analysis output
            if not self.has_analysis_output(simulation_name):
                log.warning("No analysis output for simulation '" + simulation_name + "': skipping ...")
                continue

            # Cache
            self.cache_simulation_analysis(simulation_name, cache_path=config.path)

    # -----------------------------------------------------------------

    def cache_simulation_analysis(self, simulation_name, cache_path=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :return:
        """

        # Debugging
        log.debug("Caching analysis output of simulation '" + simulation_name + "' ...")

        # Extraction
        if not self.is_cached_extraction(simulation_name): self.cache_simulation_extraction(simulation_name, cache_path=cache_path)
        else: log.warning("Extraction output of simulation '" + simulation_name + "' is already cached: skipping ...")

        # Plotting
        if not self.is_cached_plotting(simulation_name): self.cache_simulation_plotting(simulation_name, cache_path=cache_path)
        else: log.warning("Plotting output of simulation '" + simulation_name + "' is already cached: skipping ...")

        # Misc
        if not self.is_cached_misc(simulation_name): self.cache_simulation_misc(simulation_name, cache_path=cache_path)
        else: log.warning("Misc output of simulation '" + simulation_name + "' is already cached: skipping ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_extraction_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("output_types", "string_list", "extraction output types to cache", choices=extraction_output_type_choices)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_extraction_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _extraction_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_extraction_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether there is extraction output
            if not self.has_extraction_output(simulation_name):
                log.warning("No extraction output for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_extraction(simulation_name):
                log.warning("Extraction output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_extraction(simulation_name, cache_path=config.path, output_types=config.output_types)

    # -----------------------------------------------------------------

    @memoize_method
    def get_cache_extraction_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Root path is given
        if self.config.cache_root is not None:

            # Get the relative simulation extraction path
            simulation_extraction_path = self.get_extraction_path(simulation_name)
            relative_extraction_path = fs.relative_to(simulation_extraction_path, self.config.cache_root)

            # Determine the extraction cache path
            extraction_cache_path = fs.create_directory_in(self.config.cache_path, relative_extraction_path, recursive=True)

        # Root path is not given: create directory for the simulation in the given cache path
        else:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(self.config.cache_path, simulation_name)

            # Create extraction cache path
            extraction_cache_path = fs.create_directory_in(simulation_cache_path, "extr")

        # Return the path
        return extraction_cache_path

    # -----------------------------------------------------------------

    def cache_simulation_extraction(self, simulation_name, cache_path=None, output_types=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching extraction output of simulation '" + simulation_name + "' ...")

        # Get the simulation extraction path
        extraction_path = self.get_extraction_path(simulation_name)

        # Get the simulation extraction output
        extraction = self.get_extraction_output(simulation_name)

        # Create the cache extraction path
        if cache_path is not None:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(cache_path, simulation_name)

            # Create extraction cache path
            cache_extraction_path = fs.create_directory_in(simulation_cache_path, "extr")

        # Get the cache extraction path
        else: cache_extraction_path = self.get_cache_extraction_path_for_simulation(simulation_name)

        # Cache
        self._cache_output(extraction, extraction_path, cache_extraction_path, output_types=output_types)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_plotting_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("output_types", "string_list", "plotting output types to cache", choices=plotting_output_type_choices)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_plotting_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _plotting_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_plotting_definition, name=name)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether there is plotting output
            if not self.has_plotting_output(simulation_name):
                log.warning("No plotting output for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_plotting(simulation_name):
                log.warning("Plotting output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_plotting(simulation_name, cache_path=config.path, output_types=config.output_types)

    # -----------------------------------------------------------------

    @memoize_method
    def get_cache_plotting_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Root path is given
        if self.config.cache_root is not None:

            # Get the relative simulation plotting path
            simulation_plotting_path = self.get_plotting_path(simulation_name)
            relative_plotting_path = fs.relative_to(simulation_plotting_path, self.config.cache_root)

            # Determine the plotting cache path
            plotting_cache_path = fs.create_directory_in(self.config.cache_path, relative_plotting_path, recursive=True)

        # Root path is not given: create directory for the simulation in the given cache path
        else:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(self.config.cache_path, simulation_name)

            # Create plotting cache path
            plotting_cache_path = fs.create_directory_in(simulation_cache_path, "plot")

        # Return the path
        return plotting_cache_path

    # -----------------------------------------------------------------

    def cache_simulation_plotting(self, simulation_name, cache_path=None, output_types=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching plotting output of simulation '" + simulation_name + "' ...")

        # Get the simulation plotting path
        plotting_path = self.get_plotting_path(simulation_name)

        # Get the simulation output
        plotting = self.get_output(simulation_name)

        # Create the cache plotting path
        if cache_path is not None:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(cache_path, simulation_name)

            # Create plotting cache path
            cache_plotting_path = fs.create_directory_in(simulation_cache_path, "plot")

        # Get the cache plotting path
        else: cache_plotting_path = self.get_cache_plotting_path_for_simulation(simulation_name)

        # Cache
        self._cache_output(plotting, plotting_path, cache_plotting_path, output_types=output_types)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_misc_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("output_types", "string_list", "misc output types to cache", choices=misc_output_type_choices)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_misc_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _misc_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_misc_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether there is misc output
            if not self.has_misc_output(simulation_name):
                log.warning("No miscellaneous output for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_misc(simulation_name):
                log.warning("Miscellaneous output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_misc(simulation_name, cache_path=config.path, output_types=config.output_types)

    # -----------------------------------------------------------------

    @memoize_method
    def get_cache_misc_path_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Root path is given
        if self.config.cache_root is not None:

            # Get the relative simulation misc path
            simulation_misc_path = self.get_misc_path(simulation_name)
            relative_misc_path = fs.relative_to(simulation_misc_path, self.config.cache_root)

            # Determine the misc cache path
            misc_cache_path = fs.create_directory_in(self.config.cache_path, relative_misc_path, recursive=True)

        # Root path is not given: create directory for the simulation in the given cache path
        else:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(self.config.cache_path, simulation_name)

            # Create misc cache path
            misc_cache_path = fs.create_directory_in(simulation_cache_path, "misc")

        # Return the path
        return misc_cache_path

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_simulations_images_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("path", "directory_path", "caching path")
        return definition

    # -----------------------------------------------------------------

    def cache_simulations_images_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _cache_command_name + " " + _images_command_name

        # Parse
        splitted, simulation_names, config = self.parse_simulations_command(command, self.cache_simulations_images_definition, name=name, **kwargs)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check whether the simulation has images
            if not self.has_images(simulation_name):
                log.warning("No images for simulation '" + simulation_name + "': skipping ...")
                continue

            # Check whether already cached
            if self.is_cached_misc(simulation_name):
                log.warning("Miscellaneous output for simulation '" + simulation_name + "' is already cached: skipping ...")
                continue

            # Cache
            self.cache_simulation_images(simulation_name, cache_path=config.path)

    # -----------------------------------------------------------------

    def cache_simulation_images(self, simulation_name, cache_path=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :return:
        """

        # Cache
        self.cache_simulation_misc(simulation_name, cache_path=cache_path, output_types=image_types)

    # -----------------------------------------------------------------

    def cache_simulation_misc(self, simulation_name, cache_path=None, output_types=None):

        """
        This function ...
        :param simulation_name:
        :param cache_path:
        :param output_types:
        :return:
        """

        # Debugging
        log.debug("Caching misc output of simulation '" + simulation_name + "' ...")

        # Get the simulation misc path
        misc_path = self.get_misc_path(simulation_name)

        # Get the simulation output
        misc = self.get_misc_output(simulation_name)

        # Create the cache output path
        if cache_path is not None:

            # Create simulation cache path
            simulation_cache_path = fs.create_directory_in(cache_path, simulation_name)

            # Create misc cache path
            cache_misc_path = fs.create_directory_in(simulation_cache_path, "misc")

        # Get the cache output path
        else: cache_misc_path = self.get_cache_misc_path_for_simulation(simulation_name)

        # Cache
        self._cache_output(misc, misc_path, cache_misc_path, output_types=output_types)

    # -----------------------------------------------------------------

    @lazyproperty
    def relaunch_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("finished", "relaunch already finished simulations", False)
        definition.add_positional_optional("parallelization", "parallelization", "new parallelization scheme (by default, previous one is used)")
        definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)

        # Return
        return definition

    # -----------------------------------------------------------------

    def relaunch_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.relaunch_simulation_definition, **kwargs)

        # Check simulation
        if not self.is_failed(simulation_name):
            if config.finished and self.is_finished(simulation_name):
                log.info("Simulation '" + simulation_name + "' is already finished, removing local and remote output and undoing analysis ...")
                if self.is_retrieved(simulation_name): self.unretrieve_simulation(simulation_name)
                self.clear_simulation_output_remote(simulation_name)
                self.unanalyse_simulation(simulation_name)
            else: raise ValueError("Simulation '" + simulation_name + "' is running, finished or still queued")

        # Relaunch the simulation
        self.relaunch_simulation(simulation_name, parallelization=config.parallelization, logging_options=config.logging, scheduling_options=config.scheduling)

    # -----------------------------------------------------------------

    def relaunch_simulation(self, simulation_name, parallelization=None, logging_options=None, scheduling_options=None):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :return:
        """

        # Debugging
        log.debug("Relaunching simulation '" + simulation_name + "' ...")

        # Screen
        if self.is_screen_execution(simulation_name): self.relaunch_simulation_screen(simulation_name, parallelization=parallelization, logging_options=logging_options)

        # Job
        elif self.is_job_execution(simulation_name): self.relaunch_simulation_job(simulation_name, parallelization=parallelization, logging_options=logging_options, scheduling_options=scheduling_options)

        # Not supported
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    def relaunch_simulation_screen(self, simulation_name, parallelization=None, logging_options=None):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :param logging_options:
        :return:
        """

        # Get the screen session name
        screen_name = self.get_screen_name(simulation_name)

        # Debugging
        log.debug("Relaunching the screen session '" + screen_name + "' to relaunch simulation '" + simulation_name + "' ...")

        # Get the remote
        remote = self.get_remote_for_simulation(simulation_name)

        # Check if the screen is active: this will only be the case for the first time the relaunch_simulation_screen function is called (for the same screen)
        if remote.is_active_screen(screen_name):

            # Give warning
            log.warning("The screen '" + screen_name + "' is active: currently running simulations will have to be stopped and relaunched for simulation '" + simulation_name + "'")

            # Stop the screen
            remote.kill_screen(screen_name)

            # Check which simulations are still queued or running in this screen

            # Stop/cancel and re-add the simulations to the queue (they will be launched together with the relaunched simulations in a new screen)

            # Remove their remote output
            self.unfinish_simulation(simulation_name)

        # Set the parallelization scheme
        if parallelization is None:
            screen_script = self.get_screen_script(simulation_name)
            if screen_script is None: raise RuntimeError("Cannot determine parallelization scheme: original screen script not found")
            parallelization = screen_script.get_parallelization(simulation_name)

        # Set the logging options
        if logging_options is None:
            screen_script = self.get_screen_script(simulation_name)
            if screen_script is None: raise RuntimeError("Cannot determine logging options: original screen script not found")
            logging_options = screen_script.get_logging_options(simulation_name)

        # Get the simulation
        simulation = self.get_simulation(simulation_name)
        host_id = self.get_host_id_for_simulation(simulation_name)

        # Unanalyse: remove analysis output and unset analysis flags
        self.unanalyse_simulation(simulation_name)

        # Unretrieve: remove local simulation output and unset retrieved flag
        self.unretrieve_simulation(simulation_name)

        # Unfinish: remove remote simulation output and unset finished flag
        self.unfinish_simulation(simulation_name)

        # Add the simulation to the queue
        self.launcher.add_to_queue(simulation.definition, simulation_name, host_id=host_id, parallelization=parallelization,
                                   logging_options=logging_options, analysis_options=simulation.analysis)

        # Add simulation to the relaunched table
        self.relaunched.add_simulation(simulation)

        # Remove original simulation?
        if not self.config.dry: self.remove_simulation(simulation_name)  # backup is done here

        # Don't remove
        else: log.warning("[DRY] Not removing the simulation from '" + simulation.path + "' ...")

    # -----------------------------------------------------------------

    def relaunch_simulation_job(self, simulation_name, parallelization=None, logging_options=None, scheduling_options=None):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :return:
        """

        # Debugging
        log.debug("Relaunching simulation '" + simulation_name + "' as a new job ...")

        # Set the parallelization scheme
        if parallelization is None:
            job_script = self.get_job_script(simulation_name)
            if job_script is None: raise RuntimeError("Cannot determine parallelization scheme: original job script not found")
            parallelization = job_script.parallelization

        # Set the logging options
        if logging_options is None:
            job_script = self.get_job_script(simulation_name)
            if job_script is None: raise RuntimeError("Cannot determine logging options: original job script not found")
            logging_options = job_script.logging_options

        # Set the scheduling options
        if scheduling_options is None:
            job_script = self.get_job_script(simulation_name)
            if job_script is None: raise RuntimeError("Cannot determine scheduling options: original job script not found")
            scheduling_options = job_script.scheduling_options

        # Get the simulation
        simulation = self.get_simulation(simulation_name)
        host_id = self.get_host_id_for_simulation(simulation_name)

        # TODO: actually support giving the specific cluster it is supposed to be run on

        # Add the simulation to the queue
        self.launcher.add_to_queue(simulation.definition, simulation_name, host_id=host_id,
                                   parallelization=parallelization, logging_options=logging_options,
                                   scheduling_options=scheduling_options, analysis_options=simulation.analysis)

        # Add simulation to the relaunched table
        self.relaunched.add_simulation(simulation)

        # Remove original simulation?
        if not self.config.dry: self.remove_simulation(simulation_name)  # backup is done here

        # Don't remove
        else: log.warning("[DRY] Not removing the simulation from '" + simulation.path + "' ...")

    # -----------------------------------------------------------------

    @property
    def nrelaunched(self):

        """
        This function ...
        :return:
        """

        return len(self.relaunched)

    # -----------------------------------------------------------------

    @property
    def has_relaunched(self):

        """
        This function ...
        :return:
        """

        return self.nrelaunched > 0

    # -----------------------------------------------------------------

    @property
    def relaunched_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.relaunched.simulation_names

    # -----------------------------------------------------------------

    def is_relaunched(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.relaunched_simulation_names

    # -----------------------------------------------------------------

    def set_simulation_options_from_original(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Get the original simulation
        original_simulation = self.get_simulation(simulation.name)

        # Analyser paths
        simulation.analyser_paths = original_simulation.analyser_paths

        # Options for retrieval
        simulation.retrieve_types = original_simulation.retrieve_types

        # Options for removing remote or local input and output
        simulation.remove_remote_input = original_simulation.remove_remote_input
        simulation.remove_remote_output = original_simulation.remove_remote_output  # After retrieval
        simulation.remove_remote_simulation_directory = original_simulation.remove_remote_simulation_directory  # After retrieval
        simulation.remove_local_output = original_simulation.remove_local_output  # After analysis

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching simulations ...")

        # Get remotes
        remotes = []
        for host_id in self.launcher.queue_host_ids: remotes.append(self.get_remote(host_id))

        # Debugging
        log.debug("Running the batch simulation launcher ...")

        # Run the launcher
        self.launcher.run(remotes=remotes, timing=self.timing, memory=self.memory)

        # Set moved IDs
        self.set_moved_simulation_ids()

        # Set relaunched simulation IDs
        self.set_relaunched_simulation_ids()

        # Debugging
        log.debug("Creating backup of assignment scheme ...")

        # Backup assignment?
        if self.config.backup_assignment: self.assignment.saveto(self.backup_assignment_path)

        # Debugging
        log.debug("Updating the assignment scheme ...")

        # Change assignment
        for simulation in self.launched_simulations: self.assignment.update_simulation(simulation)

        # Set that the assignment has changed
        self._adapted_assignment = True

        # Reset status?
        if self.has_moved or self.has_relaunched: self.reset_status()

    # -----------------------------------------------------------------

    def set_moved_simulation_ids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the new simulation IDs for the moved simulations ...")

        # Set moved IDs
        for simulation in self.moved_simulations:

            # Set options
            self.set_simulation_options_from_original(simulation)

            # Save the new simulation
            if not self.config.dry: simulation.save()

            # Set the new ID of the moved simulation
            self.moved.set_new_id(simulation.name, simulation.id)

    # -----------------------------------------------------------------

    def set_relaunched_simulation_ids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the new simulation IDs for the relaunched simulations ...")

        # Set relaunched IDs
        for simulation in self.relaunched_simulations:

            # Set options
            self.set_simulation_options_from_original(simulation)

            # Save the new simulation
            if not self.config.dry: simulation.save()

            # Set the new ID of the relaunched simulation
            self.relaunched.set_new_id(simulation.name, simulation.id)

    # -----------------------------------------------------------------

    @property
    def launched_simulations(self):

        """
        This function ...
        :return:
        """

        return self.launcher.launched_simulations

    # -----------------------------------------------------------------

    @property
    def nlaunched_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.launched_simulations)

    # -----------------------------------------------------------------

    @lazyproperty
    def launched_simulation_names(self):

        """
        This function ...
        :return:
        """

        return [simulation.name for simulation in self.launched_simulations]

    # -----------------------------------------------------------------

    @lazyproperty
    def moved_simulations(self):

        """
        This function ...
        :return:
        """

        simulations = []
        for simulation in self.launched_simulations:
            if not self.is_moved(simulation.name): continue
            simulations.append(simulation)
        return simulations

    # -----------------------------------------------------------------

    @property
    def nmoved_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.moved_simulations)

    # -----------------------------------------------------------------

    @lazyproperty
    def relaunched_simulations(self):

        """
        This function ...
        :return:
        """

        simulations = []
        for simulation in self.launched_simulations:
            if not self.is_relaunched(simulation.name): continue
            simulations.append(simulation)
        return simulations

    # -----------------------------------------------------------------

    @property
    def nrelaunched_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.relaunched_simulations)

    # -----------------------------------------------------------------

    @memoize_method
    def has_timing_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.has_timing and self.timing.has_simulation(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def has_memory_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.has_memory and self.memory.has_simulation(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def timing_index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.timing.index_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def memory_index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.memory.index_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            total_time = self.timing.total_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add
            times[host][parallelization].append(total_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            setup_time = self.timing.setup_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(setup_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            stellar_time = self.timing.stellar_emission_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(stellar_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictinoary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            spectra_time = self.timing.spectra_calculation_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(spectra_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            dust_time = self.timing.dust_emission_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(dust_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            writing_time = self.timing.writing_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(writing_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            waiting_time = self.timing.waiting_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(waiting_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            communication_time = self.timing.communication_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(communication_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def intermediate_times(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        times = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_timing_for_simulation(simulation_name): continue

            # Get timing
            intermediate_time = self.timing.intermediate_time_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the runtime
            times[host][parallelization].append(intermediate_time)

        # Return the runtimes
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def total_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            total_memory = self.memory.total_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the memory usage
            memories[host][parallelization].append(total_memory)

        # Return the memory usages
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            setup_memory = self.memory.setup_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the memory usage
            memories[host][parallelization].append(setup_memory)

        # Return the memory usages
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            stellar_memory = self.memory.stellar_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add memory usage
            memories[host][parallelization].append(stellar_memory)

        # Return the memory usages
        return memories

    #  -----------------------------------------------------------------

    @lazyproperty
    def spectra_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            spectra_memory = self.memory.spectra_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the memory usage
            memories[host][parallelization].append(spectra_memory)

        # Return
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            dust_memory = self.memory.dust_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the memory usage
            memories[host][parallelization].append(dust_memory)

        # Return
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_memories(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the different hosts
        memories = create_nested_defaultdict(2, list)

        # Loop over the simulations
        for simulation_name in self.simulation_names:

            # Check
            if not self.has_memory_for_simulation(simulation_name): continue

            # Get memory usage
            writing_memory = self.memory.writing_memory_for_simulation(simulation_name)

            # Get the host for the simulation
            host = self.get_host_for_simulation(simulation_name)

            # Get the parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)

            # Add the memory usage
            memories[host][parallelization].append(writing_memory)

        # Return
        return memories

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Assignment scheme
        if self.config.show_assignment: self.show_assignment()

        # Status
        if self.config.show_status: self.show_status()

        # Runtimes
        if self.config.show_runtimes: self.show_runtimes()

        # Memory
        if self.config.show_memory: self.show_memory()

    # -----------------------------------------------------------------

    def get_scalar_times(self, times_dict):

        """
        This function ...
        :param times_dict:
        :return:
        """

        times = create_nested_defaultdict(2, list)
        for host in times_dict:
            for parallelization in times_dict[host]:
                for runtime in times_dict[host][parallelization]: times[host][parallelization].append(runtime.to("min").value)
        return times

    # -----------------------------------------------------------------

    @lazyproperty
    def total_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.total_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.setup_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.stellar_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.spectra_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.dust_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.writing_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.waiting_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.communication_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def intermediate_times_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_times(self.intermediate_times)

    # -----------------------------------------------------------------

    def get_scalar_memories(self, memories_dict):

        """
        This function ...
        :param memories_dict:
        :return:
        """

        memories = create_nested_defaultdict(2, list)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                for memory in memories_dict[host][parallelization]: memories[host][parallelization].append(memory.to("GB").value)
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def total_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.total_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.setup_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.stellar_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.spectra_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.dust_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_memories_scalar(self):

        """
        This function ...
        :return:
        """

        return self.get_scalar_memories(self.writing_memories)

    # -----------------------------------------------------------------

    def get_ntimes(self, times_dict):

        """
        This function ...
        :param times_dict:
        :return:
        """

        ntimes = defaultdict(dict)
        for host in times_dict:
            for parallelization in times_dict[host]:
                ntimes[host][parallelization] = len(times_dict[host][parallelization])
        return ntimes

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.total_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.setup_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.stellar_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.spectra_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.dust_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.writing_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwaiting_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.waiting_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def ncommunication_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.communication_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def nintermediate_times(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.intermediate_times)

    # -----------------------------------------------------------------

    def get_nmemories(self, memories_dict):

        """
        This function ...
        :param memories_dict:
        :return:
        """

        nmemories = defaultdict(dict)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                nmemories[host][parallelization] = len(memories_dict[host][parallelization])
        return nmemories

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.total_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.setup_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.stellar_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.spectra_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.dust_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.writing_memories)

    # -----------------------------------------------------------------

    def get_times_clipped(self, times_dict):

        """
        This function ...
        :param times_dict:
        :return:
        """

        clipped = defaultdict(dict)
        for host in times_dict:
            for parallelization in times_dict[host]:
                clipped[host][parallelization] = numbers.sigma_clip(times_dict[host][parallelization])
        return clipped

    # -----------------------------------------------------------------

    @lazyproperty
    def total_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.total_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.setup_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.stellar_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.spectra_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.dust_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.writing_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.waiting_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.communication_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def intermediate_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_times_clipped(self.intermediate_times_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.total_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.ntotal_times, self.ntotal_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.setup_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nsetup_times, self.nsetup_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.stellar_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nstellar_times, self.nstellar_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.spectra_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nspectra_times, self.nspectra_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.dust_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.ndust_times, self.ndust_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.writing_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nwriting_times, self.nwriting_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwaiting_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.waiting_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwaiting_times_outliers(self):

        """
        Thisf unction ...
        :return:
        """

        return get_noutliers(self.nwaiting_times, self.nwaiting_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ncommunication_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.communication_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ncommunication_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.ncommunication_times, self.ncommunication_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nintermediate_times_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_ntimes(self.intermediate_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nintermediate_times_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nintermediate_times, self.nintermediate_times_clipped)

    # -----------------------------------------------------------------

    def get_memories_clipped(self, memories_dict):

        """
        This function ...
        :param memories_dict:
        :return:
        """

        clipped = defaultdict(dict)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                clipped[host][parallelization] = numbers.sigma_clip(memories_dict[host][parallelization])
        return clipped

    # -----------------------------------------------------------------

    @lazyproperty
    def total_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.total_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.setup_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.stellar_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.spectra_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.dust_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_clipped(self.writing_memories_scalar)

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.total_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ntotal_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.ntotal_memories, self.ntotal_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.setup_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nsetup_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nsetup_memories, self.nsetup_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.stellar_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nstellar_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nstellar_memories, self.nstellar_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.spectra_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectra_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nspectra_memories, self.nspectra_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.dust_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.ndust_memories, self.ndust_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_memories_clipped(self):

        """
        This function ...
        :return:
        """

        return self.get_nmemories(self.writing_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwriting_memories_outliers(self):

        """
        This function ...
        :return:
        """

        return get_noutliers(self.nwriting_memories, self.nwriting_memories_clipped)

    # -----------------------------------------------------------------

    def show_assignment(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the assignment scheme ...")

        # Show
        print(self.assignment)

    # -----------------------------------------------------------------

    def get_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Add host ID
        if host_id not in self.remotes: self.remotes.add_host_id(host_id)

        # Return the remote
        return self.remotes[host_id]

    # -----------------------------------------------------------------

    def get_home_path(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_remote(host_id).home_directory

    # -----------------------------------------------------------------

    def set_remote(self, remote, name=None):

        """
        This function ...
        :param remote:
        :param name:
        :return:
        """

        # Add
        self.remotes.set_or_add_remote(remote, name=name)

    # -----------------------------------------------------------------

    def get_display_name(self, simulation, add_quotes=False, id_size=None, host_id_size=None):

        """
        This function ...
        :param simulation:
        :param add_quotes:
        :param id_size:
        :param host_id_size:
        :return:
        """

        # Set simulation ID string
        if id_size is not None: id_string = strings.integer(simulation.id, id_size)
        else: id_string = simulation.id

        # Set host ID string
        if host_id_size is not None: host_string = strings.to_length(tostr(simulation.host), host_id_size)
        else: host_string = tostr(simulation.host)

        # Create name and return
        if add_quotes: return "'" + simulation.name + "' (" + host_string + " " + id_string + ")"
        else: return simulation.name + " (" + host_string + " " + id_string + ")"

    # -----------------------------------------------------------------

    @lazyproperty
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def status_column_names(self):

        """
        This function ...
        :return:
        """

        names = ["", "", "Name", "Host", "ID"]
        if self.has_info: names.extend(self.info_names)
        names.append("Status")
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def status_column_units(self):

        """
        This function ...
        :return:
        """

        units = ["", "", "", "", ""]
        if self.has_info: units.extend(self.info_unit_strings)
        units.append("")
        return units

    # -----------------------------------------------------------------

    @lazyproperty
    def show_status_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("extra_columns", "string_list", "extra columns to show", choices=extra_columns)
        definition.add_optional("path", "string", "save the status information as a table at this path")
        definition.add_flag("refresh", "refresh the status info", False)

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
        self.show_status(extra=config.extra_columns, path=config.path, refresh=config.refresh)

    # -----------------------------------------------------------------

    def show_status(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing the simulation status ...")

        # Get settings
        extra = kwargs.pop("extra", None)
        path = kwargs.pop("path", None)
        refresh = kwargs.pop("refresh", False)

        # Refresh if requested
        if refresh: self.reset_status()

        # Show
        print("")
        print(fmt.bold + "Total number of simulations: " + fmt.reset + str(self.nsimulations))
        print(fmt.bold + "Number of finished simulations: " + fmt.reset + str(self.nfinished) + " (" + tostr(self.percentage_nfinished, round=True, ndigits=2) + "%)")
        print(fmt.bold + "Number of retrieved simulations: " + fmt.reset + str(self.nretrieved) + " (" + tostr(self.percentage_nretrieved, round=True, ndigits=2) + "%)")
        print(fmt.bold + "Number of analysed simulations: " + fmt.reset + str(self.nanalysed) + " (" + tostr(self.percentage_nanalysed, round=True, ndigits=2) + "%)")
        print("")

        # If path is given, initialize table
        #if path is not None: table = SmartTable()
        #else: table = None
        table = None

        # Print in columns
        with fmt.print_in_columns() as print_row:

            # Set the column names and units
            column_names = []
            column_units = []
            for name in self.status_column_names: column_names.append(name)
            for unit in self.status_column_units: column_units.append(unit)
            if extra is not None:
               for col in extra:
                   column_names.append(extra_column_names[col])
                   column_units.append(extra_column_units[col] if col in extra_column_units else None)

            # Show the header
            print_row(*column_names)
            if self.has_info: print_row(*column_units)

            # Initialize columns
            if path is not None: columns = [[] for _ in column_names]
            else: columns = None

            # Loop over the simulations
            for index, simulation_name in enumerate(self.simulation_names):

                # Get the simulation
                simulation = self.get_simulation(simulation_name)

                # Get the status
                status = self.get_status(simulation_name)

                # Get index string
                index_string = "[" + strings.integer(index, 3, fill=" ") + "] "

                # Set color
                if status == analysed_name: color = "green"
                elif status == retrieved_name: color = "yellow"
                elif status == finished_name: color = "yellow"
                elif is_running_status(status): color = None
                else: color = "red"

                # Create strings
                host_string = tostr(simulation.host)
                id_string = tostr(simulation.id)

                # Set parts
                parts = []
                parts.append(" - ")
                parts.append(index_string)
                parts.append(simulation_name)
                parts.append(host_string)
                parts.append(id_string)

                # Add info
                if self.has_info:
                    for name in self.info_names:

                        info = self.get_info_for_simulation(simulation_name)
                        if name not in info: string = "--"
                        else:
                            value = info[name]
                            unit = self.info_units[name]
                            if unit is not None: value = value.to(unit).value
                            string = tostr(value, scientific=self.config.info_scientific, decimal_places=self.config.info_ndecimal_places)

                        parts.append(string)

                # Add the simulation status
                parts.append(status)

                # Add extra columns
                if extra is not None:
                    for col in extra:
                        if col == _screen_extra_name: value = self.get_screen_name(simulation_name)
                        elif col == _job_extra_name: value = self.get_job_id_string(simulation_name)
                        elif col == _disk_extra_name: value = self.get_disk_size(simulation_name)
                        elif col == _runtime_extra_name: value = self.get_runtime(simulation_name)
                        elif col == _memory_extra_name: value = self.get_peak_memory(simulation_name)
                        else: raise ValueError("Invalid extra column name: '" + col + "'")
                        if col in extra_column_units and value is not None: value = value.to(extra_column_units[col]).value
                        string = tostr(value) if value is not None else "--"
                        parts.append(string)

                # Print the row
                print_row(*parts, color=color)

                # Add to the columns
                if columns is not None:
                    for j, part in enumerate(parts): columns[j].append(part)

            # Debugging
            log.debug("Creating status info table ...")

            # Create table from the columns
            if columns is not None: table = SmartTable.from_columns(*columns, names=column_names, units=column_units)

        # End with empty line
        print("")

        # Save the table
        if table is not None and path is not None: table.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def find_simulations_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("host_id", "string", "host ID for which to find simulation(s)")
        return definition

    # -----------------------------------------------------------------

    def find_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Parse the command
        config = self.get_config_from_command(command, self.find_simulations_definition, **kwargs)

        # No info?
        if not self.has_info: raise NotImplementedError("Not yet implemented without info: based on any choosen set of parameters")

        # Find simulations
        self.find_simulations(host_id=config.host_id)

    # -----------------------------------------------------------------

    @memoize_method # memoize because anyway self.simulation_names is a lazyproperty
    def get_index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.simulation_names.index(simulation_name)

    # -----------------------------------------------------------------

    def find_simulations(self, host_id=None):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Inform the user
        log.info("Finding simulations based on parameters ...")

        # Get subset of simulations to start with
        if host_id is not None: simulation_names = self.get_simulation_names_for_host_id(host_id)
        else: simulation_names = None

        # Loop over the different parameters
        for name in self.info_names:

            # Get all info for all simulations (dict)
            values = self.get_info_values_scalar(name)

            # Get subset of values (still dict)
            if simulation_names is not None: values_simulations = create_subdict(values, simulation_names)
            else: values_simulations = values

            # Get unique values
            unique_values = list(sorted(sequences.unique_values(values_simulations.values(), ignore_none=True)))

            # Create description
            description = "choice for '" + name + "'"
            if self.info_units[name] is not None: description += " [" + tostr(self.info_units[name]) + "]"

            # Prompt for value
            choice = prompt_choice(name, description, unique_values)

            # Find simulation names for the choice
            simulation_names = [] # reset
            for simulation_name in values_simulations:
                value = values[simulation_name]
                if value != choice: continue
                simulation_names.append(simulation_name)

            # Show the simulation names
            nsimulations = len(simulation_names)
            if nsimulations == 0:
                log.warning("No simulations with the chosen parameter")
                break
            elif nsimulations == 1:
                log.success("One matching simulation: '" + simulation_names[0] + "'")
                break
            else:
                log.success("Matching simulations:")
                indices = []
                for simulation_name in simulation_names:
                    index = self.get_index_for_simulation(simulation_name)
                    indices.append(index)
                for index in sorted(indices):
                    index_string = "[" + strings.integer(index, 3, fill=" ") + "] "
                    simulation_name = self.simulation_names[index]
                    log.success(" - " + index_string + simulation_name)

    # -----------------------------------------------------------------

    def get_average_times(self, times_dict):

        """
        This function ...
        :param times_dict:
        :return:
        """

        averages = defaultdict(dict)
        for host in times_dict:
            for parallelization in times_dict[host]:
                averages[host][parallelization] = numbers.arithmetic_mean(*times_dict[host][parallelization])
        return averages

    # -----------------------------------------------------------------

    def get_stddev_times(self, times_dict, averages=None):

        """
        This function ...
        :param times_dict:
        :param averages:
        :return:
        """

        stddevs = defaultdict(dict)
        for host in times_dict:
            for parallelization in times_dict[host]:
                if averages is not None and host in averages and parallelization in averages[host]: average = averages[host][parallelization]
                else: average = None
                stddevs[host][parallelization] = numbers.standard_deviation(*times_dict[host][parallelization], mean=average)
        return stddevs

    # -----------------------------------------------------------------

    @lazyproperty
    def average_total_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.total_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_total_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.total_times_clipped, self.average_total_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_setup_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.setup_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_setup_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.setup_times_clipped, self.average_setup_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_stellar_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.stellar_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_stellar_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.stellar_times_clipped, self.average_stellar_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_spectra_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.spectra_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_spectra_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.spectra_times_clipped, self.average_spectra_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_dust_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.dust_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_dust_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.dust_times_clipped, self.average_dust_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_writing_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.writing_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_writing_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.writing_times_clipped, self.average_writing_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_waiting_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.waiting_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_waiting_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.waiting_times_clipped, self.average_waiting_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_communication_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.communication_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_communication_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.communication_times_clipped, self.average_communication_times)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_intermediate_times(self):

        """
        This function ...
        :return:
        """

        return self.get_average_times(self.intermediate_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_intermediate_times(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_times(self.intermediate_times_clipped, self.average_intermediate_times)

    # -----------------------------------------------------------------

    def show_runtimes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, **kwargs)

        # Show
        if host is None: self.show_runtimes()
        elif parallelization is None: self.show_runtimes_host(host)
        else: self.show_runtimes_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def show_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing runtimes ...")

        # Loop over the hosts
        for host in self.average_total_times:

            # Show
            self.show_runtimes_host(host)

    # -----------------------------------------------------------------

    def show_runtimes_host(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        # Loop over the parallelization schemes
        for parallelization in self.average_total_times[host]:

            # Show
            self.show_runtimes_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def show_runtimes_host_parallelization(self, host, parallelization):

        """
        Thisf unction ...
        :param host:
        :param parallelization:
        :return:
        """

        # Debugging
        log.debug("Showing the runtimes for host '" + tostr(host) + "' and parallelization scheme '" + tostr(parallelization) + "' ...")

        # Get average runtimes
        total = self.average_total_times[host][parallelization]
        setup = self.average_setup_times[host][parallelization]
        stellar = self.average_stellar_times[host][parallelization]
        spectra = self.average_spectra_times[host][parallelization]
        dust = self.average_dust_times[host][parallelization]
        writing = self.average_writing_times[host][parallelization]
        waiting = self.average_waiting_times[host][parallelization]
        communication = self.average_communication_times[host][parallelization]
        intermediate = self.average_intermediate_times[host][parallelization]

        # Get standard deviations
        total_err = self.stddev_total_times[host][parallelization]
        setup_err = self.stddev_setup_times[host][parallelization]
        stellar_err = self.stddev_stellar_times[host][parallelization]
        spectra_err = self.stddev_spectra_times[host][parallelization]
        dust_err = self.stddev_dust_times[host][parallelization]
        writing_err = self.stddev_writing_times[host][parallelization]
        waiting_err = self.stddev_waiting_times[host][parallelization]
        communication_err = self.stddev_communication_times[host][parallelization]
        intermediate_err = self.stddev_intermediate_times[host][parallelization]

        # Get numer of measurements
        ntotal = self.ntotal_times[host][parallelization]
        nsetup = self.nsetup_times[host][parallelization]
        nstellar = self.nstellar_times[host][parallelization]
        nspectra = self.nspectra_times[host][parallelization]
        ndust = self.ndust_times[host][parallelization]
        nwriting = self.nwriting_times[host][parallelization]
        nwaiting = self.nwaiting_times[host][parallelization]
        ncommunication = self.ncommunication_times[host][parallelization]
        nintermediate = self.nintermediate_times[host][parallelization]

        # Get number of outliers
        total_noutliers = self.ntotal_times_outliers[host][parallelization]
        setup_noutliers = self.nsetup_times_outliers[host][parallelization]
        stellar_noutliers = self.nstellar_times_outliers[host][parallelization]
        spectra_noutliers = self.nspectra_times_outliers[host][parallelization]
        dust_noutliers = self.ndust_times_outliers[host][parallelization]
        writing_noutliers = self.nwriting_times_outliers[host][parallelization]
        waiting_noutliers = self.nwaiting_times_outliers[host][parallelization]
        communication_noutliers = self.ncommunication_times_outliers[host][parallelization]
        intermediate_noutliers = self.nintermediate_times_outliers[host][parallelization]

        # Show
        print("")
        print(fmt.bold + "Runtimes:" + fmt.reset)
        print("")
        print(" - Total time: (" + tostr(total, round=True, ndigits=3) + " ± " + tostr(total_err, round=True, ndigits=3) + ") minutes [" + str(total_noutliers) + " outliers out of " + str(ntotal) + " data points]")
        print(" - Setup time: (" + tostr(setup, round=True, ndigits=3) + " ± " + tostr(setup_err, round=True, ndigits=3) + ") minutes [" + str(setup_noutliers) + " outliers out of " + str(nsetup) + " data points]")
        print(" - Stellar time: (" + tostr(stellar, round=True, ndigits=3) + " ± " + tostr(stellar_err, round=True, ndigits=3) + ") minutes [" + str(stellar_noutliers) + " outliers out of " + str(nstellar) + " data points]")
        print(" - Spectra time: (" + tostr(spectra, round=True, ndigits=3) + " ± " + tostr(spectra_err, round=True, ndigits=3) + ") minutes [" + str(spectra_noutliers) + " outliers out of " + str(nspectra) + " data points]")
        print(" - Dust time: (" + tostr(dust, round=True, ndigits=3) + " ± " + tostr(dust_err, round=True, ndigits=3) + ") minutes [" + str(dust_noutliers) + " outliers out of " + str(ndust) + " data points]")
        print(" - Writing time: (" + tostr(writing, round=True, ndigits=3) + " ± " + tostr(writing_err, round=True, ndigits=3) + ") minutes [" + str(writing_noutliers) + " outliers out of " + str(nwriting) + " data points]")
        print(" - Waiting time: (" + tostr(waiting, round=True, ndigits=3) + " ± " + tostr(waiting_err, round=True, ndigits=3) + ") minutes [" + str(waiting_noutliers) + " outliers out of " + str(nwaiting) + " data points]")
        print(" - Communication time: (" + tostr(communication, round=True, ndigits=3) + " ± " + tostr(communication_err, round=True, ndigits=3) + ") minutes [" + str(communication_noutliers) + " outliers out of " + str(ncommunication) + " data points]")
        print(" - Intermediate time: (" + tostr(intermediate, round=True, ndigits=3) + " ± " + tostr(intermediate_err, round=True, ndigits=3) + ") minutes [" + str(intermediate_noutliers) + " outliers out of " + str(nintermediate) + " data points]")

    # -----------------------------------------------------------------

    def get_average_memories(self, memories_dict):

        """
        This function ...
        :param memories_dict:
        :return:
        """

        averages = defaultdict(dict)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                averages[host][parallelization] = numbers.arithmetic_mean(*memories_dict[host][parallelization])
        return averages

    # -----------------------------------------------------------------

    def get_stddev_memories(self, memories_dict, averages=None):

        """
        This function ...
        :param memories_dict:
        :param averages:
        :return:
        """

        stddevs = defaultdict(dict)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                if averages is not None and host in averages and parallelization in averages[host]: average = averages[host][parallelization]
                else: average = None
                stddevs[host][parallelization] = numbers.standard_deviation(*memories_dict[host][parallelization], mean=average)
        return stddevs

    # -----------------------------------------------------------------

    @lazyproperty
    def average_total_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.total_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_total_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.total_memories_clipped, self.average_total_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_setup_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.setup_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_setup_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.setup_memories_clipped, self.average_setup_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_stellar_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.stellar_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_stellar_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.stellar_memories_clipped, self.average_stellar_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_spectra_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.spectra_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_spectra_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.spectra_memories_clipped, self.average_spectra_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_dust_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.dust_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_dust_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.dust_memories_clipped, self.average_dust_memories)

    # -----------------------------------------------------------------

    @lazyproperty
    def average_writing_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_average_memories(self.writing_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stddev_writing_memories(self):

        """
        This function ...
        :return:
        """

        return self.get_stddev_memories(self.writing_memories_clipped, self.average_writing_memories)

    # -----------------------------------------------------------------

    def show_memory_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, **kwargs)

        # Show
        if host is None: self.show_memory()
        elif parallelization is None: self.show_memory_host(host)
        else: self.show_memory_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def show_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing memory usage ...")

        # Loop over the hosts
        for host in self.average_total_memories:

            # Show
            self.show_memory_host(host)

    # -----------------------------------------------------------------

    def show_memory_host(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        # Loop over the parallelization schemes
        for parallelization in self.average_total_memories[host]:

            # Show
            self.show_memory_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def show_memory_host_parallelization(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Debugging
        log.debug("Showing the memory usage for host '" + tostr(host) + "' and parallelization scheme '" + tostr(parallelization) + "' ...")

        # Get average values
        total = self.average_total_memories[host][parallelization]
        setup = self.average_setup_memories[host][parallelization]
        stellar = self.average_stellar_memories[host][parallelization]
        spectra = self.average_spectra_memories[host][parallelization]
        dust = self.average_dust_memories[host][parallelization]
        writing = self.average_writing_memories[host][parallelization]

        # Get stddevs
        total_err = self.stddev_total_memories[host][parallelization]
        setup_err = self.stddev_setup_memories[host][parallelization]
        stellar_err = self.stddev_stellar_memories[host][parallelization]
        spectra_err = self.stddev_spectra_memories[host][parallelization]
        dust_err = self.stddev_dust_memories[host][parallelization]
        writing_err = self.stddev_writing_memories[host][parallelization]

        # Get number of measurements
        ntotal = self.ntotal_memories[host][parallelization]
        nsetup = self.nsetup_memories[host][parallelization]
        nstellar = self.nstellar_memories[host][parallelization]
        nspectra = self.nspectra_memories[host][parallelization]
        ndust = self.ndust_memories[host][parallelization]
        nwriting = self.nwriting_memories[host][parallelization]

        # Get number of outliers
        total_noutliers = self.ntotal_memories_outliers[host][parallelization]
        setup_noutliers = self.nsetup_memories_outliers[host][parallelization]
        stellar_noutliers = self.nstellar_memories_outliers[host][parallelization]
        spectra_noutliers = self.nspectra_memories_outliers[host][parallelization]
        dust_noutliers = self.ndust_memories_outliers[host][parallelization]
        writing_noutliers = self.nwriting_memories_outliers[host][parallelization]

        # Show
        print("")
        print(fmt.bold + "Memory usage:" + fmt.reset)
        print("")
        print(" - Total memory: (" + tostr(total, round=True, ndigits=3) + " ± " + tostr(total_err, round=True, ndigits=3) + ") GB [" + str(total_noutliers) + " outliers out of " + str(ntotal) + " data points]")
        print(" - Setup memory: (" + tostr(setup, round=True, ndigits=3) + " ± " + tostr(setup_err, round=True, ndigits=3) + ") GB [" + str(setup_noutliers) + " outliers out of " + str(nsetup) + " data points]")
        print(" - Stellar memory: (" + tostr(stellar, round=True, ndigits=3) + " ± " + tostr(stellar_err, round=True,ndigits=3) + ") GB [" + str(stellar_noutliers) + " outliers out of " + str(nstellar) + " data points]")
        print(" - Spectra memory: (" + tostr(spectra, round=True, ndigits=3) + " ± " + tostr(spectra_err, round=True, ndigits=3) + ") GB [" + str(spectra_noutliers) + " outliers out of " + str(nspectra) + " data points]")
        print(" - Dust memory: (" + tostr(dust, round=True, ndigits=3) + " ± " + tostr(dust_err, round=True, ndigits=3) + ") GB [" + str(dust_noutliers) + " outliers out of " + str(ndust) + " data points]")
        print(" - Writing memory: (" + tostr(writing, round=True, ndigits=3) + " ± " + tostr(writing_err, round=True, ndigits=3) + ") GB [" + str(writing_noutliers) + " outliers out of " + str(nwriting) + " data points]")

    # -----------------------------------------------------------------

    @property
    def do_write_assignment(self):

        """
        This function ...
        :return:
        """

        if self.config.dry:
            log.warning("[DRY] Not writing new assignment scheme in dry mode")
            return False

        if self.config.write_assignment is not None: return self.config.write_assignment
        else: return self._adapted_assignment or self._new_assignment

    # -----------------------------------------------------------------

    @property
    def do_write_status(self):

        """
        This function ...
        :return:
        """

        # use self.has_status?
        return self.config.write_status

    # -----------------------------------------------------------------

    @property
    def do_write_moved(self):

        """
        This function ...
        :return:
        """

        return self.config.write_moved and self.has_moved

    # -----------------------------------------------------------------

    @property
    def do_write_relaunched(self):

        """
        This function ...
        :return:
        """

        return self.config.write_relaunched and self.has_relaunched

    # -----------------------------------------------------------------

    @property
    def do_write_commands(self):

        """
        This function ...
        :return:
        """

        return self.config.write_commands and self.has_commands

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the assignment scheme
        if self.do_write_assignment: self.write_assignment()

        # Write the status
        if self.do_write_status: self.write_status()

        # Write the moved
        if self.do_write_moved: self.write_moved()

        # Write the relaunched
        if self.do_write_relaunched: self.write_relaunched()

        # Write the commands
        if self.do_write_commands: self.write_commands()

    # -----------------------------------------------------------------

    def write_assignment(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the assignment scheme ...")

        # Determine path
        path = self.output_path_file("assignment.dat")

        # Overwrite the original assignment table
        if self.assignment.path is not None: self.assignment.save()

        # Write the assignment to the output directory
        self.assignment.saveto(path)

    # -----------------------------------------------------------------

    def write_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulation status table ...")

        # Determine path
        path = self.output_path_file("status.dat")

        # Write the table
        self.status.saveto(path)

    # -----------------------------------------------------------------

    def write_moved(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the moved simulations ...")

        # Determine path
        path = self.output_path_file("moved.dat")

        # Write
        self.moved.saveto(path)

    # -----------------------------------------------------------------

    def write_relaunched(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the relaunched simulations ...")

        # Determine path
        path = self.output_path_file("relaunched.dat")

        # Write
        self.relaunched.saveto(path)

    # -----------------------------------------------------------------

    def write_commands(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the commands ...")

        # Determine the path
        path = self.output_path_file("commands.dat")

        # Write
        fs.write_lines(path, self.commands)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Runtimes
        if self.config.plot_runtimes: self.plot_runtimes()

        # Memory
        if self.config.plot_memory: self.plot_memory()

    # -----------------------------------------------------------------

    def get_times_distributions(self, times_dict):

        """
        This function ...
        :param times_dict:
        :return:
        """

        distributions = defaultdict(dict)
        for host in times_dict:
            for parallelization in times_dict[host]:
                distributions[host][parallelization] = Distribution.from_values("Runtime", times_dict[host][parallelization], unit="min")
        return distributions

    # -----------------------------------------------------------------

    @lazyproperty
    def total_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.total_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.setup_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.stellar_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.spectra_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.dust_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.writing_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def waiting_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.waiting_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def communication_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.communication_times_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def intermediate_times_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_times_distributions(self.intermediate_times_clipped)

    # -----------------------------------------------------------------

    def plot_runtimes_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, **kwargs)

        # Plot
        if host is None: self.plot_runtimes()
        elif parallelization is None: self.plot_runtimes_host(host)
        else: self.plot_runtimes_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        # Loop over the hosts
        for host in self.total_times_distributions:

            # Plot for the host
            self.plot_runtimes_host(host)

    # -----------------------------------------------------------------

    def plot_runtimes_host(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        # Loop over the parallelization schemes
        for parallelization in self.total_times_distributions[host]:

            # Plot
            self.plot_runtimes_host_parallelization(host, parallelization, phases=self.config.plot_runtimes_phases)

    # -----------------------------------------------------------------

    def plot_runtimes_host_parallelization(self, host, parallelization, phases=None):

        """
        This function ...
        :param host:
        :param parallelization:
        :param phases:
        :return:
        """

        # Debugging
        log.debug("Plotting the runtimes for host '" + tostr(host) + "' and parallelization scheme '" + tostr(parallelization) + "' ...")

        # Plot for phases
        if phases is None or "total" in phases: self.plot_runtimes_total(host, parallelization)
        if phases is None or "setup" in phases: self.plot_runtimes_setup(host, parallelization)
        if phases is None or "stellar" in phases: self.plot_runtimes_stellar(host, parallelization)
        if phases is None or "spectra" in phases: self.plot_runtimes_spectra(host, parallelization)
        if phases is None or "dust" in phases: self.plot_runtimes_dust(host, parallelization)
        if phases is None or "writing" in phases: self.plot_runtimes_writing(host, parallelization)
        if phases is None or "waiting" in phases: self.plot_runtimes_waiting(host, parallelization)
        if phases is None or "communication" in phases: self.plot_runtimes_communication(host, parallelization)
        if phases is None or "intermediate" in phases: self.plot_runtimes_intermediate(host, parallelization)

    # -----------------------------------------------------------------

    def plot_runtimes_total(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        total = self.total_times_distributions[host][parallelization]

        # Get the average runtime
        mean_total = self.average_total_times[host][parallelization]

        # Plot
        plot_distribution(total, title="Total runtime", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_total)

    # -----------------------------------------------------------------

    def plot_runtimes_setup(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        setup = self.setup_times_distributions[host][parallelization]

        # Get the average runtime
        mean_setup = self.average_setup_times[host][parallelization]

        # Plot
        plot_distribution(setup, title="Setup runtime", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_setup)

    # -----------------------------------------------------------------

    def plot_runtimes_stellar(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        stellar = self.stellar_times_distributions[host][parallelization]

        # Get the average runtime
        mean_stellar = self.average_stellar_times[host][parallelization]

        # Plot
        plot_distribution(stellar, title="Stellar emission runtime", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_stellar)

    # -----------------------------------------------------------------

    def plot_runtimes_spectra(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        spectra = self.spectra_times_distributions[host][parallelization]

        # Get the average runtime
        mean_spectra = self.average_spectra_times[host][parallelization]

        # Plot
        plot_distribution(spectra, title="Spectra calculation runtime", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_spectra)

    # -----------------------------------------------------------------

    def plot_runtimes_dust(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        dust = self.dust_times_distributions[host][parallelization]

        # Get the average runtime
        mean_dust = self.average_dust_times[host][parallelization]

        # Plot
        plot_distribution(dust, title="Dust emission runtime", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_dust)

    # -----------------------------------------------------------------

    def plot_runtimes_writing(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        writing = self.writing_times_distributions[host][parallelization]

        # Get the average runtime
        mean_writing = self.average_writing_times[host][parallelization]

        # Plot
        plot_distribution(writing, title="Writing runtime", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_writing)

    # -----------------------------------------------------------------

    def plot_runtimes_waiting(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        waiting = self.waiting_times_distributions[host][parallelization]

        # Get the average runtime
        mean_waiting = self.average_waiting_times[host][parallelization]

        # Plot
        plot_distribution(waiting, title="Waiting runtime", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_waiting)

    # -----------------------------------------------------------------

    def plot_runtimes_communication(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        communication = self.communication_times_distributions[host][parallelization]

        # Get the average runtime
        mean_communication = self.average_communication_times[host][parallelization]

        # Plot
        plot_distribution(communication, title="Communication runtime", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_communication)

    # -----------------------------------------------------------------

    def plot_runtimes_intermediate(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get the distribution
        intermediate = self.intermediate_times_distributions[host][parallelization]

        # Get the average runtime
        mean_intermediate = self.average_intermediate_times[host][parallelization]

        # Plot
        plot_distribution(intermediate, title="Intermediate runtime", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_intermediate)

    # -----------------------------------------------------------------

    def get_memories_distributions(self, memories_dict):

        """
        This function ...
        :param memories_dict:
        :return:
        """

        distributions = defaultdict(dict)
        for host in memories_dict:
            for parallelization in memories_dict[host]:
                distributions[host][parallelization] = Distribution.from_values("Memory usage", memories_dict[host][parallelization], unit="GB")
        return distributions

    # -----------------------------------------------------------------

    @lazyproperty
    def total_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.total_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def setup_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.setup_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.stellar_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectra_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.spectra_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.dust_memories_clipped)

    # -----------------------------------------------------------------

    @lazyproperty
    def writing_memories_distributions(self):

        """
        This function ...
        :return:
        """

        return self.get_memories_distributions(self.writing_memories_clipped)

    # -----------------------------------------------------------------

    def plot_memory_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, **kwargs)

        # Plot
        if host is None: self.plot_memory()
        elif parallelization is None: self.plot_memory_host(host)
        else: self.plot_memory_host_parallelization(host, parallelization)

    # -----------------------------------------------------------------

    def plot_memory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the memory usage ...")

        # Loop over the hosts
        for host in self.total_memories_distributions:

            # Plot
            self.plot_memory_host(host)

    # -----------------------------------------------------------------

    def plot_memory_host(self, host):

        """
        This function ...
        :param host:
        :return:
        """

        # Loop over the parallelization schemes
        for parallelization in self.total_memories_distributions[host]:

            # Plot
            self.plot_memory_host_parallelization(host, parallelization, phases=self.config.plot_memory_phases)

    # -----------------------------------------------------------------

    def plot_memory_host_parallelization(self, host, parallelization, phases=None):

        """
        This function ...
        :param host:
        :param parallelization:
        :param phases:
        :return:
        """

        # Debugging
        log.debug("Plotting the memory usages for host '" + tostr(host) + "' and parallelization scheme '" + tostr(parallelization) + "' ...")

        # Plot for phases
        if phases is None or "total" in phases: self.plot_memory_total(host, parallelization)
        if phases is None or "setup" in phases: self.plot_memory_setup(host, parallelization)
        if phases is None or "stellar" in phases: self.plot_memory_stellar(host, parallelization)
        if phases is None or "spectra" in phases: self.plot_memory_spectra(host, parallelization)
        if phases is None or "dust" in phases: self.plot_memory_dust(host, parallelization)
        if phases is None or "writing" in phases: self.plot_memory_writing(host, parallelization)

    # -----------------------------------------------------------------

    def plot_memory_total(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        total = self.total_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_total = self.average_total_memories[host][parallelization]

        # Plot
        plot_distribution(total, title="Total memory usage", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_total)

    # -----------------------------------------------------------------

    def plot_memory_setup(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        setup = self.setup_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_setup = self.average_setup_memories[host][parallelization]

        # Plot
        plot_distribution(setup, title="Setup memory usage", x_limits=plot_x_limits, soft_xmin=True, statistics=False,
                          show_mean=True, mean=mean_setup)

    # -----------------------------------------------------------------

    def plot_memory_stellar(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        stellar = self.stellar_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_stellar = self.average_stellar_memories[host][parallelization]

        # Plot
        plot_distribution(stellar, title="Stellar emission memory usage", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_stellar)

    # -----------------------------------------------------------------

    def plot_memory_spectra(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        spectra = self.spectra_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_spectra = self.average_spectra_memories[host][parallelization]

        # Plot
        plot_distribution(spectra, title="Spectra calculation memory usage", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_spectra)

    # -----------------------------------------------------------------

    def plot_memory_dust(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        dust = self.dust_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_dust = self.average_dust_memories[host][parallelization]

        # Plot
        plot_distribution(dust, title="Dust emission memory usage", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_dust)

    # -----------------------------------------------------------------

    def plot_memory_writing(self, host, parallelization):

        """
        This function ...
        :param host:
        :param parallelization:
        :return:
        """

        # Get distribution
        writing = self.writing_memories_distributions[host][parallelization]

        # Get the average memory usage
        mean_writing = self.average_writing_memories[host][parallelization]

        # Plot
        plot_distribution(writing, title="Writing memory usage", x_limits=plot_x_limits, soft_xmin=True,
                          statistics=False, show_mean=True, mean=mean_writing)

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_timeline_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_optional("path", "string", "output path")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def plot_timeline_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Set command name
        name = _plot_command_name + " " + _timeline_command_name

        # Get simulation name
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.plot_timeline_definition, name=name, **kwargs)

        # Plot
        self.plot_timeline(simulation_name, path=config.path)

    # -----------------------------------------------------------------

    def plot_timeline(self, simulation_name, path=None):

        """
        This function ...
        :param simulation_name:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Plotting timeline for simulation '" + simulation_name + "' ...")

        # Get the timeline of the simulation
        timeline = self.get_timeline_for_simulation(simulation_name)

        # Plot the timeline
        plot_timeline(timeline, path=path)

    # -----------------------------------------------------------------

    def reanalyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Re-analysing the simulations ...")

        # Get the config for analysis
        config = self.config.analysis.copy()
        config.batch.replace = True

        # Loop over the analysed simulations
        for simulation in self.all_retrieved_simulations:

            # Check remote
            if self.config.reanalyse_remotes is not None and simulation.host_id not in self.config.reanalyse_remotes: continue

            # Check simulation name
            if self.config.reanalyse_simulations is not None and simulation.name not in self.config.reanalyse_simulations: continue

            # Check whether analysis is already done (so it is actually a re-analysis)
            if not has_analysed(simulation, self.config.reanalyse, self.config.features_reanalysis): continue

            # Set display name
            display_name = self.get_display_name(simulation, add_quotes=True)

            # Re-analyse?
            if self.config.prompt_simulations_reanalysis and not prompt_proceed("re-analyse simulation " + display_name + "?"): continue

            # Re-analyse
            self.reanalyse_simulation(simulation.name, steps=self.config.reanalyse, features=self.config.features_reanalysis, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def reanalyse_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
        definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")
        definition.add_optional("not_steps", "string_list", "don't analyse these steps", choices=all_steps)
        definition.add_optional("not_features", "string_list", "don't analyse these features (if a single not_step is defined)")
        definition.import_section("analysis", "options for the simulation analyser", analyse_simulation_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def reanalyse_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, command_definition=self.reanalyse_simulation_definition, **kwargs)
        steps = config.steps
        features = config.features
        not_steps = config.not_steps
        not_features = config.not_features

        # Reanalyse the simulation
        self.reanalyse_simulation(simulation_name, steps=steps, features=features, not_steps=not_steps, not_features=not_features, config=config.analysis)

    # -----------------------------------------------------------------

    def reanalyse_simulation(self, simulation_name, steps=None, features=None, not_steps=None, not_features=None, config=None):

        """
        This function ...
        :param simulation_name:
        :param steps:
        :param features:
        :param not_steps:
        :param not_features:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Re-analysing simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Reanalyse simulation
        reanalyse_simulation(simulation, steps, features, not_steps=not_steps, not_features=not_features, config=config)

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the simulations ...")

        # Loop over the retrieved simulations
        for simulation in self.all_retrieved_simulations:

            # Check remote
            if self.config.analyse_remotes is not None and simulation.host_id not in self.config.analyse_remotes: continue

            # Check simulation name
            if self.config.analyse_simulations is not None and simulation.name not in self.config.analyse_simulations: continue

            # Check if already analysed
            if simulation.analysed: continue

            # Set display name
            display_name = self.get_display_name(simulation, add_quotes=True)

            # Re-analyse?
            if self.config.prompt_simulations_analysis and not prompt_proceed("analyse simulation " + display_name + "?"): continue

            # Show steps that will be performed
            show_analysis_steps(simulation)

            # Analyse
            self.analyse_simulation(simulation.name, config=self.config.analysis)

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_simulation_definition(self):

        """
        This function ...
        :return:
        """

        return analyse_simulation_definition.copy(pos_optional=False)

    # -----------------------------------------------------------------

    def analyse_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, command_definition=self.analyse_simulation_definition, **kwargs)

        # Analyse
        self.analyse_simulation(simulation_name, config=config)

    # -----------------------------------------------------------------

    def analyse_simulation(self, simulation_name, config=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Analysing simulation " + simulation_name + " ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Analyse the simulation
        analyse_simulation(simulation, config=config)

    # -----------------------------------------------------------------

    def cache(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Caching output of the simulations ...")

        # Loop over the retrieved simulations
        for simulation in self.all_retrieved_simulations:

            # Cache output
            if self.config.cache_output: self.cache_simulation_output(simulation.name)

            # Cache datacubes
            if self.config.cache_datacubes: self.cache_simulation_datacubes(simulation.name)

            # Cache misc output
            if self.config.cache_misc and simulation.analysed_all_misc: self.cache_simulation_misc(simulation.name)

            # Cache images
            if self.config.cache_images and simulation.analysed_all_misc: self.cache_simulation_images(simulation.name)

# -----------------------------------------------------------------

def get_noutliers(nvalues_dict, nclipped_dict):

    """
    This function ...
    :param nvalues_dict:
    :param nclipped_dict:
    :return:
    """

    noutliers = defaultdict(dict)
    for host in nvalues_dict:
        for parallelization in nvalues_dict[host]:
            noutliers = nvalues_dict[host][parallelization] - nclipped_dict[host][parallelization]
    return noutliers

# -----------------------------------------------------------------
