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
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import InteractiveConfigurable, InvalidCommandError
from ..basics.configuration import ConfigurationDefinition, parse_arguments, prompt_choice, prompt_settings, prompt_yn
from ..basics.log import log
from ..tools.stringify import tostr
from .tables import SimulationAssignmentTable, SimulationStatusTable
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
from ..remote.host import load_host
from ..basics.containers import create_nested_defaultdict, create_subdict
from ..tools import sequences
from ..simulation.remote import is_running_status, finished_name, is_invalid_or_unknown_status
from ..simulation.remote import is_analysing_or_analysed_status, is_retrieved_status, is_finished_status
from ..basics.configuration import prompt_variable, prompt_string
from ..remote.host import find_host_ids
from ..simulation.shower import show_simulation, show_analysis, compare_simulations, compare_analysis
from ..simulation.shower import select_simulation_settings, select_analysis_options
from ..simulation.adapter import adapt_simulation, adapt_analysis, adapt_simulations, adapt_analysis_simulations
from ..config.show_simulation_settings import definition as show_simulation_definition
from ..config.adapt_simulation_settings import definition as adapt_simulations_definition
from ..config.show_analysis_options import definition as show_analysis_definition
from ..config.adapt_analysis_options import definition as adapt_analysis_definition
from .batch import BatchLauncher
from ..basics.table import SmartTable
from ..tools import tables
from ..tools import time
from ..config.analyse_simulation import definition as analyse_simulation_definition
from .analyser import all_steps, all_steps_and_extra
from .options import LoggingOptions, SchedulingOptions, AnalysisOptions
from ..remote.load import show_status
from ..remote.mounter import mount_remote
from ..simulation.skifile import show_normalizations
from ..prep.summarize import show_instrument, show_stellar_component, show_dust_component
from ..units.unit import get_common_unit
from ..plot.sed import plot_seds
from ..data.sed import load_multiple_seds
from ..misc.fluxes import get_sed_instrument_name
from ..misc.images import get_datacube_instrument_name
from ..simulation.status import show_log_summary
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
from ..remote.ensemble import SKIRTRemotesEnsemble
from ..simulation.simulation import SkirtSimulation, RemoteSimulation
from ..simulation.definition import SingleSimulationDefinition
from ..simulation.remote import get_simulation_for_host
from ..tools import introspection
from ..simulation.logfile import LogFile
from ..simulation.remote import aborted_name, crashed_name, cancelled_name
from ..simulation.remote import analysed_name, retrieved_name, finished_name

# -----------------------------------------------------------------

# Define types
datacube_types = [output_type_names.total_images, output_type_names.count_images, output_type_names.direct_images, output_type_names.transparent_images, output_type_names.scattered_images, output_type_names.dust_images, output_type_names.dust_scattered_images]
image_types = [misc_type_names.images, misc_type_names.images_for_fluxes]

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

plot_x_limits = (0., None)

# -----------------------------------------------------------------

clear_analysis_steps = ["extraction", "plotting", "misc"]

# -----------------------------------------------------------------

failed_name = "failed" # all failed simulations
exceeded_walltime_name = "exceeded_walltime"
exceeded_memory_name = "exceeded_memory"

# Cases of failed simulations
failed_cases = [failed_name, aborted_name, crashed_name, cancelled_name, exceeded_walltime_name, exceeded_memory_name]

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
_unlaunch_command_name = "unlaunch"
_unretrieve_command_name = "unretrieve"
_unanalyse_command_name = "unanalyse"
_relaunch_command_name = "relaunch"
_retry_command_name = "retry"
_log_command_name = "log"
_error_command_name = "error"
_settings_command_name = "settings"
_analysis_command_name = "analysis"
_steps_command_name = "steps"
_adapt_command_name = "adapt"
_compare_command_name = "compare"
_steal_command_name = "steal"
_retrieve_command_name = "retrieve"
_analyse_command_name = "analyse"
_allretrieve_command_name = "allretrieve"
_allanalyse_command_name = "allanalyse"
_reanalyse_command_name = "reanalyse"
_mimic_command_name = "mimic"
_launch_command_name = "launch"
_assignment_command_name = "assignment"
_runtimes_command_name = "runtimes"
_memory_command_name = "memory"
_cpu_command_name = "cpu"
_timeline_command_name = "timeline"
_scaling_command_name = "scaling"
_base_command_name = "base"
_simulation_command_name = "simulation"

# -----------------------------------------------------------------

# Define show commands
show_commands = OrderedDict()
show_commands.description = "show"
show_commands[_assignment_command_name] = ("show_assignment", False, "show the simulation assignment scheme", None)
show_commands[_status_command_name] = ("show_status", False, "show the simulation status", None)
show_commands[_runtimes_command_name] = ("show_runtimes_command", True, "show the simulation runtimes", "host_parallelization")
show_commands[_memory_command_name] = ("show_memory_command", True, "show the simulation memory usages", "host_parallelization")
show_commands[_cpu_command_name] = ("show_cpu_command", True, "show the CPU usage of the simulations", "simulations")
show_commands[_timeline_command_name] = ("show_timeline_command", True, "show the timeline of a simulation", "simulation")

# -----------------------------------------------------------------

# Define plot commands
plot_commands = OrderedDict()
plot_commands.description = "plot"
plot_commands[_runtimes_command_name] = ("plot_runtimes_command", True, "plot simulation runtimes", "host_parallelization")
plot_commands[_memory_command_name] = ("plot_memory_command", True, "plot simulation memory usages", "host_parallelization")
plot_commands[_timeline_command_name] = ("plot_timeline_command", True, "plot simulation timeline", "simulation")
plot_commands[_scaling_command_name] = ("plot_scaling_command", True, "plot scaling of simulations", None)

# -----------------------------------------------------------------

# Define open commands
open_commands = OrderedDict()
open_commands.description = "open input, output or base simulation directory"
open_commands[_base_command_name] = ("open_base_command", True, "open simulation base directory", "simulation")
open_commands[_input_command_name] = ("open_input_command", True, "open simulation input directory", "simulation")
open_commands[_output_command_name] = ("open_output_command", True, "open simulation output directory", "simulation")
open_commands[_extraction_command_name] = ("open_extraction_command", True, "open simulation extraction output directory", "simulation")
open_commands[_plotting_command_name] = ("open_plotting_command", True, "open simulation plotting output directory", "simulation")
open_commands[_misc_command_name] = ("open_misc_command", True, "open simulation miscellaneous output directory", "simulation")

# -----------------------------------------------------------------

# Define adapt commands
adapt_commands = OrderedDict()
adapt_commands.description = "adapt simulation settings or analysis options"
adapt_commands[_simulation_command_name] = ("adapt_simulation_settings_command", True, "adapt simulation settings", "simulations")
adapt_commands[_analysis_command_name] = ("adapt_analysis_options_command", True, "adapt analysis options", "simulations")

# -----------------------------------------------------------------

# Define compare commands
compare_commands = OrderedDict()
compare_commands.description = "compare simulation settings or analysis options between simulations"
compare_commands[_simulation_command_name] = ("compare_simulation_settings_command", True, "compare simulation settings", "simulations")
compare_commands[_analysis_command_name] = ("compare_analysis_options_command", True, "compare analysis options", "simulations")

# -----------------------------------------------------------------

# Define steal commands
steal_commands = OrderedDict()
steal_commands.description = "take on simulation settings or analysis options from one simulation to another"
steal_commands[_simulation_command_name] = ("steal_simulation_settings_command", True, "steal simulation settings", "simulations_simulation")
steal_commands[_analysis_command_name] = ("steal_analysis_options_command", True, "steal analysis options", "simulations_simulation")

# -----------------------------------------------------------------

# Define clear commands
clear_commands = OrderedDict()
clear_commands.subject = "simulation"
clear_commands.description = "clear simulation output/input/analysis"
clear_commands[_input_command_name] = ("clear_simulation_input_command", True, "clear simulation input", "simulation")
clear_commands[_output_command_name] = ("clear_simulation_output_command", True, "clear simulation output", "simulation")
clear_commands[_analysis_command_name] = ("clear_simulation_analysis_command", True, "clear simulation analysis output", "simulation")
clear_commands[_extraction_command_name] = ("clear_extraction_command", True, "clear simulation extraction output", "simulation")
clear_commands[_plotting_command_name] = ("clear_plotting_command", True, "clear simulation plotting output", "simulation")
clear_commands[_misc_command_name] = ("clear_misc_command", True, "clear simulation misc output", "simulation")

# -----------------------------------------------------------------

# Define cache commands
cache_commands = OrderedDict()
cache_commands.subject = "simulations"
cache_commands.description = "cache simulation/analysis output to another directory/filesystem"
cache_commands[_output_command_name] = ("cache_simulations_output_command", True, "cache simulation output", "simulations")
cache_commands[_analysis_command_name] = ("cache_simulations_analysis_command", True, "cache simulation analysis output", "simulations")
cache_commands[_extraction_command_name] = ("cache_simulations_extraction_command", True, "cache simulation extraction output", "simulations")
cache_commands[_plotting_command_name] = ("cache_simulations_plotting_command", True, "cache simulation plotting output", "simulations")
cache_commands[_misc_command_name] = ("cache_simulations_misc_command", True, "cache simulation misc output", "simulations")
cache_commands[_datacube_command_name] = ("cache_simulations_datacubes_command", True, "cache datacubes", "simulations")
cache_commands[_images_command_name] = ("cache_simulations_images_command", True, "cache images", "simulations")

# -----------------------------------------------------------------

# Define commands
commands = OrderedDict()
commands[_help_command_name] = ("show_help", False, "show help", None)
commands[_history_command_name] = ("show_history_command", True, "show history of executed commands", None)
commands[_status_command_name] = ("show_status_command", True, "show simulations status", None)
commands[_find_command_name] = ("find_simulations_command", True, "find particular simulation(s) based on info or other parameters", None)
commands[_hosts_command_name] = ("show_hosts_command", True, "show remote hosts of the simulations", "hosts")
commands[_parallelizations_command_name] = ("show_parallelizations_command", True, "show parallelization schemes used per host", "host")
commands[_info_command_name] = ("show_info_command", True, "show info about the simulation (if defined in info tables)", "simulation")
commands[_open_command_name] = open_commands #(None, None, "open input, output or base simulation directory", "simulation")
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
commands[_show_command_name] = show_commands
commands[_plot_command_name] = plot_commands
commands[_move_command_name] = ("move_simulations_command", True, "move simulations from one remote to another", "simulations")
commands[_stop_command_name] = ("stop_simulations_command", True, "stop running simulations", "simulation")
commands[_cancel_command_name] = ("cancel_simulations_command", True, "cancel queued simulations", "simulations")
commands[_remove_command_name] = ("remove_simulation_command", True, "remove simulation", "simulations")
commands[_clear_command_name] = clear_commands #(None, None, "clear simulation output/input/analysis", "simulation")
commands[_cache_command_name] = cache_commands #(None, None, "cache simulation/analysis output to another directory/filesystem", "simulations")
commands[_unfinish_command_name] = ("unfinish_simulation_command", True, "remove remote output of a simulation and unset finished flag", "simulation")
commands[_unlaunch_command_name] = ("unlaunch_simulation_command", True, "undo launching a simulation", "simulation")
commands[_unretrieve_command_name] = ("unretrieve_simulation_command", True, "remove local (retrieved) output of a simulation and unset retrieved flag", "simulation")
commands[_unanalyse_command_name] = ("unanalyse_simulation_command", True, "remove analysis output of a simulation and unset analysis flags", "simulation")
commands[_relaunch_command_name] = ("relaunch_simulation_command", True, "relaunch simulations on the original remote host", "simulation")
commands[_retry_command_name] = ("retry_simulations_command", True, "try launching certain failed simulations again", None)
commands[_log_command_name] = ("show_simulation_log_command", True, "show log output of a simulation", "simulation")
commands[_error_command_name] = ("show_simulation_errors_command", True, "show error output of a simulation", "simulation")
commands[_settings_command_name] = ("show_simulation_settings_command", True, "show simulation settings", "simulation")
commands[_analysis_command_name] = ("show_analysis_options_command", True, "show analysis options", "simulation")
commands[_steps_command_name] = ("show_analysis_steps_command", True, "show analysis steps", "simulation")
commands[_adapt_command_name] = adapt_commands #(None, None, "adapt simulation settings or analysis options", "simulations")
commands[_compare_command_name] = compare_commands #(None, None, "compare simulation settings or analysis options between simulations", "simulations")
commands[_steal_command_name] = steal_commands #(None, None, "take on simulation settings or analysis options from one simulation to another", "simulations_simulation")
commands[_retrieve_command_name] = ("retrieve_simulations_command", True, "retrieve one or multiple simulation(s) from the remote host", "simulations")
commands[_analyse_command_name] = ("analyse_simulations_command", True, "analyse one or multiple simulation(s)", "simulations")
commands[_allretrieve_command_name] = ("retrieve_all", False, "retrieve all finished simulations", None) # -> REPLACE BY USING 'RETRIEVE ALL'? -> NO, (not yet), this function loops over all simulations that CAN be retrieved
commands[_allanalyse_command_name] = ("analyse_all_command", True, "analyse all retrieved simulations", None) # -> REPLACE BY USING 'ANALYSE ALL'? -> NO, (not yet), this function loops over all simulations that CAN be analysed
commands[_reanalyse_command_name] = ("reanalyse_simulations_command", True, "re-analyse a simulation", "simulations")
commands[_mimic_command_name] = ("mimic_simulation_command", True, "mimic a simulation", "simulation")
commands[_launch_command_name] = ("launch_simulation_command", True, "launch a new simulation", None)

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

class NewSimulationsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Simulation name"] = (str, None, "name of the simulation")
    _column_info["Ski filepath"] = (str, None, "path of the skifile")
    _column_info["Host ID"] = (str, None, "host ID")
    _column_info["ID"] = (str, None, "simulation ID")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(NewSimulationsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):
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

    def add_simulation(self, name, ski_path, host_id=None, simulation_id=None):

        """
        This function ...
        :param name:
        :param ski_path:
        :param host_id:
        :param simulation_id:
        :return:
        """

        # Set values
        values = [name, ski_path, host_id, simulation_id]

        # Add the row
        self.add_row(values)

    # -----------------------------------------------------------------

    def set_id(self, simulation_name, id):

        """
        This function ...
        :param simulation_name:
        :param id:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        self.set_value("ID", index, str(id))

# -----------------------------------------------------------------

# Extra info
_screen_extra_name = "screen"
_job_extra_name = "job"
_disk_extra_name = "disk"
_runtime_extra_name = "runtime"
_memory_extra_name = "memory"
_elapsed_extra_name = "elapsed"
_cached_output_name = "cached_out"
_cached_extraction_name = "cached_extr"
_cached_plotting_name = "cached_plot"
_cached_misc_name = "cached_misc"

# -----------------------------------------------------------------

# Define extra columns
extra_columns = OrderedDict()
extra_columns[_screen_extra_name] = "screen session name"
extra_columns[_job_extra_name] = "job ID"
extra_columns[_disk_extra_name] = "simulation output disk size"
extra_columns[_runtime_extra_name] = "total simulation runtime"
extra_columns[_memory_extra_name] = "peak simulation memory usage"
extra_columns[_elapsed_extra_name] = "elapsed time since start of simulation"
extra_columns[_cached_output_name] = "whether the output of the simulation has been cached"
extra_columns[_cached_extraction_name] = "whether the extraction output has been cached"
extra_columns[_cached_plotting_name] = "whether the plotting output has been cached"
extra_columns[_cached_misc_name] = "whether the miscellaneous output has been cached"

# Define extra column names
extra_column_names = dict()
extra_column_names[_screen_extra_name] = "Screen name"
extra_column_names[_job_extra_name] = "Job ID"
extra_column_names[_disk_extra_name] = "Disk size"
extra_column_names[_runtime_extra_name] = "Runtime"
extra_column_names[_memory_extra_name] = "Peak memory"
extra_column_names[_elapsed_extra_name] = "Elapsed time"
extra_column_names[_cached_output_name] = "Cached output"
extra_column_names[_cached_extraction_name] = "Cached extraction"
extra_column_names[_cached_plotting_name] = "Cached plotting"
extra_column_names[_cached_misc_name] = "Cached misc output"

# Define extra column units
extra_column_units = dict()
extra_column_units[_disk_extra_name] = "GB"
extra_column_units[_runtime_extra_name] = "h"
extra_column_units[_memory_extra_name] = "GB"
extra_column_units[_elapsed_extra_name] = "h"

# -----------------------------------------------------------------

class SimulationManager(InteractiveConfigurable):

    """
    This class ...
    """

    _commands = commands
    _log_section = "SIMULATION MANAGER"

    # -----------------------------------------------------------------

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

        # The new simulations
        self.new = NewSimulationsTable()

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
    def has_assignment(self):

        """
        This function ...
        :return:
        """

        return self.assignment is not None

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
        return self.config.move

    # -----------------------------------------------------------------

    @property
    def has_showing(self):
        return self.config.show and (self.config.show_assignment or self.config.show_status or self.config.show_runtimes or self.config.show_memory)

    # -----------------------------------------------------------------

    @property
    def has_plotting(self):
        return self.config.plot and (self.config.plot_runtimes or self.config.plot_memory)

    # -----------------------------------------------------------------

    @property
    def has_writing(self):
        return self.config.write and (self.config.write_assignment or self.config.write_status or self.config.write_moved or self.config.write_relaunched or self.config.write_commands)

    # -----------------------------------------------------------------

    @property
    def has_analysis(self):
        return self.config.analyse

    # -----------------------------------------------------------------

    @property
    def has_reanalysis(self):
        return self.config.reanalysis is not None

    # -----------------------------------------------------------------

    @property
    def has_any(self):
        return self.has_moving or self.has_showing or self.has_plotting or self.has_writing or self.has_analysis or self.has_reanalysis

    # -----------------------------------------------------------------

    @property
    def do_commands(self):
        return self.config.commands is not None and len(self.config.commands) > 0

    # -----------------------------------------------------------------

    @property
    def do_interactive(self):
        if self.config.interactive is None: return not self.has_any and not self.do_commands
        else: return self.config.interactive

    # -----------------------------------------------------------------

    @property
    def do_moving(self):
        return self.config.move

    # -----------------------------------------------------------------

    @property
    def do_launch(self):
        return self.launcher.has_queued

    # -----------------------------------------------------------------

    @property
    def do_showing(self):
        return self.config.show

    # -----------------------------------------------------------------

    @property
    def do_writing(self):
        return self.config.write

    # -----------------------------------------------------------------

    @property
    def do_plotting(self):
        return self.config.plot

    # -----------------------------------------------------------------

    @property
    def do_reanalysis(self):
        return self.config.reanalyse is not None

    # -----------------------------------------------------------------

    @property
    def do_analysis(self):
        return self.config.analyse

    # -----------------------------------------------------------------

    @property
    def do_caching(self):
        return self.config.cache_output or self.config.cache_datacubes or self.config.cache_misc or self.config.cache_images

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

        # 12. Write the history
        if self.has_commands: self.write_history()

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
        if kwargs.get("status", None) is not None: self.status = kwargs.pop("status")
        elif self.config.status is not None: self.status = SimulationStatusTable.from_file(self.config.status)

        # Initialize simulations and assignment scheme
        self.initialize(**kwargs)

        # Get timing table
        self.get_timing_table(**kwargs)

        # Get memory table
        self.get_memory_table(**kwargs)

        # Get the info tables
        self.set_info(**kwargs)

        # Set the remotes
        self.set_remotes(**kwargs)

        # Set launcher options
        self.set_launcher_options()

    # -----------------------------------------------------------------

    def set_info(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Setting info tables for the simulations ...")

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

    # -----------------------------------------------------------------

    def set_remotes(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Setting the remotes ...")

        # Get remotes
        if kwargs.get("remotes", None) is not None:

            remotes = kwargs.pop("remotes")

            # Input is already an ensemble
            if isinstance(remotes, SKIRTRemotesEnsemble):
                for name in remotes.names:
                    if remotes.is_used(name): self.set_remote(remotes[name], name)
                    else: self.add_host_id(remotes.get_host_id(name), name=name)

            # Regular dict
            elif types.is_dictionary(remotes):
                for name in remotes: self.set_remote(remotes[name], name)

            # Sequence
            elif types.is_sequence(remotes):
                for remote in remotes: self.set_remote(remote)

            # Invalid
            else: raise ValueError("Invalid type for 'remotes'")

        # Set remotes from config
        if self.config.remotes is not None: self.add_host_ids(self.config.remotes)

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
    def all_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_queued_simulations(self):
        
        """
        This function ...
        :return: 
        """

        for simulation in self.all_simulations:
            if not self.is_queued(simulation.name): continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_queued_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_queued_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_running_simulations(self):
        
        """
        This function ...
        :return: 
        """

        for simulation in self.all_simulations:
            if not self.is_running(simulation.name): continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_running_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_running_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_queued_or_running_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if not (self.is_queued(simulation.name) or self.is_running(simulation.name)): continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_queued_or_running_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_queued_or_running_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_finished_simulations(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_simulations:
            if not simulation.finished: continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_finished_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_finished_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_not_finished_simulations(self):
        
        """
        This function ...
        :return: 
        """
        
        for simulation in self.all_simulations:
            if simulation.finished: continue
            yield simulation
        
    # -----------------------------------------------------------------

    @property
    def all_not_finished_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_not_finished_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_finished_not_retrieved_simulations(self):
        
        """
        Thisn function ...
        :return: 
        """
        
        for simulation in self.all_simulations:
            if not simulation.finished: continue
            if simulation.retrieved: continue
            yield simulation

    # -----------------------------------------------------------------

    @property
    def all_finished_not_retrieved_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_finished_not_retrieved_simulations: yield simulation.name

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
    def all_retrieved_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_retrieved_simulations: yield simulation.name

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
    def all_not_retrieved_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_not_retrieved_simulations: yield simulation.name

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
    def all_analysed_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_analysed_simulations: yield simulation.name

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

    @property
    def all_not_analysed_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_not_analysed_simulations: yield simulation.name

    # -----------------------------------------------------------------

    @property
    def all_retrieved_not_analysed_simulations(self):

        """
        This function ...
        :return:
        """
        
        for simulation in self.all_retrieved_simulations:
            if not simulation.analysed: yield simulation

    # -----------------------------------------------------------------

    @property
    def all_retrieved_not_analysed_simulation_names(self):

        """
        This function ...
        :return:
        """

        for simulation in self.all_retrieved_not_analysed_simulations: yield simulation.name

    # -----------------------------------------------------------------

    def get_timing_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the timing table ...")

        if "timing" in kwargs: self.timing = kwargs.pop("timing")
        elif self.config.timing is not None: self.timing = TimingTable.from_file(self.config.timing)
        else:

            table_paths = []
            for simulation in self.all_simulations: table_paths.append(simulation.analysis.timing_table_path)
            if not sequences.all_none(table_paths):
                timing_table_path = sequences.get_all_equal_value(table_paths, ignore_none=True, return_none=True)
                if timing_table_path is not None: self.timing = TimingTable.from_file(timing_table_path)

    # -----------------------------------------------------------------

    def get_memory_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the memory table ...")

        if "memory" in kwargs: self.memory = kwargs.pop("memory")
        elif self.config.memory is not None: self.memory = MemoryTable.from_file(self.config.memory)
        else:

            table_paths = []
            for simulation in self.all_simulations: table_paths.append(simulation.analysis.memory_table_path)
            if not sequences.all_none(table_paths):
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

        # Check passed arguments
        #has_assignment = kwargs.get("assignment", None) is not None
        #has_simulations = kwargs.get("simuations", None) is not None
        #has_simulation_names = kwargs.get("simulation_names", None) is not None
        #has_simulation_ids = kwargs.get("simulation_ids", None) is not None
        #if has_simulations and has_simulation_names: raise ValueError("Cannot specify both simulations and simulation names")
        # Set the assignment scheme if present
        #if has_assignment: self.assignment = kwargs.pop("assignment")

        # Get assignment table
        assignment = kwargs.pop("assignment", None)
        if self.config.assignment is not None:
            if assignment is not None: raise ValueError("Assignment is passed as an argument: cannot specify assignment table path")
            assignment = SimulationAssignmentTable.from_file(self.config.assignment)

        # Get simulations
        simulations = kwargs.pop("simulations", None)

        # Get simulation names
        simulation_names = kwargs.pop("simulation_names", None)
        if self.config.simulation_names is not None:
            if simulation_names is not None: raise ValueError("Simulation names are passed as an argument: cannot specify simulation names in configuration")
            simulation_names = self.config.simulation_names

        # Get simulation IDs
        simulation_ids = kwargs.pop("simulation_ids", None)
        if self.config.simulation_ids is not None:
            if simulation_ids is not None: raise ValueError("Simulation IDs are passed as an argument: cannot specify simulation IDs in configuration")
            simulation_ids = self.config.simulation_ids
        if simulation_ids is not None:
            # Check whether only one remote host is specified
            if self.config.remotes is None: raise ValueError("Remote host is not specified")
            nremotes = len(self.config.remotes)
            if nremotes != 1: raise ValueError("Only one remote host can be specified when giving simulation IDs")

        # CHECKS
        if simulations is not None and simulation_names is not None: raise ValueError("Cannot give simulations and simulation names")
        if simulations is not None and simulation_ids is not None: raise ValueError("Cannot give simulations and simulation IDs")
        if simulation_names is not None and simulation_ids is not None: raise ValueError("Cannot give simulation names and simulation IDs")

        # Set the assignment table
        if assignment is not None: self.assignment = assignment

        ## LOAD SIMULATIONS

        # Load simulations
        if simulations is not None: self.initialize_from_simulations(simulations)

        # Load simulations from names
        elif simulation_names is not None: self.initialize_from_simulation_names(simulation_names)

        # Load simulations from IDs
        elif simulation_ids is not None: self.initialize_from_simulation_ids(simulation_ids)

        # Load simulations from directory names
        elif self.config.from_directories: self.initialize_from_simulation_names(fs.directories_in_path(self.config.path, returns="name"))

        # Load simulations from assignment scheme, POSSIBLE?
        elif self.has_assignment: self.initialize_from_assignment()

        # NOT POSSIBLE
        else: raise ValueError("Not enough input to initialize simulations and/or assignment")

        # Set remotes now that assignment is present and correct: TEMPORARY HACK?
        self.config.remotes = self.host_ids

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

        # Debugging
        #log.debug("Trying to find the simulation '" + simulation_name + "' ...")

        the_host_id = None

        # Loop over the remotes
        for host_id in self.simulations_for_hosts:

            # Debugging
            #log.debug("Trying to find the simulation amongst the simulations of host '" + host_id + "' ...")

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
        
        #return self.find_host_id_for_simulation(simulation_name)

        # Loop over the simulations dictionary
        for host_id in self.simulations:            
            for name in self.simulations[host_id]:
                if name == simulation_name: return host_id
        
        # Not found
        return None

    # -----------------------------------------------------------------

    def get_simulation(self, simulation_name, host_id=None):

        """
        This function ...
        :param simulation_name:
        :param host_id:
        :return:
        """

        # Get the host ID if necessary:
        # if there is not a group under 'None' (host ID unknown) in the simulations dictionary
        if host_id is None and None not in self.simulations:
            host_id = self.host_id_for_simulation(simulation_name)
            if host_id is None: raise ValueError("Host ID for simulation '" + simulation_name + "' is not found")

        # Checks
        if host_id not in self.simulations: raise ValueError("No simulations for host '" + host_id + "'")
        if simulation_name not in self.simulations[host_id]: raise ValueError("Simulation '" + simulation_name + "' not in simulations for host '" + host_id + "'")

        # Get the simulation
        return self.simulations[host_id][simulation_name]

    # -----------------------------------------------------------------

    def add_simulation(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        self.simulations[simulation.host_id][simulation.name] = simulation

    # -----------------------------------------------------------------

    def is_remote_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        sim = self.get_simulation(simulation_name)
        return isinstance(sim, RemoteSimulation)

    # -----------------------------------------------------------------

    def get_simulations(self, names):
        
        """
        This function ...
        :param names: 
        :return: 
        """

        return [self.get_simulation(name) for name in names]

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

    def has_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        #return self.get_simulation(simulation_name).has_logfile
        return self.has_local_log(simulation_name) or self.has_remote_log(simulation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def get_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        #return self.get_simulation(simulation_name).log_file

        if self.has_local_log(simulation_name): return LogFile(self.get_local_log_path(simulation_name))
        elif self.has_remote_log(simulation_name): return LogFile.from_remote_file(self.get_remote_log_path(simulation_name), self.get_remote_for_simulation(simulation_name))
        else: return None

    # -----------------------------------------------------------------

    @memoize_method
    def get_logfiles(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # return self.get_simulation(simulation_name).logfiles()

        # LOCAL
        if self.has_local_log(simulation_name): return self.get_simulation(simulation_name).logfiles()

        # REMOTE
        elif self.has_remote_log(simulation_name):

            remote = self.get_remote_for_simulation(simulation_name)
            return [LogFile.from_remote_file(logfilepath, remote) for logfilepath in self.get_remote_log_paths(simulation_name)]

        # Nothing
        else: return []

    # -----------------------------------------------------------------

    def get_runtime(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_logfile(simulation_name):
            runtime = self.get_logfile(simulation_name).total_runtime
            if runtime is None: return None
            else: return runtime * u("s")
        else: return None

    # -----------------------------------------------------------------

    def get_peak_memory(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_logfile(simulation_name):
            memory = self.get_logfile(simulation_name).peak_memory
            if memory is None: return None
            else: return memory * u("GB")
        else: return None

    # -----------------------------------------------------------------

    def get_elapsed_time(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.has_logfile(simulation_name): return self.get_elapsed_time_logfile(simulation_name)
        else: return None

    # -----------------------------------------------------------------

    def get_elapsed_time_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Has local logfile
        if self.has_local_log(simulation_name):

            # Try the quick method
            try: return self.get_elapsed_time_local_logfile_quick(simulation_name)
            except ValueError: return self.get_logfile(simulation_name).elapsed_time * u("s")

        # Has remote logfile
        elif self.has_remote_log(simulation_name):

            # Try the quick method (doesn't read all of the log lines), but if it fails, use the normal method (reading the entire log file remotely)
            try: return self.get_elapsed_time_remote_logfile_quick(simulation_name)
            except ValueError: return self.get_logfile(simulation_name).elapsed_time * u("s")

        # Has no logfile
        else: raise IOError("Log file not present for simulation '" + simulation_name + "'")

    # -----------------------------------------------------------------

    def get_elapsed_time_local_logfile_quick(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the local log file path
        filepath = self.get_local_log_path(simulation_name)

        # Get first and last line
        first = fs.get_first_line(filepath)
        last = fs.get_last_line(filepath)

        # Get the elapsed time
        return self.get_elapsed_time_from_first_and_last_line(first, last)

    # -----------------------------------------------------------------

    def get_elapsed_time_remote_logfile_quick(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        filepath = self.get_remote_log_path(simulation_name)
        remote = self.get_remote_for_simulation(simulation_name)

        # Get first and last line
        first = remote.get_first_line(filepath)
        last = remote.get_last_line(filepath)

        # Get the elapsed time
        return self.get_elapsed_time_from_first_and_last_line(first, last)

    # -----------------------------------------------------------------

    def get_elapsed_time_from_first_and_last_line(self, first, last):

        """
        This function ...
        :param first:
        :param last:
        :return:
        """

        # Get the first time
        if not time.has_valid_timestamp(first): raise ValueError("Invalid timestamp")
        first_time = time.parse_line(first)

        # Get the last time
        if not time.has_valid_timestamp(last): raise ValueError("Invalid timestamp")
        last_time = time.parse_line(last)

        # Return the elapsed time
        seconds = (last_time - first_time).total_seconds()
        return seconds * u("s")

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

    def get_local_log_path(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Determine the path to the simulation log file
        local_log_file_path = simulation.log_file_path

        # Return the filepath        
        return local_log_file_path

    # -----------------------------------------------------------------

    def has_local_log(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        filepath = self.get_local_log_path(simulation_name)
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def get_remote_log_path(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return:
        """

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Get the remote
        #remote = self.get_remote(simulation.host_id)

        # The path to the simulation log file
        remote_log_file_path = simulation.remote_log_file_path

        # Check whether log file exists remotely?

        # Return the filepath
        return remote_log_file_path

    # -----------------------------------------------------------------

    def get_remote_log_paths(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Get paths
        remote = self.get_remote_for_simulation(simulation_name)
        return simulation.get_remote_log_file_paths(remote)

    # -----------------------------------------------------------------

    def has_remote_log(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        filepath = self.get_remote_log_path(simulation_name)
        return self.get_remote_for_simulation(simulation_name).is_file(filepath)

    # -----------------------------------------------------------------

    def get_log_lines(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # Simulation is retrieved
        if self.is_retrieved(simulation_name): lines = fs.read_lines(self.get_local_log_path(simulation_name))

        # Not yet retrieved, has remote log output already?
        elif self.has_remote_log(simulation_name): lines = self.get_remote_for_simulation(simulation_name).read_lines(self.get_remote_log_path(simulation_name))

        # No output yet
        else: lines = []

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    def get_execution_handle(self, simulation_name):
        return self.get_simulation(simulation_name).handle

    # -----------------------------------------------------------------

    def has_execution_handle(self, simulation_name):
        return self.is_remote_simulation(simulation_name) and self.get_execution_handle(simulation_name) is not None

    # -----------------------------------------------------------------

    def is_screen_execution(self, simulation_name):
        return self.get_execution_handle(simulation_name).is_screen

    # -----------------------------------------------------------------

    def is_job_execution(self, simulation_name):
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

    def get_queued_and_running_simulation_names_for_screen(self, screen_name):

        """
        This function ...
        :param screen_name:
        :return:
        """

        names = []
        for simulation in self.all_queued_or_running_simulations:
            if self.get_screen_name(simulation.name) != screen_name: continue
            #yield simulation.name
            names.append(simulation.name)
        return names

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

    def get_screen_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the screenlog path
        filepath = self.get_screen_log_path(simulation_name)
        remote_ski_path = self.get_remote_skifile_path(simulation_name)

        # Show the remote logfile path
        log.debug("Remote screen log path: '" + filepath + "'")

        # Get the lines
        lines = self.get_remote_for_simulation(simulation_name).get_lines_between(filepath, remote_ski_path, "Finished simulation")

        # Return the lines
        return lines

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
        handle = self.get_execution_handle(simulation_name)
        return handle.remote_job_script_path if hasattr(handle, "remote_job_script_path") else None

    # -----------------------------------------------------------------

    def has_remote_job_script(self, simulation_name):
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
        return self.get_local_job_script_path(simulation_name) is not None

    # -----------------------------------------------------------------

    def find_remote_job_script_path(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get base path
        base_path = self.get_remote_base_path(simulation_name)

        # Look for sh files
        filepaths = self.get_remote_for_simulation(simulation_name).files_in_path(base_path, extension="sh")

        # Return unique
        return sequences.get_single(filepaths, method="none") # return None if there are multiple: we can't be sure of which one is the actual latest?

    # -----------------------------------------------------------------

    @memoize_method
    def get_job_script(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get host ID and cluster name
        host_id = self.get_host_id_for_simulation(simulation_name)
        #cluster_name = self.get_cluster_name_for_simulation(simulation_name)
        cluster = self.get_cluster_for_simulation(simulation_name)
        remote = self.get_remote(host_id)

        # Local job script is found
        if self.has_local_job_script(simulation_name):

            filepath = self.get_local_job_script_path(simulation_name)
            return SKIRTJobScript.from_file(filepath, host_id=host_id, cluster=cluster, remote=remote)

        # Remote job script is found
        elif self.has_remote_job_script(simulation_name):

            filepath = self.get_remote_job_script_path(simulation_name)
            return SKIRTJobScript.from_remote_file(filepath, remote, host_id=host_id, cluster=cluster)

        # No job script is found
        else: #return None

            filepath = self.find_remote_job_script_path(simulation_name)
            if filepath is None: return None
            else: return SKIRTJobScript.from_remote_file(filepath, remote, host_id=host_id, cluster=cluster)

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
        return self.get_simulation(simulation_name).extraction_path

    # -----------------------------------------------------------------

    def get_simulation_plotting_path(self, simulation_name):
        return self.get_simulation(simulation_name).plotting_path

    # -----------------------------------------------------------------

    def get_simulation_misc_path(self, simulation_name):
        return self.get_simulation(simulation_name).misc_path

    # -----------------------------------------------------------------

    @memoize_method
    def get_input(self, simulation_name):
        return self.get_simulation(simulation_name).input

    # -----------------------------------------------------------------

    @memoize_method
    def get_remote_input(self, simulation_name):
        simulation = self.get_simulation(simulation_name)
        return simulation.get_remote_input(remote=self.get_remote(simulation.host_id))

    # -----------------------------------------------------------------

    def get_base_path(self, simulation_name):
        return self.get_simulation(simulation_name).base_path

    # -----------------------------------------------------------------

    def get_remote_base_path(self, simulation_name):
        return self.get_simulation(simulation_name).remote_simulation_path

    # -----------------------------------------------------------------

    def get_output_path(self, simulation_name):
        return self.get_simulation(simulation_name).output_path

    # -----------------------------------------------------------------

    def get_remote_output_path(self, simulation_name):
        return self.get_simulation(simulation_name).remote_output_path

    # -----------------------------------------------------------------

    def get_extraction_path(self, simulation_name):
        return self.get_simulation(simulation_name).extraction_path

    # -----------------------------------------------------------------

    def has_extraction_directory(self, simulation_name):
        path = self.get_extraction_path(simulation_name)
        return path is not None and fs.is_directory(path)

    # -----------------------------------------------------------------

    def get_plotting_path(self, simulation_name):
        return self.get_simulation(simulation_name).plotting_path

    # -----------------------------------------------------------------

    def has_plotting_directory(self, simulation_name):
        path = self.get_plotting_path(simulation_name)
        return path is not None and fs.is_directory(path)

    # -----------------------------------------------------------------

    def get_misc_path(self, simulation_name):
        return self.get_simulation(simulation_name).misc_path

    # -----------------------------------------------------------------

    def has_misc_directory(self, simulation_name):
        path = self.get_misc_path(simulation_name)
        return path is not None and fs.is_directory(path)

    # -----------------------------------------------------------------

    @memoize_method_reset
    def get_output(self, simulation_name):
        return self.get_simulation(simulation_name).output

    # -----------------------------------------------------------------

    def reset_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        self.get_output._reset_for_args(simulation_name)

    # -----------------------------------------------------------------

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
    # EXTRACTION
    # -----------------------------------------------------------------

    @memoize_method_reset
    def get_extraction_output(self, simulation_name):
        return self.get_simulation(simulation_name).extraction_output

    # -----------------------------------------------------------------

    def reset_extraction_output(self, simulation_name):
        self.get_extraction_output._reset_for_args(simulation_name)

    # -----------------------------------------------------------------

    def get_nextraction_output(self, simulation_name):
        return len(self.get_extraction_output(simulation_name))

    # -----------------------------------------------------------------

    def has_extraction_output(self, simulation_name):
        if not self.has_extraction_directory(simulation_name): return False
        return self.get_nextraction_output(simulation_name) > 0

    # -----------------------------------------------------------------

    def get_nextraction_files(self, simulation_name, otypes):
        return self.get_extraction_output(simulation_name).get_nfiles_types(otypes)

    # -----------------------------------------------------------------

    def get_nextraction_directories(self, simulation_name, otypes):
        return self.get_extraction_output(simulation_name).get_ndirectories_types(otypes)

    # -----------------------------------------------------------------

    def get_nextraction_files_and_directories(self, simulation_name, otypes):
        return self.get_nextraction_files(simulation_name, otypes) + self.get_nextraction_directories(simulation_name, otypes)

    # -----------------------------------------------------------------

    def is_cached_extraction(self, simulation_name):
        return self.get_extraction_output(simulation_name).has_cached

    # -----------------------------------------------------------------
    # ANALYSIS
    # -----------------------------------------------------------------

    def has_analysis_output(self, simulation_name):
        return self.has_extraction_output(simulation_name) or self.has_plotting_output(simulation_name) or self.has_misc_output(simulation_name)

    # -----------------------------------------------------------------

    def has_timeline_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if not self.has_extraction_output(simulation_name): return False
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
        else: return extract_timeline(log_files=self.get_logfiles(simulation_name))

        #return extract_timeline(self.get_simulation(simulation_name))

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

    @memoize_method_reset
    def get_plotting_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).plotting_output

    # -----------------------------------------------------------------

    def reset_plotting_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        self.get_plotting_output._reset_for_args(simulation_name)

    # -----------------------------------------------------------------

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

    @memoize_method_reset
    def get_misc_output(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).misc_output

    # -----------------------------------------------------------------

    def reset_misc_output(self, simulation_name):

        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        self.get_misc_output._reset_for_args(simulation_name)
        
    # -----------------------------------------------------------------

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

    def get_cluster_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_host_for_simulation(simulation_name).cluster

    # -----------------------------------------------------------------

    @memoize_method
    def get_parallelization_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Check if the parallelization is defined
        if simulation.parallelization is not None: return simulation.parallelization
        elif self.has_execution_handle(simulation_name):
            if self.is_screen_execution(simulation_name): return self.get_parallelization_for_simulation_screen(simulation_name)
            elif self.is_job_execution(simulation_name): return self.get_parallelization_for_simulation_job(simulation_name)
            else: raise NotImplementedError("Execution handle not supported")
        elif self.has_logfile(simulation_name): return self.get_parallelization_for_simulation_logfile(simulation_name)
        else: raise ValueError("Parallelization cannot be determined for simulation '" + simulation_name + "'")

    # -----------------------------------------------------------------

    def get_parallelization_for_simulation_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_logfile(simulation_name).parallelization

    # -----------------------------------------------------------------

    def get_ncores_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_parallelization_for_simulation(simulation_name).ncores

    # -----------------------------------------------------------------

    def get_parallelization_for_simulation_screen(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        screen_script = self.get_screen_script(simulation_name)
        if screen_script is None: raise RuntimeError("Cannot determine parallelization scheme: original screen script not found")
        parallelization = screen_script.get_parallelization(simulation_name)
        return parallelization

    # -----------------------------------------------------------------

    def get_parallelization_for_simulation_job(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        job_script = self.get_job_script(simulation_name)
        if job_script is None: raise RuntimeError("Cannot determine parallelization scheme: original job script not found")
        parallelization = job_script.parallelization
        return parallelization

    # -----------------------------------------------------------------

    def get_logging_options_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Is an execution handle defined for the simulation?
        if self.has_execution_handle(simulation_name):
    
            # Screen
            if self.is_screen_execution(simulation_name): options = self.get_logging_options_for_simulation_screen(simulation_name, return_none=True)

            # Job
            elif self.is_job_execution(simulation_name): options = self.get_logging_options_for_simulation_job(simulation_name, return_none=True)

            # Not supported
            else: raise NotImplementedError("Execution handle not supported")

        # Not defined
        else: options = None

        # None?
        if options is None:

            if self.has_logfile(simulation_name): options = self.get_logging_options_for_simulation_logfile(simulation_name)
            else: raise RuntimeError("Cannot determine logging options: original screen or job script not found, and remote or local log file not yet created")

        # Return
        return options

    # -----------------------------------------------------------------
        
    def get_logging_options_for_simulation_screen(self, simulation_name, return_none=False):

        """
        This function ...
        :param simulation_name:
        :param return_none:
        :return:
        """

        screen_script = self.get_screen_script(simulation_name)

        if screen_script is None:
            if return_none: return None
            else: raise RuntimeError("Cannot determine logging options: original screen script not found")

        logging_options = screen_script.get_logging_options(simulation_name)
        return logging_options

    # -----------------------------------------------------------------

    def get_logging_options_for_simulation_job(self, simulation_name, return_none=False):
        
        """
        This function ...
        :param simulation_name:
        :param return_none:
        :return: 
        """

        job_script = self.get_job_script(simulation_name)

        if job_script is None:
            if return_none: return None
            else: raise RuntimeError("Cannot determine logging options: original job script not found")

        logging_options = job_script.logging_options
        return logging_options

    # -----------------------------------------------------------------

    def get_logging_options_for_simulation_logfile(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the logfile
        logfile = self.get_logfile(simulation_name)

        # Return the logging options
        return logfile.logging_options

    # -----------------------------------------------------------------

    def get_analysis_options_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysis

    # -----------------------------------------------------------------

    def get_scheduling_options_for_simulation(self, simulation_name):

        """
        THis function ...
        :param simulation_name:
        :return:
        """

        # Screen: none
        if self.is_screen_execution(simulation_name): return None

        # Job
        elif self.is_job_execution(simulation_name): return self.get_scheduling_options_for_simulation_job(simulation_name)

        # Not supported
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    def get_scheduling_options_for_simulation_job(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        job_script = self.get_job_script(simulation_name)
        if job_script is None: raise RuntimeError("Cannot determine scheduling options: original job script not found")
        scheduling_options = job_script.scheduling_options
        return scheduling_options

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

    def initialize_from_assignment(self):

        """
        This function ...
        :return:
        """

        # Set remotes?
        if self.config.remotes is None: self.config.remotes = self.assignment.unique_host_ids

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

        # Initialize assignment table?
        if not self.has_assignment:
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

            # Remote simulation?
            if isinstance(simulation, RemoteSimulation):
                simulation_id = simulation.id
                host_id = simulation.host_id
                cluster_name = simulation.cluster_name

            # Regular simulation object
            elif isinstance(simulation, SkirtSimulation): simulation_id = host_id = cluster_name = None

            # Invalid
            else: raise ValueError("Invalid type for simulation '" + simulation_name + "'")

            # Add to assignment
            if self._new_assignment: self.assignment.add_simulation(simulation_name, host_id=host_id, cluster_name=cluster_name,
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

    def initialize_from_simulation_ids(self, ids):

        """
        This function ...
        :param ids:
        :return:
        """

        # Get the host ID
        host_id = self.config.remotes[0]

        # Initialize assignment table
        self.assignment = SimulationAssignmentTable()
        self._new_assignment = True

        # Get the simulation IDs for the remote host
        simulation_ids = introspection.simulation_ids_for_host(host_id)

        # Loop over the simulation IDs
        for simulation_id in ids:

            # Check whether the simulation exists
            if simulation_id not in simulation_ids: raise ValueError("There is no simulation with ID '" + str(simulation_id) + "' for remote host '" + host_id + "'")

            # Load the simulation
            simulation = get_simulation_for_host(host_id, simulation_id)

            # Add to assignment
            self.assignment.add_simulation_object(simulation, success=self.config.success)

            # Add the simulation
            self.simulations[host_id][simulation.name] = simulation

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

    def get_simulation_status(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Getting status of simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Analysed
        if simulation.analysed: simulation_status = analysed_name

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
        elif isinstance(simulation, RemoteSimulation) and simulation.retrieved: simulation_status = retrieved_name

        # Finished
        elif isinstance(simulation, RemoteSimulation) and simulation.finished: simulation_status = finished_name

        # Remote simulation that is not yet retrieved
        elif isinstance(simulation, RemoteSimulation):

            host_id = simulation.host_id
            screen_states = self.screens[host_id]
            jobs_status = self.jobs[host_id]
            with log.no_debugging(): simulation_status = self.get_remote(host_id).get_simulation_status(simulation, screen_states=screen_states, jobs_status=jobs_status)

        else: simulation_status = "unknown"

        # Check success flag in assignment
        if self.config.fix_success and not self.assignment.is_launched(simulation.name) and not is_invalid_or_unknown_status(simulation_status):
            log.warning("Settting the launch of simulation '" + simulation.name + "' as succesful in the assignment table as this was not yet done")
            self.assignment.set_success_for_simulation(simulation.name)

        # Retrieve finished simulations?
        if simulation_status == finished_name and self.config.retrieve:
            self.retrieve_simulation(simulation_name)
            simulation_status = "retrieved"

        # Return the status
        return simulation_status

    # -----------------------------------------------------------------

    @lazyproperty
    def status(self):

        """
        This function ...
        :return:
        """

        # Initialize lists
        status_list = []

        # Debugging
        log.debug("Creating simulation status table ...")

        # TODO: show progress bar!

        # Loop over the simulations
        for index, simulation_name in enumerate(self.simulation_names):

            # Get the simulation status
            simulation_status = self.get_simulation_status(simulation_name)

            # Add the status
            status_list.append(simulation_status)

        # Create the table and return
        return SimulationStatusTable.from_columns(self.simulation_names, status_list)

    # -----------------------------------------------------------------

    def reset_status_for_simulation(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # If status is not created yet, return
        if not self.has_status: return

        # Get the new status
        simulation_status = self.get_simulation_status(simulation_name)

        # Set the new status
        self.status.reset_for_simulation(simulation_name, simulation_status)

    # -----------------------------------------------------------------

    def retrieve_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the simulation name
        simulation_names = self.get_simulation_names_from_command(command, **kwargs)

        # Retrieve all simulations
        for simulation_name in simulation_names:

            # Check
            if not self.is_finished(simulation_name): raise ValueError("Simulation '" + simulation_name + "' is not finished")

            # Retrieve the simulation
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

    def retrieve_all(self, **kwargs):
        
        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Loop over all finished but not yet retrieved simulations
        for simulation_name in self.all_finished_not_retrieved_simulation_names:

            # Retrieve the simulation
            self.retrieve_simulation(simulation_name)

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
        return "status" in self.__dict__

    # -----------------------------------------------------------------

    @property
    def screens(self):
        return self.remotes.screens

    # -----------------------------------------------------------------

    @property
    def jobs(self):
        return self.remotes.jobs

    # -----------------------------------------------------------------

    def get_status(self, simulation_name):
        return self.status.get_status(simulation_name)

    # -----------------------------------------------------------------

    def is_queued(self, simulation_name):
        return self.status.is_queued(simulation_name)

    # -----------------------------------------------------------------

    def is_running(self, simulation_name):
        return self.status.is_running(simulation_name)

    # -----------------------------------------------------------------

    def is_queued_or_running(self, simulation_name):
        return self.status.is_queued_or_running(simulation_name)

    # -----------------------------------------------------------------

    def is_finished(self, simulation_name):
        return self.status.is_finished(simulation_name)

    # -----------------------------------------------------------------

    def is_running_or_finished(self, simulation_name):
        return self.status.is_running_or_finished(simulation_name)

    # -----------------------------------------------------------------

    def is_retrieved(self, simulation_name):
        return self.status.is_retrieved(simulation_name)

    # -----------------------------------------------------------------

    def is_analysed(self, simulation_name):
        return self.status.is_analysed(simulation_name)

    # -----------------------------------------------------------------

    def is_cancelled(self, simulation_name):
        return self.status.is_cancelled(simulation_name)

    # -----------------------------------------------------------------

    def is_aborted(self, simulation_name):
        return self.status.is_aborted(simulation_name)

    # -----------------------------------------------------------------

    def is_crashed(self, simulation_name):
        return self.status.is_crashed(simulation_name)

    # -----------------------------------------------------------------

    def is_running_or_finished_or_aborted_or_crashed(self, simulation_name):
        return self.status.is_running_or_finished_or_aborted_or_crashed(simulation_name)

    # -----------------------------------------------------------------

    def is_failed(self, simulation_name):
        return self.status.is_failed(simulation_name)

    # -----------------------------------------------------------------

    def has_unknown_status(self, simulation_name):
        return self.status.has_unknown_status(simulation_name)

    # -----------------------------------------------------------------

    @property
    def nfinished(self):
        return self.status.nfinished

    # -----------------------------------------------------------------

    @property
    def relative_nfinished(self):
        return self.status.relative_nfinished

    # -----------------------------------------------------------------

    @property
    def percentage_nfinished(self):
        return self.status.percentage_nfinished

    # -----------------------------------------------------------------

    @property
    def nretrieved(self):
        return self.status.nretrieved

    # -----------------------------------------------------------------

    @property
    def relative_nretrieved(self):
        return self.status.relative_nretrieved

    # -----------------------------------------------------------------

    @property
    def percentage_nretrieved(self):
        return self.status.percentage_nretrieved

    # -----------------------------------------------------------------

    @property
    def nanalysed(self):
        return self.status.nanalysed

    # -----------------------------------------------------------------

    @property
    def relative_nanalysed(self):
        return self.status.relative_nanalysed

    # -----------------------------------------------------------------

    @property
    def percentage_nanalysed(self):
        return self.status.percentage_nanalysed

    # -----------------------------------------------------------------

    @property
    def nrunning(self):
        return self.status.nrunning

    # -----------------------------------------------------------------

    @property
    def relative_nrunning(self):
        return self.status.relative_nrunning

    # -----------------------------------------------------------------

    @property
    def percentage_nrunning(self):
        return self.status.percentage_nrunning

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

        # Set the SEDs to plot
        seds = OrderedDict()
        seds.update(self.reference_seds)

        # Plot the SEDs
        for path in output.seds:

            # Get instrument name
            instr_name = get_sed_instrument_name(path, prefix)

            # Load SEDs
            sed_contributions = load_multiple_seds(path, as_dict=True)

            # Loop over the SEDs
            for contribution in sed_contributions:

                # Determine label
                label = instr_name + " " + contribution

                # Get the SED
                sed = seds[contribution]

                # Add the SED to the dictionary
                seds[label] = sed

        # Make (and show) the plot
        plot_seds(seds, path=path, show_file=True)

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
        definition.add_flag("share_normalization", "share normalization between the different frames", False)

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
        self.plot_simulation_datacubes(simulation_name, contributions=config.contributions, instruments=config.instruments, share_normalization=config.share_normalization)

    # -----------------------------------------------------------------

    def plot_simulation_datacubes(self, simulation_name, contributions=None, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param contributions:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting simulated datacubes for simulation '" + simulation_name + "' ...")

        # Total
        if contributions is None or "total" in contributions: self.plot_total_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

        # Direct
        if contributions is None or "direct" in contributions: self.plot_direct_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

        # Transparent
        if contributions is None or "transparent" in contributions: self.plot_transparent_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

        # Scattered
        if contributions is None or "scattered" in contributions: self.plot_scattered_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

        # Dust
        if contributions is None or "dust" in contributions: self.plot_dust_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

        # Dust-scattered
        if contributions is None or "dustscattered" in contributions: self.plot_dustscattered_datacubes(simulation_name, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def _plot_datacubes(self, simulation_name, paths, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param paths:
        :param instruments:
        :param share_normalization:
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
            plotting.plot_datacube(datacube, title=instr_name, share_normalization=share_normalization, show_axes=False)

    # -----------------------------------------------------------------

    def plot_total_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Showing total datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.total_images, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def plot_direct_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting direct datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.direct_images, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def plot_transparent_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting transparent datacubes for simulation '" + simulation_name + "' ...")

        # Get simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.transparent_images, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def plot_scattered_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting scattered datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.scattered_images, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def plot_dust_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting dust datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Plot the datacubes
        self._plot_datacubes(simulation_name, output.dust_images, instruments=instruments, share_normalization=share_normalization)

    # -----------------------------------------------------------------

    def plot_dustscattered_datacubes(self, simulation_name, instruments=None, share_normalization=False):

        """
        This function ...
        :param simulation_name:
        :param instruments:
        :param share_normalization:
        :return:
        """

        # Debugging
        log.debug("Plotting dust-scattered datacubes for simulation '" + simulation_name + "' ...")

        # Get the simulation output
        output = self.get_output(simulation_name)

        # Show
        self._plot_datacubes(simulation_name, output.dust_scattered_images, instruments=instruments, share_normalization=share_normalization)

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

        # Set SEDs
        seds = OrderedDict()
        seds.update(self.reference_seds)
        seds["Mock observation"] = fluxes

        # Make (and show) the plot
        plot_seds(seds, path=path, show_file=True)

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

        # Set SEDs
        seds = OrderedDict()
        seds.update(self.reference_seds)
        seds["Mock observation"] = fluxes

        # Make (and show) the plot
        plot_seds(seds, path=path, show_file=True)

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
        print("[LOCAL DIRECTORY: " + self.get_output_path(simulation_name) + "]")
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
        print("[REMOTE HOST: " + self.get_host_id_for_simulation(simulation_name) + "]")
        print("[DIRECTORY: " + self.get_remote_output_path(simulation_name) + "]")
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
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get the host
        host = config.pop("host")

        # Return
        return splitted, host, config

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
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Return
        hosts = config.pop("hosts")
        return splitted, hosts, config

    # -----------------------------------------------------------------

    def get_hosts_from_command(self, command, name=None, index=1, required=False, choices=None):

        """
        This function ...
        :param command:
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

    def get_host_parallelization_command_definition(self, command_definition=None, required_to_optional=True):

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

    def parse_host_parallelization_command(self, command, command_definition=None, name=None, index=1,
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
        definition = self.get_host_parallelization_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get host and parallelization
        host = config.pop("host")
        parallelization = config.pop("parallelization")

        # Return
        return splitted, host, parallelization, config

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
        splitted, host, parallelization, config = self.parse_host_parallelization_command(command, name=name, index=index, interactive=interactive)

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
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation name
        if types.is_integer_type(config.simulation): simulation_name = self.simulation_names[config.pop("simulation")]
        else:
            simulation_name = config.pop("simulation")
            if simulation_name not in self.simulation_names: raise InvalidCommandError("Invalid simulation name: '" + simulation_name + "'", command)

        # Return
        return splitted, simulation_name, config

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
        definition.add_required("simulations", "integer_and_string_list", "simulation indices or names")
        #definition.add_positional_optional("simulations", "integer_and_string_list", "simulation indices or names", default=["all"])

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
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation names
        nsimulations = len(config.simulations)
        if nsimulations == 1 and config.simulations[0] == "all": simulation_names = self.simulation_names
        else:
            simulation_names = []
            for index_or_name in config.simulations:
                if types.is_string_type(index_or_name):
                    simulation_name = index_or_name
                    if simulation_name not in self.simulation_names: raise InvalidCommandError("Invalid simulation name '" + simulation_name + "'", command)
                elif types.is_integer_type(index_or_name): simulation_name = self.simulation_names[index_or_name]
                else: raise ValueError("Invalid type")
                simulation_names.append(simulation_name)

        # Return
        return splitted, simulation_names, config

    # -----------------------------------------------------------------

    def get_simulations_simulation_command_definition(self, command_definition=None, required_to_optional=True):

        """
        This function ...
        :param command_definition:
        :param required_to_optional:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_required("simulations", "integer_list", "simulation indices")
        definition.add_required("simulation", "integer_or_string", "simulation index or name")

        # Add definition settings
        if command_definition is not None:
            if required_to_optional: definition.import_settings(command_definition, required_to="optional")
            else: definition.import_settings(command_definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def parse_simulations_simulation_command(self, command, command_definition=None, name=None, index=1, required_to_optional=True, interactive=False):
        
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

        # Get the definition
        definition = self.get_simulations_simulation_command_definition(command_definition, required_to_optional=required_to_optional)

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=splitted[index:], error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation names
        simulation_names = [self.simulation_names[index] for index in config.pop("simulations")]

        # Get simulation name
        if types.is_integer_type(config.simulation): simulation_name = self.simulation_names[config.pop("simulation")]
        else:
            simulation_name = config.pop("simulation")
            if simulation_name not in self.simulation_names: raise InvalidCommandError("Invalid simulation name '" + simulation_name + "'", command)

        # Return
        return splitted, simulation_names, simulation_name, config

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
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=splitted[index:], error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Get simulation_a name
        if types.is_integer_type(config.simulation_a): simulation_a_name = self.simulation_names[config.pop("simulation_a")]
        else:
            simulation_a_name = config.pop("simulation_a")
            if simulation_a_name not in self.simulation_names: raise InvalidCommandError("Invalid simulation name '" + simulation_a_name + "'", command)

        # Get simulation_b name
        if types.is_integer_type(config.simulation_b): simulation_b_name = self.simulation_names[config.pop("simulation_b")]
        else:
            simulation_b_name = config.pop("simulation_b")
            if simulation_b_name not in self.simulation_names: raise InvalidCommandError("Invalid simulation name '" + simulation_b_name + "'", command)

        # Return
        return splitted, simulation_a_name, simulation_b_name, config

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

    def get_simulation_names_and_config_from_command(self, command, command_definition, name=None, interactive=False):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :return:
        """

        # Parse the command
        splitted, simulation_names, config = self.parse_simulations_command(command, command_definition=command_definition, name=name, interactive=interactive)

        # Return the simulation names and config
        return simulation_names, config

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

    def get_simulation_name_and_config_from_command(self, command, command_definition, name=None, interactive=False, required_to_optional=True):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :param interactive:
        :param required_to_optional:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, command_definition=command_definition, name=name, interactive=interactive, required_to_optional=required_to_optional)

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

        # Get the log lines
        lines = self.get_log_lines(simulation_name)

        # Show the lines
        self.show_simulation_output_lines(lines, summarize=summarize, short=short)

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

        # Get the screen output lines for the simulation
        lines = self.get_screen_output(simulation_name)
        nlines = len(lines)
        if nlines == 0:
            log.warning("No screen output for simulation '" + simulation_name + "'")
            return

        # Show the lines
        self.show_simulation_output_lines(lines, summarize=summarize, short=short)

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
        self.show_simulation_output_lines(lines, summarize=summarize, short=short)

    # -----------------------------------------------------------------

    def show_simulation_output_lines(self, lines, summarize=False, short=False):

        """
        This function ...
        :param lines:
        :param summarize:
        :param short:
        :return:
        """

        # Show
        print("")
        if summarize: show_log_summary(lines, debug_output=True)
        elif short: show_log_summary(lines, debug_output=False)
        else:
            for line in lines:

                # Show the line
                if not time.has_valid_timestamp(line): print(fmt.red + line + fmt.reset)
                else: print(line)

        print("")

    # -----------------------------------------------------------------

    def get_simulation_error_lines(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Executed in screen
        if self.is_screen_execution(simulation_name): return self.get_simulation_error_lines_screen(simulation_name)

        # Executed in job
        elif self.is_job_execution(simulation_name): return self.get_simulation_error_lines_job(simulation_name)

        # Invalid
        else: raise NotImplementedError("Execution handle not supported")

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_error_lines_screen(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get the screen output lines for the simulation
        all_lines = self.get_screen_output(simulation_name)

        # Initialize for error lines
        lines = []

        # Loop over the lines
        break_after = False
        for line in all_lines:

            # Skip regular SKIRT output messages
            if time.has_valid_timestamp(line): continue

            # Add the line
            lines.append(line)
            if break_after: break

            # Skip after certain MPI message
            if "YOU CAN IGNORE THE BELOW CLEANUP MESSAGES" in line: break_after = True

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_error_lines_job(self, simulation_name):

        """
        Thisf unction ...
        :param simulation_name:
        :return:
        """

        # Get the error lines
        all_lines = self.get_job_error(simulation_name)

        # Initialize for error lines
        lines = []

        # Loop over the lines
        break_after = False
        for line in all_lines:

            # Skip regular SKIRT output messages
            if time.has_valid_timestamp(line): continue

            # Add the line
            lines.append(line)
            if break_after: break

            # Skip after certain MPI message
            if "YOU CAN IGNORE THE BELOW CLEANUP MESSAGES" in line: break_after = True

        # Return the lines
        return lines

    # -----------------------------------------------------------------

    def has_exceed_walltime(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Loop over the error lines
        for line in self.get_simulation_error_lines(simulation_name):

            # Check for walltime exceeded message
            if "killed: walltime" in line and "exceeded limit" in line: return True

        # Error message not encountered
        return False

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

        # Get the screen output lines for the simulation
        lines = self.get_screen_output(simulation_name)

        # Show the lines
        self.show_simulation_error_lines(lines)

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

        # Show lines
        self.show_simulation_error_lines(lines)

    # -----------------------------------------------------------------

    def show_simulation_error_lines(self, lines):
        
        """
        This function ...
        :param lines: 
        :return: 
        """

        # Check number of lines
        nlines = len(lines)
        if nlines == 0:
            print(fmt.red + "no error output" + fmt.reset)
            return

        # Show
        print("")
        break_after = False
        for line in lines:

            # Skip regular SKIRT output messages
            if time.has_valid_timestamp(line): continue

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

        definition = show_simulation_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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

        definition = show_analysis_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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

    def show_analysis_steps_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show analysis steps
        self.show_analysis_steps(simulation_name)

    # -----------------------------------------------------------------

    def show_analysis_steps(self, simulation_name):
        
        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # Debugging
        log.debug("Showing analysis steps for simulation '" + simulation_name + "' ...")
        
        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Show steps that will be performed
        show_analysis_steps(simulation)
        
    # -----------------------------------------------------------------

    @lazyproperty
    def adapt_simulation_settings_definition(self):

        """
        This function ...
        :return:
        """

        definition = adapt_simulations_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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
        simulation_names, config = self.get_simulation_names_and_config_from_command(command, self.adapt_simulation_settings_definition, name=name, **kwargs)
        nsimulations = len(simulation_names)

        # Adapt single simulation
        if nsimulations == 1: self.adapt_simulation_settings(simulation_names[0], config=config)

        # Adapt multiple simulations
        else: self.adapt_simulations_settings(simulation_names, config=config)

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

    def adapt_simulations_settings(self, simulation_names, config=None):

        """
        This function ...
        :param simulation_names:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Adapting simulations '" + ",".join(simulation_names) + "' ...")

        # Get simulations
        simulations = self.get_simulations(simulation_names)

        # Adapt
        adapt_simulations(*simulations, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def adapt_analysis_options_definition(self):

        """
        This function ...
        :return:
        """

        definition = adapt_analysis_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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
        simulation_names, config = self.get_simulation_names_and_config_from_command(command, self.adapt_analysis_options_definition, name=name, **kwargs)
        nsimulations = len(simulation_names)

        # Adapt single simulation
        if nsimulations == 1: self.adapt_analysis_options(simulation_names[0], config=config)

        # Adapt multiple simulations
        else: self.adapt_analysis_options_simulations(simulation_names, config=config)

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

    def adapt_analysis_options_simulations(self, simulation_names, config=None):

        """
        This function ...
        :param simulation_names:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Adapt analysis options of simulations '" + ", ".join(simulation_names) + "' ...")

        # Get the simulations
        simulations = self.get_simulations(simulation_names)

        # Adapt
        adapt_analysis_simulations(*simulations, config=config)

    # -----------------------------------------------------------------

    @lazyproperty
    def compare_simulation_settings_definition(self):

        """
        This function ...
        :return:
        """

        definition = show_simulation_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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
        #splitted, simulation_a_name, simulation_b_name, config = self.parse_two_simulations_command(command, self.compare_simulation_settings_definition, name=name, **kwargs)
        # Get the simulation names
        splitted, simulation_names, config = self.parse_simulations_command(command, self.compare_simulation_settings_definition, name=name, **kwargs)

        # Compare
        #self.compare_simulation_settings(simulation_a_name, simulation_b_name, config=config)
        self.compare_simulation_settings(simulation_names, config=config)

    # -----------------------------------------------------------------

    def compare_two_simulation_settings(self, simulation_a_name, simulation_b_name, config=None):

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

    def compare_simulation_settings(self, simulation_names, config=None):

        """
        This function ....
        :param simulation_names:
        :param config:
        :return:
        """

        # Debugging
        nsimulations = len(simulation_names)
        log.debug("Comparing " + str(nsimulations) + "' simulations ...")

        # Get the simulations
        simulations = self.get_simulations(simulation_names)

        # Compare
        compare_simulations(*simulations, config=config)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def compare_analysis_options_definition(self):

        """
        This function ...
        :return:
        """

        definition = show_analysis_definition.copy(pos_optional=False)
        definition.remove_optional("names")
        definition.remove_flag("from_directories")
        return definition

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
        #splitted, simulation_a_name, simulation_b_name, config = self.parse_two_simulations_command(command, self.compare_analysis_options_definition, name=name, **kwargs)
        splitted, simulation_names, config = self.parse_simulations_command(command, self.compare_analysis_options_definition, name=name, **kwargs)

        # Compare
        #self.compare_analysis_options(simulation_a_name, simulation_b_name, config=config)
        self.compare_analysis_options(simulation_names, config=config)

    # -----------------------------------------------------------------

    def compare_two_analysis_options(self, simulation_a_name, simulation_b_name, config=None):

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

    def compare_analysis_options(self, simulation_names, config=None):
        
        """
        This function ...
        :param simulation_names: 
        :param config: 
        :return: 
        """

        # Debugging
        nsimulations = len(simulation_names)
        log.debug("Comparing analysis options between " + str(nsimulations) + " simulations ...")

        # Get the simulations
        simulations = self.get_simulations(simulation_names)

        # Compare
        compare_analysis(*simulations, config=config)
        print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def steal_simulation_settings_definition(self):
        
        """
        This function ...
        :return: 
        """
        
        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Select certain settings
        definition.add_positional_optional("matching", "string", "only select settings with a name matching this string")
        definition.add_optional("contains", "string", "only select settings containing this string in their name")
        definition.add_optional("not_contains", "string", "don't select settings containing this string in their name")
        definition.add_optional("exact_name", "string", "only select settings with this exact string as their name")
        definition.add_optional("exact_not_name", "string", "don't select settings with this exact string as their name")
        definition.add_optional("startswith", "string", "only select settings whose name starts with this string")
        definition.add_optional("endswith", "string", "only select settings whose name starts with this string")

        # Flags
        definition.add_flag("successive", "select settings successively", False)
        definition.add_flag("confirm", "confirm before adapting", True)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def steal_simulation_settings_command(self, command, **kwargs):
        
        """
        This function ...
        :param command: 
        :param kwargs: 
        :return: 
        """

        # Set name
        name = _steal_command_name + " " + _simulation_command_name

        # Get simulation names and config
        splitted, simulation_names, from_simulation_name, config = self.parse_simulations_simulation_command(command, self.steal_simulation_settings_definition, name=name, **kwargs)

        # Check settings
        if config.matching is not None:
            if config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
            config.contains = config.matching

        # Get the settings here (so they don't have to be selected for each stealing simulation)
        settings = select_simulation_settings(successive=config.successive, contains=config.contains, not_contains=config.not_contains,
                                              exact_name=config.exact_name, exact_not_name=config.exact_not_name,
                                              startswith=config.startswith, endswith=config.endswith)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Give warning
            if self.is_retrieved(simulation_name): log.warning("Simulation '" + simulation_name + "' is already retrieved. Changing simulation settings will not have any effect on the output")
            else: log.warning("Simulation settings can be adapted but this will possibly not affect anything about the simulation output")

            # Steal
            self.steal_simulation_settings(simulation_name, from_simulation_name, confirm=config.confirm, settings=settings)

    # -----------------------------------------------------------------

    def steal_simulation_settings(self, simulation_name, from_simulation_name, successive=False, confirm=True,
                                  contains=None, not_contains=None, exact_name=None, exact_not_name=None,
                                  startswith=None, endswith=None, settings=None):

        """
        This function ...
        :param simulation_name:
        :param from_simulation_name:
        :param successive:
        :param confirm:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :param settings:
        :return:
        """

        # Get the settings if necessary
        if settings is None: settings = select_simulation_settings(successive=successive, contains=contains, not_contains=not_contains,
                                                          exact_name=exact_name, exact_not_name=exact_not_name, startswith=startswith,
                                                          endswith=endswith)

        # Get the simulations
        simulation = self.get_simulation(simulation_name)
        from_simulation = self.get_simulation(from_simulation_name)

        # Loop over the settings
        for name in settings:

            # Get the values
            old_value = getattr(simulation, name)
            value = getattr(from_simulation, name)

            # Check
            if old_value == value:
                log.success("Property '" + name + "' is equal for both simulations: skipping ...")
                continue

            # Confirm?
            if confirm:

                print("")
                print("Are you sure you want to change the property '" + name + "' of simulation '" + simulation_name + "' from:")
                print(" - " + tostr(old_value))
                print("to:")
                print(" - " + tostr(value) + "?")

                # Prompt for confirmation
                if not prompt_yn("change", "change value", default=True): continue

            # Debugging
            else: log.debug("Changing property '" + name + "' of simulation '" + simulation_name + "' from '" + tostr(old_value) + "' to '" + tostr(value) + "' ...")

            # Set the value
            setattr(simulation, name, value)

            # Show success
            log.success("Property '" + name + "' adapted")

        # Debugging
        log.debug("Saving the simulation ...")

        # Save the simulation
        simulation.save()

    # -----------------------------------------------------------------

    @lazyproperty
    def steal_analysis_options_definition(self):
        
        """
        This function ...
        :return: 
        """
        
        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Select certain settings
        definition.add_positional_optional("matching", "string", "only select options with a name matching this string")
        definition.add_optional("contains", "string", "only select options containing this string in their name")
        definition.add_optional("not_contains", "string", "don't select options containing this string in their name")
        definition.add_optional("exact_name", "string", "only select options with this exact string as their name")
        definition.add_optional("exact_not_name", "string", "don't select options with this exact string as their name")
        definition.add_optional("startswith", "string", "only select options whose name starts with this string")
        definition.add_optional("endswith", "string", "only select options whose name starts with this string")
        
        # Flags
        definition.add_flag("successive", "select options successively", False)
        definition.add_flag("hierarchic", "select properties first, then properties of sections", True)
        definition.add_flag("confirm", "confirm before adapting", True)
        
        # Return the definition
        return definition
        
    # -----------------------------------------------------------------

    def steal_analysis_options_command(self, command, **kwargs):
        
        """
        This function ...
        :param command: 
        :param kwargs: 
        :return: 
        """

        # Set command name
        name = _steal_command_name + " " + _analysis_command_name

        # Get simulation names and config
        splitted, simulation_names, from_simulation_name, config = self.parse_simulations_simulation_command(command, self.steal_analysis_options_definition, name=name, **kwargs)

        # Check settings
        if config.matching is not None:
            if config.contains is not None: raise ValueError("Cannot specify both matching string and containing string")
            config.contains = config.matching

        # Get the options here (so they don't have to be selected for each stealing simulation)
        options = select_analysis_options(successive=config.successive, hierarchic=config.hierarchic, contains=config.contains,
                                          not_contains=config.not_contains, exact_name=config.exact_name,
                                          exact_not_name=config.exact_not_name, startswith=config.startswith, endswith=config.endswith)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check
            if self.is_analysed(simulation_name): log.warning("Simulation '" + simulation_name + "' is already analysed. Simulation will have to be re-analysed")

            # Steal
            self.steal_analysis_options(simulation_name, from_simulation_name, confirm=config.confirm, options=options)

    # -----------------------------------------------------------------

    def steal_analysis_options(self, simulation_name, from_simulation_name, successive=False, hierarchic=True,
                               confirm=True, contains=None, not_contains=None, exact_name=None, exact_not_name=None,
                               startswith=None, endswith=None, options=None):

        """
        This function ...
        :param simulation_name:
        :param from_simulation_name:
        :param successive:
        :param hierarchic:
        :param confirm:
        :param contains:
        :param not_contains:
        :param exact_name:
        :param exact_not_name:
        :param startswith:
        :param endswith:
        :param options:
        :return:
        """
        
        # Get the options if necessary
        if options is None: options = select_analysis_options(successive=successive, hierarchic=hierarchic, contains=contains,
                                              not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name,
                                              startswith=startswith, endswith=endswith)

        # Get the simulations
        simulation = self.get_simulation(simulation_name)
        from_simulation = self.get_simulation(from_simulation_name)

        # Loop over the options
        for spec in options:

            # Property of section
            if types.is_string_tuple(spec):

                # Get the values
                section_name = spec[0]
                property_name = spec[1]
                old_value = simulation.analysis[section_name][property_name]
                value = from_simulation.analysis[section_name][property_name]

                # Check
                if old_value == value:
                    log.success("The " + section_name + " property '" + property_name + "' is equal for both simulations: skipping ...")
                    continue

                # Confirm?
                if confirm:

                    print("")
                    print("Are you sure you want to change the " + section_name + " property '" + property_name + "' of simulation '" + simulation_name + "' from:")
                    print(" - " + tostr(old_value))
                    print("to:")
                    print(" - " + tostr(value) + "?")

                    # Prompt for confirmation
                    if not prompt_yn("change", "change value", default=True): continue

                # Debugging
                else: log.debug("Changing the " + section_name + " property '" + property_name + "' of simulation '" + simulation_name + "' from '" + tostr(old_value) + "' to '" + tostr(value) + "' ...")

                # Set value
                simulation.analysis[section_name][property_name] = value

                # Show success
                log.success("Property '" + property_name + "' adapted")

            # Property
            elif types.is_string_type(spec):

                # Get the values
                old_value = simulation.analysis[spec]
                value = from_simulation.analysis[spec]

                # Check
                if old_value == value:
                    log.success("Property '" + spec + "' is equal for both simulations: skipping ...")
                    continue

                # Confirm?
                if confirm:

                    print("")
                    print("Are you sure you want to change the analysis property '" + spec + "' of simulation '" + simulation_name + "' from:")
                    print(" - " + tostr(old_value))
                    print("to:")
                    print(" - " + tostr(value) + "?")

                    # Prompt for confirmation
                    if not prompt_yn("change", "change value", default=True): continue

                # Debugging
                else: log.debug("Changing the analysis property '" + spec + "' of simulation '" + simulation_name + "' from '" + tostr(old_value) + "' to '" + tostr(value) + "' ...")

                # Set value
                simulation.analysis[spec] = value

                # Show success
                log.success("Property '" + spec + "' adapted")

            # Invalid
            else: raise RuntimeError("Something went wrong")

        # Debugging
        log.debug("Saving the simulation ...")

        # Save the simulation
        simulation.save()

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
        definition.add_flag("new_logging", "specify new logging options", False)
        definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)
        definition.add_flag("cancel", "cancel (or stop) the original simulations", True)

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

        # Set options
        parallelization = config.parallelization
        if config.new_logging: logging_options = config.logging
        else: logging_options = None # original logging options will be used
        scheduling_options = config.scheduling

        # Loop over the simulation names
        for simulation_name in simulation_names:

            # Check
            if simulation_name in self.moved_simulation_names: raise ValueError("Simulation '" + simulation_name + "' is already moved in this session")

            # Check status
            if self.is_finished(simulation_name): raise ValueError("Simulation is already finished")
            if self.is_running(simulation_name): log.warning("Simulation is already running")

            # Get host and check
            if config.host == self.get_simulation(simulation_name).host: raise ValueError("Simulation '" + simulation_name + "' is already queued/running on host '" + tostr(config.host) + "'")

            # Move simulation
            self.move_simulation(simulation_name, config.host, parallelization, logging_options, scheduling_options, cancel=config.cancel)

    # -----------------------------------------------------------------

    def move_simulation(self, simulation_name, host, parallelization=None, logging_options=None, scheduling_options=None, cancel=True):

        """
        This function ...
        :param simulation_name:
        :param host:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :param cancel: cancel or stop the original simulation
        :return:
        """

        # Debugging
        log.debug("Moving simulation '" + simulation_name + "' to host '" + tostr(host) + "' ...")

        # Get logging options from previous launch
        if logging_options is None: logging_options = self.get_logging_options_for_simulation(simulation_name)
        # Parallelization and scheduling options will have to be determined automatically for the new host if not specified

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Cancel or stop the simulation
        if cancel: self.cancel_or_stop_simulation(simulation_name)

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

    def cancel_or_stop_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        if self.is_running(simulation_name): self.stop_simulation(simulation_name)
        elif self.is_queued(simulation_name): self.cancel_simulation(simulation_name)
        else: log.warning("Simulation '" + simulation_name + "' is not running or queued: cancelling is not required")

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

        # Check whether the remote ski file exists
        if remote.is_file(ski_path):

            # Debugging
            log.debug("Removing the ski file '" + ski_path + "' from remote host '" + remote.host_id + "' ...")

            # Remove the remote ski file
            remote.remove_file(ski_path)

        # Ski file is not present (anymore)
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

    def unlaunch_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command)

        # Unlaunch
        self.unlaunch_simulation(simulation_name)

    # -----------------------------------------------------------------

    def unlaunch_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Unlaunching simulation '" + simulation_name + "' ...")

        # Unfinish: remove remote simulation output and unset finished flag
        self.unfinish_simulation(simulation_name)

        # Unretrieve: remove local simulation output and unset retrieved flag
        self.unretrieve_simulation(simulation_name)

        # Unanalyse: remove analysis output and unset analysis flags
        self.unanalyse_simulation(simulation_name)

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
        if self.has_output(simulation_name): self.clear_simulation_output_local(simulation_name)

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

            # Check whether the directory exists
            if not fs.is_directory(dirpath): raise RuntimeError("The directory '" + dirpath + "' does not exist")

            # Determine relative directory path
            relpath = fs.relative_to(dirpath, directory_path)

            # Determine new path
            new_path = fs.join(cache_path, relpath)

            # Get directory name and containing directory
            dirname = fs.name(dirpath)
            cache_directory_path = fs.directory_of(new_path)
            if not fs.is_directory(cache_directory_path): fs.create_directory(cache_directory_path, recursive=True)

            # Debugging
            from_name = fs.directory_of(relpath)
            if from_name == "": from_name = fs.directory_of(dirpath)
            log.debug("Caching '" + dirname + "' directory from '" + from_name + "' to '" + cache_directory_path + "' ...")

            # Copy the directory
            fs.copy_directory(dirpath, cache_directory_path, replace_files=True, replace_directories=False) # avoid replacing filed directories with empty ones

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
            cache_directory_path = fs.directory_of(new_path)
            if not fs.is_directory(cache_directory_path): fs.create_directory(cache_directory_path, recursive=True)

            # Debugging
            from_name = fs.directory_of(relpath)
            if from_name == "": from_name = fs.directory_of(filepath)
            log.debug("Caching '" + filename + "' file from '" + from_name + "' to '" + cache_directory_path + "' ...")

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
        definition.add_flag("new_logging", "adapt the logging options", False)
        definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)
        definition.add_optional("screen_or_job", "string", "specifies whether to launch in screen or job if this information is unavailable (missing execution handle)")

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

        # Check
        if simulation_name in self.relaunched_simulation_names: raise ValueError("Simulation '" + simulation_name + "' is already scheduled for relaunching")

        # Check simulation
        if not self.is_failed(simulation_name):
            #print(self.get_status(simulation_name))
            if config.finished and self.is_finished(simulation_name):
                log.info("Simulation '" + simulation_name + "' is already finished, removing local and remote output and undoing analysis ...")
                if self.is_retrieved(simulation_name): self.unretrieve_simulation(simulation_name)
                self.clear_simulation_output_remote(simulation_name)
                self.unanalyse_simulation(simulation_name)
            elif self.has_unknown_status(simulation_name): log.warning("The current status of simulation '" + simulation_name + "' is unknown")
            else: raise ValueError("Simulation '" + simulation_name + "' is running, finished or still queued")

        # Set the options
        parallelization = config.parallelization
        if config.new_logging: logging_options = config.logging
        else: logging_options = None
        scheduling_options = config.scheduling

        # Relaunch the simulation
        self.relaunch_simulation(simulation_name, parallelization=parallelization, logging_options=logging_options, 
                                 scheduling_options=scheduling_options, screen_or_job=config.screen_or_job)

    # -----------------------------------------------------------------

    def relaunch_simulation(self, simulation_name, parallelization=None, logging_options=None, scheduling_options=None,
                            screen_or_job=None):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :param screen_or_job:
        :return:
        """

        # Debugging
        log.debug("Relaunching simulation '" + simulation_name + "' ...")

        # Screen or job?
        if screen_or_job is None:
            if self.is_screen_execution(simulation_name): screen_or_job = "screen"
            elif self.is_job_execution(simulation_name): screen_or_job = "job"
            else: raise NotImplementedError("Execution handle not supported")

        # Screen
        if screen_or_job == "screen": self.relaunch_simulation_screen(simulation_name, parallelization=parallelization, logging_options=logging_options)

        # Job
        elif screen_or_job == "job": self.relaunch_simulation_job(simulation_name, parallelization=parallelization, logging_options=logging_options, scheduling_options=scheduling_options)

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
            if not self.config.dry:

                # Debugging
                log.debug("Stopping the screen session '" + screen_name + "' ...")

                # Stop the screen
                remote.kill_screen(screen_name)

            # Dry mode: don't stop
            else: log.warning("[DRY] Not stopping the screen session '" + screen_name + "' ...")

            # Check which simulations are still queued or running in this screen
            screen_simulation_names = self.get_queued_and_running_simulation_names_for_screen(screen_name)

            # Loop over the simulations
            for screen_simulation_name in screen_simulation_names:

                # Stop/cancel and re-add the simulations to the queue (they will be launched together with the relaunched simulations in a new screen)
                #self.cancel_or_stop_simulation(screen_simulation_name): NO, THE SCREEN IS ALREADY STOPPED

                # Add the simulation to the queue
                self._relaunch_simulation_screen_impl(screen_simulation_name) # use original parallelization and logging options

        # Relaunch
        self._relaunch_simulation_screen_impl(simulation_name, parallelization=parallelization, logging_options=logging_options)
        
    # -----------------------------------------------------------------

    def _relaunch_simulation_screen_impl(self, simulation_name, parallelization=None, logging_options=None):

        """
        This function ...
        :param simulation_name:
        :param parallelization:
        :param logging_options:
        :return: 
        """

        # Set the parallelization scheme
        if parallelization is None: parallelization = self.get_parallelization_for_simulation(simulation_name)

        # Set the logging options
        if logging_options is None: logging_options = self.get_logging_options_for_simulation(simulation_name)

        # Unlaunch: remove remote and local output and analysis output
        self.unlaunch_simulation(simulation_name)

        # Get the simulation
        simulation = self.get_simulation(simulation_name)
        host_id = self.get_host_id_for_simulation(simulation_name)

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
        if parallelization is None: parallelization = self.get_parallelization_for_simulation(simulation_name)

        # Set the logging options
        if logging_options is None: logging_options = self.get_logging_options_for_simulation(simulation_name)

        # Set the scheduling options
        if scheduling_options is None: scheduling_options = self.get_scheduling_options_for_simulation(simulation_name)

        # Unlaunch: remove remote and local output and analysis output
        self.unlaunch_simulation(simulation_name)

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

    @property
    def nnew(self):

        """
        This function ...
        :return:
        """

        return len(self.new)

    # -----------------------------------------------------------------

    @property
    def has_new(self):

        """
        This function ...
        :return:
        """

        return self.nnew > 0

    # -----------------------------------------------------------------

    @property
    def new_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.new.simulation_names

    # -----------------------------------------------------------------

    def is_relaunched(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.relaunched_simulation_names

    # -----------------------------------------------------------------

    def is_new(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self.new_simulation_names

    # -----------------------------------------------------------------

    def set_simulation_options_from_original(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Get the original simulation
        original_simulation = self.get_simulation(simulation.name)

        # Analyser paths: NO, this is in the analysis options now!
        #simulation.analyser_paths = original_simulation.analyser_paths

        # Options for retrieval
        simulation.retrieve_types = original_simulation.retrieve_types

        # Options for removing remote or local input and output
        simulation.remove_remote_input = original_simulation.remove_remote_input
        simulation.remove_remote_output = original_simulation.remove_remote_output  # After retrieval
        simulation.remove_remote_simulation_directory = original_simulation.remove_remote_simulation_directory  # After retrieval
        simulation.remove_local_output = original_simulation.remove_local_output  # After analysis

    # -----------------------------------------------------------------

    @lazyproperty
    def retry_simulations_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Which simulations?
        definition.add_required("case", "string", "which simulations to relaunch (based on their status or reason for failing)", choices=failed_cases)

        # Increase walltime
        definition.add_optional("walltime_factor", "positive_real", "factor with which to increase the walltime for simulations which were aborted because they exceeded their specified walltime", 1.2)
        definition.add_optional("walltime", "duration", "new walltime for the simulations")

        # Return
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def retry_simulations_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    @property
    def failed_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.status.failed_names

    # -----------------------------------------------------------------

    @property
    def aborted_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.status.aborted_names

    # -----------------------------------------------------------------

    @property
    def crashed_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.status.crashed_names

    # -----------------------------------------------------------------

    @property
    def cancelled_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.status.cancelled_names

    # -----------------------------------------------------------------

    @lazyproperty
    def exceeded_walltime_simulation_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.aborted_simulation_names:
            if self.has_exceed_walltime(simulation_name): names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def exceed_memory_simulation_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for simulation_name in self.aborted_simulation_names: pass
        return names

    # -----------------------------------------------------------------

    def retry_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Update
        #kwargs.update(self.retry_simulations_kwargs)

        # Get the configuration
        config = self.get_config_from_command(command, self.retry_simulations_definition, **kwargs)

        # The scheduling options for the different simulations
        scheduling = dict()

        # All failed simulations
        if config.case == failed_name: simulation_names = self.failed_simulation_names

        # Aborted simulations
        elif config.case == aborted_name:

            # Get the names of the simulations
            simulation_names = self.aborted_simulation_names

            # Check for simulations aborted because they exceed the walltime
            for simulation_name in simulation_names:
                
                if self.has_exceed_walltime(simulation_name):

                    if config.walltime is not None:

                        # Check the new walltime
                        log.debug("Simulation '" + simulation_name + "' was aborted because the walltime was exceeded")
                        scheduling_options = self.get_scheduling_options_for_simulation(simulation_name).copy()
                        if scheduling_options.walltime >= config.walltime: raise ValueError("Invalid walltime: must be higher than the previous walltime of " + tostr(scheduling_options.walltime) + " seconds")
                        scheduling_options.walltime = config.walltime

                    else:

                        # Determine the new walltime
                        log.debug("Simulation '" + simulation_name + "' was aborted because the walltime was exceeded, increasing the walltime by a factor of " + str(config.walltime_factor) + " ...")
                        scheduling_options = self.get_scheduling_options_for_simulation(simulation_name).copy()
                        log.debug("Previous walltime: " + tostr(scheduling_options.walltime))
                        scheduling_options.walltime *= config.walltime_factor
                        log.debug("New walltime: " + tostr(scheduling_options.walltime))

                    # Set the scheduling options for this simulation
                    scheduling[simulation_name] = scheduling_options

        # Crashed simulations
        elif config.case == crashed_name: simulation_names = self.crashed_simulation_names

        # Cancelled simulations
        elif config.case == cancelled_name: simulation_names = self.cancelled_simulation_names

        # Exceeded the walltime
        elif config.case == exceeded_walltime_name:

            # Get the names of the simulations
            simulation_names = self.exceeded_walltime_simulation_names

            # Increase the walltimes
            for simulation_name in simulation_names:

                # New walltime is given, check
                if config.walltime is not None:

                    # Check the new walltime
                    scheduling_options = self.get_scheduling_options_for_simulation(simulation_name).copy()
                    log.debug("Simulation '" + simulation_name + "' exceeded the walltime of " + tostr(scheduling_options.walltime) + " seconds")
                    if scheduling_options.walltime >= config.walltime: raise ValueError("Invalid walltime: must be higher than the previous walltime of " + tostr(scheduling_options.walltime) + " seconds")
                    scheduling_options.walltime = config.walltime

                # Determine new walltime
                else:

                    scheduling_options = self.get_scheduling_options_for_simulation(simulation_name).copy()
                    log.debug("Simulation '" + simulation_name + "' exceeded the walltime of " + tostr(scheduling_options.walltime) + ", increasing by a factor of " + str(config.walltime_factor) + " ...")
                    scheduling_options.walltime *= config.walltime_factor
                    log.debug("New walltime: " + tostr(scheduling_options.walltime))

                # Set the scheduling options for this simulation
                scheduling[simulation_name] = scheduling_options

        # Exceeded memory
        elif config.case == exceeded_memory_name: raise NotImplementedError("Not implemented")

        # Invalid
        else: raise ValueError("Invalid case: '" + config.case + "'")

        # Relaunch each simulation
        for simulation_name in simulation_names:

            # Get scheduling options
            if simulation_name in scheduling: scheduling_options = scheduling[simulation_name]
            else: scheduling_options = None

            # Relaunch
            self.relaunch_simulation(simulation_name, scheduling_options=scheduling_options)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching simulations ...")

        # Set the remote input paths
        if self.config.shared_input: self.set_shared_remote_input_paths()

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

        # Set new simulation IDs
        if self.has_new_simulations: self.set_new_simulation_ids()

        # Add the new simulations
        if self.has_new_simulations: self.add_new_simulations()

        # Debugging
        log.debug("Creating backup of assignment scheme ...")

        # Backup assignment?
        if self.config.backup_assignment: self.assignment.saveto(self.backup_assignment_path)

        # Debugging
        log.debug("Updating the assignment scheme ...")

        # Change assignment
        for simulation in self.launched_simulations: self.assignment.add_or_update_simulation(simulation, success=True)

        # Set that the assignment has changed
        self._adapted_assignment = True

        # Debugging
        log.debug("Updating the simulation status table ...")

        # Add the new simulation names to the status table
        if self.has_new_simulations: self.add_new_simulations_to_status()

        # Reset status?
        for simulation in self.launched_simulations: self.reset_status_for_simulation(simulation.name)

    # -----------------------------------------------------------------

    def add_new_simulations(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adding the new simulations to the set of simulations ...")

        # Loop over only the new simulations
        for simulation in self.launched_simulations:
            if simulation.name not in self.new_simulation_names: continue

            # Add the simulation
            self.add_simulation(simulation)

    # -----------------------------------------------------------------

    def add_new_simulations_to_status(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adding the new simulations to the status table ...")

        # Create new status table
        self.status = SimulationStatusTable.from_previous(self.status, new_simulations=self.new_simulation_names)

    # -----------------------------------------------------------------

    def set_shared_remote_input_paths(self):

        """
        This function ...
        :return: 
        """

        # Debugging
        log.debug("Setting remote input paths for the simulations to be launched ...")

        # Loop over the launch host IDs
        for host_id in self.launcher.queue_host_ids:

            # Get the remote
            remote_input_path = None
            remote = self.get_remote(host_id)

            # Find remote input path
            # Loop over the simulations
            for simulation_name in self.get_simulation_names_for_host_id(host_id):

                # Get the simulation
                simulation = self.get_simulation(simulation_name)
                if simulation.remote_input_path is None: continue

                # Check remote input path
                if remote.is_directory(simulation.remote_input_path) and not remote.is_empty(simulation.remote_input_path):
                    remote_input_path = simulation.remote_input_path
                    break

            # Show
            if remote_input_path is None: continue
            else: log.debug("The remote simulation input path for host '" + host_id + "' is '" + remote_input_path + "'")

            # Set remote input path
            self.set_remote_input_path_for_host(host_id, remote_input_path)

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

    def set_new_simulation_ids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the simulation IDs of the new simulations ...")

        # Set new IDs
        for simulation in self.new_simulations:

            # Set the new ID of the new simulation
            self.new.set_id(simulation.name, simulation.id)

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

    @lazyproperty
    def new_simulations(self):
        
        """
        This function ...
        :return: 
        """
        
        simulations = []
        for simulation in self.launched_simulations:
            if not self.is_new(simulation.name): continue
            simulations.append(simulation)
        return simulations
        
    # -----------------------------------------------------------------

    @property
    def nnew_simulations(self):
        
        """
        This function ...
        :return: 
        """

        return len(self.new_simulations)

    # -----------------------------------------------------------------

    @property
    def has_new_simulations(self):

        """
        This function ...
        :return:
        """

        return self.nnew_simulations > 0

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
        #print(self.assignment)
        print("")
        fmt.print_table(self.assignment)
        print("")

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

    def add_host_id(self, host_id, name=None):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        # Add
        self.remotes.add_host_id(host_id, name=name)

    # -----------------------------------------------------------------

    def add_host_ids(self, host_ids):

        """
        This function ...
        :param host_ids:
        :return:
        """

        # Add multiple
        self.remotes.add_host_ids(host_ids)

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
        else: id_string = str(simulation.id)

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

    @property
    def has_simulations(self):
        
        """
        This function ...
        :return: 
        """
        
        return self.nsimulations > 0
        
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
        print(fmt.bold + "Number of running simulations: " + fmt.reset + str(self.nrunning) + " (" + tostr(self.percentage_nrunning, round=True, ndigits=2) + "%)")
        print("")

        # If path is given, initialize table
        #if path is not None: table = SmartTable()
        #else: table = None
        table = None
        #print(self.status)
        #for index in range(len(self.status)):
        #    status = self.status.get_value("Status", index)
        #    print(status)

        # Print in columns
        with fmt.print_in_columns() as print_row:

            # Set the column names and units
            column_names = []
            column_units = []
            for name in self.status_column_names: column_names.append(fmt.bold + name + fmt.reset)
            for unit in self.status_column_units: column_units.append("[" + unit + "]" if unit != "" else "")
            if extra is not None:
               for col in extra:
                   column_names.append(fmt.bold + extra_column_names[col] + fmt.reset)
                   if col in extra_column_units: column_units.append("[" + extra_column_units[col] + "]" if extra_column_units[col] != "" else "")
                   else: column_units.append("")

            # Show the header
            print_row(*column_names)
            #if self.has_info: print_row(*column_units)
            if not sequences.all_none(column_units): print_row(*column_units)

            # Initialize columns
            if path is not None: columns = [[] for _ in column_names]
            else: columns = None

            # Loop over the simulations
            for index, simulation_name in enumerate(self.simulation_names):

                # Get the simulation
                simulation = self.get_simulation(simulation_name)

                # Get the status
                status = self.get_status(simulation_name)
                #print(simulation_name, status)

                # Get index string
                index_string = "[" + strings.integer(index, 3, fill=" ") + "] "

                # Set color
                if is_analysing_or_analysed_status(status): color = "green"
                elif is_retrieved_status(status): color = "yellow"
                elif is_finished_status(status): color = "yellow"
                elif is_running_status(status): color = None
                else: color = "red"

                # Create strings
                if isinstance(simulation, RemoteSimulation):
                    host_string = tostr(simulation.host)
                    id_string = tostr(simulation.id)
                elif isinstance(simulation, SkirtSimulation): host_string = id_string = "--"
                else: raise ValueError("Invalid simulation object")

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
                        value = self.get_extra_info_for_simulation(simulation_name, col)
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

    def get_extra_info_for_simulation(self, simulation_name, name):

        """
        This function ...
        :param simulation_name:
        :param name:
        :return:
        """

        # Get the value
        if name == _screen_extra_name: value = self.get_screen_name(simulation_name)
        elif name == _job_extra_name: value = self.get_job_id_string(simulation_name)
        elif name == _disk_extra_name: value = self.get_disk_size(simulation_name)
        elif name == _runtime_extra_name: value = self.get_runtime(simulation_name)
        elif name == _memory_extra_name: value = self.get_peak_memory(simulation_name)
        elif name == _elapsed_extra_name: value = self.get_elapsed_time(simulation_name)
        elif name == _cached_output_name: value = self.is_cached_output(simulation_name)
        elif name == _cached_extraction_name: value = self.is_cached_extraction(simulation_name)
        elif name == _cached_plotting_name: value = self.is_cached_plotting(simulation_name)
        elif name == _cached_misc_name: value = self.is_cached_misc(simulation_name)
        else: raise ValueError("Invalid extra column name: '" + name + "'")

        # Add unit?
        if name in extra_column_units and value is not None: value = value.to(extra_column_units[name]).value

        # Return
        return value

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

    def show_cpu_command(self, command, **kwargs):
        
        """
        This function ...
        :param command: 
        :param kwargs: 
        :return: 
        """

        # Get simulation names
        simulation_names = self.get_simulation_names_from_command(command, **kwargs)

        # Show
        self.show_cpu(simulation_names)

    # -----------------------------------------------------------------

    def show_cpu(self, simulation_names):
        
        """
        This function ...
        :param simulation_names: 
        :return: 
        """

        # Debugging
        log.debug("Showing CPU usage for the simulations ...")

        # Keep a list of the runtimes and CPU times
        runtimes = []
        cputimes = []

        # Print in columns
        print("")
        with fmt.print_in_columns() as print_row:

            # Set the column names and units
            column_names = ["Simulation name", "Host ID", "Number of cores", "Total runtime", "Total CPU time"]
            column_units = ["", "", "", "h", "h"]

            # Show the header
            print_row(*column_names, bold=True)
            print_row(*column_units)

            # Loop over the simulations
            for simulation_name in simulation_names:

                # Get the host ID
                host_id = self.get_host_id_for_simulation(simulation_name)

                # Get the number of cores
                ncores = self.get_ncores_for_simulation(simulation_name)

                # Finished simulation?
                if self.is_finished(simulation_name):

                    # Get the runtime of the simulation, in hours
                    runtime = self.get_runtime(simulation_name).to("h").value

                    # Get the CPU time
                    cpu_time = runtime * ncores

                    # Add
                    runtimes.append(runtime)
                    cputimes.append(cpu_time)

                # Not finished
                else: runtime = cpu_time = "--"

                # Show the row for the simulation
                parts = [simulation_name, host_id, ncores, tostr(runtime), tostr(cpu_time)]
                print_row(*parts)

        # Show the averages (and total)
        print("")
        print(" - average runtime: " + tostr(numbers.arithmetic_mean(*runtimes)) + " h")
        print(" - average CPU time: " + tostr(numbers.arithmetic_mean(*cputimes)) + " h")
        print(" - total CPU time: " + tostr(sequences.sum(cputimes)) + " h")
        print("")

    # -----------------------------------------------------------------

    def show_timeline_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command, **kwargs)

        # Show
        self.show_timeline(simulation_name)

    # -----------------------------------------------------------------
    
    def show_timeline(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get timeline
        timeline = self.get_timeline_for_simulation(simulation_name)
        #print(timeline)

        # Show stuff
        print("")
        print(" - total runtime: " + tostr(timeline.total / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - setup time: " + tostr(timeline.setup / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - stellar emission time: " + tostr(timeline.stellar / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - spectra calculation time: " + tostr(timeline.spectra / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - dust emission time: " + tostr(timeline.dust / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - writing time: " + tostr(timeline.writing / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - communication time: " + tostr(timeline.communication / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - waiting time: " + tostr(timeline.waiting / 60., scientific=False, round=True, ndigits=2) + " min")
        print(" - other time spent: " + tostr(timeline.other / 60., scientific=False, round=True, ndigits=2) + " min")
        print("")

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

    @property
    def history_filename(self):

        """
        This function ...
        :return:
        """

        return "manager.dat"

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

    @lazyproperty
    def plot_scaling_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def plot_scaling_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        from ..plot.scaling import ScalingPlotter

        # Get config
        #config =

        # Create plotter
        plotter = ScalingPlotter()

        # Run the plotter
        plotter.run(timing=self.timing, memory=self.memory)

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

        # Add settings
        definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps_and_extra, default=all_steps_and_extra)
        definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")
        definition.add_optional("not_steps", "string_list", "don't analyse these steps", choices=all_steps_and_extra)
        definition.add_optional("not_features", "string_list", "don't analyse these features (if a single not_step is defined)")
        definition.import_section("analysis", "options for the simulation analyser", analyse_simulation_definition)

        # Prompt for extra?
        definition.add_flag("prompt_extra_configs", "prompt configuration settings for extra analyser classes")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def analysed_any(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_any

    # -----------------------------------------------------------------

    def analysed_any_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_any_extraction

    # -----------------------------------------------------------------

    def analysed_all_extraction(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_all_extraction

    # -----------------------------------------------------------------

    def analysed_any_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_any_plotting

    # -----------------------------------------------------------------

    def analysed_all_plotting(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_all_plotting

    # -----------------------------------------------------------------

    def analysed_any_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_any_misc

    # -----------------------------------------------------------------

    def analysed_all_misc(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_all_misc

    # -----------------------------------------------------------------

    def analysed_batch(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_batch

    # -----------------------------------------------------------------

    def analysed_scaling(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_scaling

    # -----------------------------------------------------------------

    def analysed_any_extra(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_any_extra

    # -----------------------------------------------------------------

    def analysed_all_extra(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return self.get_simulation(simulation_name).analysed_all_extra

    # -----------------------------------------------------------------

    def prompt_extra_analysis_configs(self, simulation_names):

        """
        This function ...
        :param simulation_names:
        :return:
        """

        extra_configs = dict()

        # Get all unique analyser class names
        analyser_classes = defaultdict(list)
        for simulation_name in simulation_names:
            simulation = self.get_simulation(simulation_name)
            for class_name in simulation.analyser_class_names: analyser_classes[class_name].append(simulation_name)

        # Loop over the class names
        for class_name in analyser_classes:

            # Prompt
            if not prompt_proceed("set configuration for '" + class_name + "'", default=True): continue

            # Initialize config for class
            class_config = dict()

            # Keep adding settings
            while True:

                # Get name
                name = prompt_string("name", "name of the setting (press ENTER to finish configuring)")
                if not name: break

                # Get ptype
                ptype = prompt_string("ptype", "type of the setting")

                # Get value
                value = prompt_variable(name, ptype, "value of the setting")

                # Add value
                class_config[name] = value

            # Add config for class
            extra_configs[class_name] = class_config

        # Return the configs
        return extra_configs

    # -----------------------------------------------------------------

    def reanalyse_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation names and config
        simulation_names, config = self.get_simulation_names_and_config_from_command(command, command_definition=self.reanalyse_simulation_definition, **kwargs)

        # Get settings
        steps = config.steps
        features = config.features
        not_steps = config.not_steps
        not_features = config.not_features

        # Set extra configs?
        if config.prompt_extra_configs: extra_configs = self.prompt_extra_analysis_configs(simulation_names)

        # No extra configs
        else: extra_configs = None

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Check
            if not self.analysed_any(simulation_name): raise ValueError("The simulation '" + simulation_name + "' has not been analysed yet")

            # Reanalyse the simulation
            self.reanalyse_simulation(simulation_name, steps=steps, features=features, not_steps=not_steps,
                                      not_features=not_features, config=config.analysis, extra_configs=extra_configs)

    # -----------------------------------------------------------------

    def reanalyse_simulation(self, simulation_name, steps=None, features=None, not_steps=None, not_features=None,
                             config=None, extra_configs=None):

        """
        This function ...
        :param simulation_name:
        :param steps:
        :param features:
        :param not_steps:
        :param not_features:
        :param config:
        :param extra_configs:
        :return:
        """

        # Debugging
        log.debug("Re-analysing simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Reanalyse simulation
        reanalyse_simulation(simulation, steps, features, not_steps=not_steps, not_features=not_features, config=config, extra_configs=extra_configs)

        # Reset the status
        self.reset_status_for_simulation(simulation_name)

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

            # Cache now?
            if self.do_caching and self.config.cache_after_analysis: self.cache_simulation(simulation.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Get the definition
        definition = analyse_simulation_definition.copy(pos_optional=False)

        # Prompt for extra?
        definition.add_flag("prompt_extra_configs", "prompt configuration settings for extra analyser classes")

        # Return
        return definition

    # -----------------------------------------------------------------

    def analyse_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get simulation names and config
        simulation_names, config = self.get_simulation_names_and_config_from_command(command, command_definition=self.analyse_simulation_definition, **kwargs)

        # Set extra configs?
        if config.get("prompt_extra_configs", False): extra_configs = self.prompt_extra_analysis_configs(simulation_names)
        else: extra_configs = None
        if "prompt_extra_configs" in config: config.pop("prompt_extra_configs")

        # Analyse simulations
        for simulation_name in simulation_names:

            # Check
            if not self.is_retrieved(simulation_name): raise ValueError("Simulation '" + simulation_name + "' is not yet retrieved")

            # Analyse
            self.analyse_simulation(simulation_name, config=config, extra_configs=extra_configs)

    # -----------------------------------------------------------------

    def analyse_simulation(self, simulation_name, config=None, extra_configs=None):

        """
        This function ...
        :param simulation_name:
        :param config:
        :param extra_configs:
        :return:
        """

        # Debugging
        log.debug("Analysing simulation " + simulation_name + " ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Analyse the simulation
        analyse_simulation(simulation, config=config, extra_configs=extra_configs)

        # Reset the status
        self.reset_status_for_simulation(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def analyse_all_definition(self):

        """
        This function ...
        :return:
        """

        return analyse_simulation_definition.copy(pos_optional=False)

    # -----------------------------------------------------------------

    def analyse_all_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get config
        config = self.get_config_from_command(command, self.analyse_all_definition, **kwargs)

        # Analyse all retrieved simulations
        self.analyse_all_simulations(config=config)

    # -----------------------------------------------------------------

    def analyse_all_simulations(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Inform the user
        log.info("Analysing all retrieved simulations ...")

        # Loop over the retrieved but not yet analysed simulations
        for simulation_name in self.all_retrieved_not_analysed_simulation_names:

            # Show steps that will be performed
            show_analysis_steps(self.get_simulation(simulation_name))

            # Analyse simulation
            self.analyse_simulation(simulation_name, config=config)

            # Cache now?
            if self.do_caching and self.config.cache_after_analysis: self.cache_simulation(simulation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def mimic_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Name
        definition.add_required("name", "string", "name for the new simulation")

        # Labeled values
        definition.add_optional("labeled", "dictionary", "new values for labeled properties")

        # Flags
        definition.add_flag("adapt_settings", "adapt simulation settings")
        definition.add_flag("adapt_analysis", "adapt analysis options")

        # Custom options
        definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulation")
        definition.add_flag("mimic_scheduling", "mimic the scheduling options", True)
        definition.add_flag("mimic_analysis", "mimic the analysis options", True)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)
        definition.import_section("analysis", "options for the simulation analyser", analyse_simulation_definition)

        # Host
        definition.add_optional("host", "host", "host for the simulation execution (default is same as original simulation)")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    @lazyproperty
    def mimic_simulation_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["required_to_optional"] = False
        return kwargs

    # -----------------------------------------------------------------

    def mimic_simulation_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Add kwargs
        kwargs.update(self.mimic_simulation_kwargs)

        # Get the simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, self.mimic_simulation_definition, **kwargs)

        # Get
        if config.mimic_scheduling: scheduling_options = None
        else: scheduling_options = config.scheduling
        if config.mimic_analysis: analysis_options = None
        else: analysis_options = config.analysis

        # Mimic the simulation
        self.mimic_simulation(simulation_name, config.name, new_values=config.labeled, parallelization=config.parallelization,
                              analysis_options=analysis_options, scheduling_options=scheduling_options, host=config.host)

    # -----------------------------------------------------------------

    def mimic_simulation(self, simulation_name, new_simulation_name, new_values=None, parallelization=None,
                         analysis_options=None, scheduling_options=None, host=None):
        
        """
        This function ...
        :param simulation_name:
        :param new_simulation_name:
        :param new_values:
        :param parallelization:
        :param analysis_options:
        :param scheduling_options:
        :param host:
        :return: 
        """

        # Check whether the new simulation name is unique
        if new_simulation_name in self.simulation_names: raise ValueError("Simulation name '" + new_simulation_name + "' already in use")

        # Create copy of the original ski file
        ski = self.get_skifile(simulation_name).copy()

        # Adapt labeled values
        if new_values is not None: ski.set_labeled_values(new_values)
        #else: # TODO: use composer to adapt the model

        # Get remote settings
        if host is not None: host_id, cluster_name = host.id, host.cluster_name
        else:
            host_id = self.get_host_id_for_simulation(simulation_name)
            cluster_name = self.get_cluster_name_for_simulation(simulation_name)
        #print(host_id, cluster_name)

        # Debugging
        host_string = host_id + ":" + cluster_name if cluster_name is not None else host_id
        log.debug("Simulation will be executed on host '" + host_string + "' ...")

        # Get parallelization
        if parallelization is None: parallelization = self.get_parallelization_for_simulation(simulation_name)

        # Get original simulation settings
        logging_options = self.get_logging_options_for_simulation(simulation_name).copy()

        # Get analysis options
        if analysis_options is not None: analysis_options = AnalysisOptions(**analysis_options)
        else:

            # Get analysis options of original simulation
            analysis_options = self.get_analysis_options_for_simulation(simulation_name).copy()

            # Fix analysis options
            from ...modeling.fitting.generation import correct_analysis_paths
            correct_analysis_paths(analysis_options, new_simulation_name)

        # Get scheduling options
        if scheduling_options is not None: scheduling_options = SchedulingOptions(**scheduling_options)
        else:

            # Get scheduling options of original simulation
            scheduling_options = self.get_scheduling_options_for_simulation(simulation_name).copy()

            # Fix scheduling options
            if scheduling_options.local_jobscript_path is not None:

                # Replace simulation name
                if simulation_name in scheduling_options.local_jobscript_path: scheduling_options.local_jobscript_path = scheduling_options.local_jobscript_path.replace(simulation_name, new_simulation_name)
                else: raise RuntimeError("Don't know how to adapt the local jobscript path")

        # Create simulation path, save the ski file
        simulation_path = self.get_base_path(simulation_name).replace(simulation_name, new_simulation_name)
        fs.create_directory(simulation_path)
        ski_path = self.get_skifile_path(simulation_name).replace(simulation_name, new_simulation_name)
        ski.saveto(ski_path)

        # Get input paths and output path
        input_paths = self.get_input(simulation_name)
        output_path = self.get_output_path(simulation_name).replace(simulation_name, new_simulation_name)
        fs.create_directory(output_path)

        # Launch simulation
        self.launch_simulation(new_simulation_name, ski_path, host_id=host_id, input_path=input_paths,
                               output_path=output_path, parallelization=parallelization, logging_options=logging_options,
                               scheduling_options=scheduling_options, analysis_options=analysis_options)

    # -----------------------------------------------------------------

    @lazyproperty
    def launch_simulation_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition()

        # Name for the simulation
        definition.add_required("name", "string", "name for the simulation")
        
        # Ski file
        definition.add_required("ski", "file_path", "name/path of the ski file")

        # Remote host
        definition.add_positional_optional("remote", "string", "remote host on which to run the simulation (if none is specified, the simulation is run locally", choices=find_host_ids())

        # Input and output
        definition.add_optional("input", "directory_path", "input directory for the simulation(s)", letter="i")
        definition.add_optional("output", "directory_path", "output directory for the simulation(s)", letter="o")

        # Simulation options
        definition.add_positional_optional("parallelization", "parallelization", "parallelization scheme for the simulation")
        definition.import_section_from_composite_class("logging", "simulation logging options", LoggingOptions)
        definition.import_section_from_composite_class("scheduling", "simulation analysis options", SchedulingOptions)

        # Analysis options
        definition.import_section("analysis", "options for the simulation analyser", analyse_simulation_definition)

        # Options from other simulations
        definition.add_optional("parallelization_from", "integer_or_string", "take parallelization from other simulation")
        definition.add_optional("logging_from", "integer_or_string", "take logging options from other simulation")
        definition.add_optional("scheduling_from", "integer_or_string", "take scheduling options from other simulation")
        definition.add_optional("analysis_from", "integer_or_string", "take the analysis options from other simulation")

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def launch_simulation_command(self, command, **kwargs):

        """
        This function ...
        """
        
        # Get the configuration
        config = self.get_config_from_command(command, self.launch_simulation_definition, **kwargs)

        # Get parallelization
        if config.parallelization_from is not None:
            if config.parallelization is not None: raise ValueError("Cannot specify parallelization and take parallelization from other simulation")
            parallelization = self.get_parallelization_for_simulation(config.parallelization_from)
        else: parallelization = config.parallelization

        # Get logging options
        if config.logging_from is not None: logging_options = self.get_logging_options_for_simulation(config.logging_from)
        else: logging_options = LoggingOptions(**config.logging)

        # Get scheduling options
        if config.scheduling_from is not None: scheduling_options = self.get_scheduling_options_for_simulation(config.scheduling_from)
        else: scheduling_options = SchedulingOptions(**config.scheduling)

        # Get analysis options
        if config.analysis_from is not None:
            analysis_options = self.get_analysis_options_for_simulation(config.analysis_from)
            from ...modeling.fitting.generation import correct_analysis_paths
            correct_analysis_paths(analysis_options, config.name)
        else: analysis_options = AnalysisOptions(**config.analysis)

        # Launch simulation
        self.launch_simulation(config.name, config.ski, host_id=config.remote, input_path=config.input,
                               output_path=config.output,
                               parallelization=parallelization, logging_options=logging_options,
                               scheduling_options=scheduling_options, analysis_options=analysis_options)

    # -----------------------------------------------------------------

    def launch_simulation(self, name, ski_path, host_id=None, input_path=None, output_path=None, parallelization=None,
                          logging_options=None, scheduling_options=None, analysis_options=None):

        """
        This function ...
        :param name: 
        :param ski_path:
        :param host_id:
        :param input_path:
        :param output_path: 
        :param parallelization: 
        :param logging_options:
        :param scheduling_options:
        :param analysis_options:
        :return: 
        """

        # Set output path if not specified
        base_path = fs.directory_of(ski_path)
        if output_path is None:
            output_path = fs.join(base_path, "out")
            if fs.is_directory(output_path): raise IOError("Output directory '" + output_path + "' already exists")
            fs.create_directory(output_path, recursive=True)

        # Make simulation definition
        definition = SingleSimulationDefinition(ski_path, output_path, input_path, name=name)

        # Set local flag (otherwise host is assigned automatically)
        local = host_id is None

        # Add the simulation to the queue
        self.launcher.add_to_queue(definition, name, host_id=host_id, parallelization=parallelization,
                                   logging_options=logging_options, scheduling_options=scheduling_options,
                                   analysis_options=analysis_options, local=local)

        # Add entry to new simulations table
        self.new.add_simulation(name, ski_path, host_id=host_id)

    # -----------------------------------------------------------------

    def cache(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Caching output of the simulations ...")

        # Loop over the analysed simulations
        for simulation in self.all_analysed_simulations:

            # Cache this simulation
            self.cache_simulation(simulation.name)

    # -----------------------------------------------------------------

    def cache_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name: 
        :return: 
        """

        # Debugging
        log.debug("Caching simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Cache output
        if self.config.cache_output:
            if self.is_cached_output(simulation.name): log.warning("Simulation output is already cached for simulation '" + simulation.name + "': skipping ...")
            else: self.cache_simulation_output(simulation.name)

        # Cache datacubes
        if self.config.cache_datacubes:
            if self.is_cached_output(simulation.name): log.warning("Simulation output is already cached for simulation '" + simulation.name + "': skipping ...")
            else: self.cache_simulation_datacubes(simulation.name)

        # Cache misc output
        if self.config.cache_misc and simulation.analysed_all_misc:
            if self.is_cached_misc(simulation.name): log.warning("Simulation miscellaneous output is already cached for simulation '" + simulation.name + "': skipping ...")
            else: self.cache_simulation_misc(simulation.name)

        # Cache images
        if self.config.cache_images and simulation.analysed_all_misc:
            if self.is_cached_misc(simulation.name): log.warning("Simulation miscellaneous output is already cached for simulation '" + simulation.name + "': skipping ...")
            else: self.cache_simulation_images(simulation.name)

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
