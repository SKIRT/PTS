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
from ..basics.configurable import Configurable
from ..basics.configuration import ConfigurationDefinition, parse_arguments
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
from ..basics.containers import create_nested_defaultdict
from ..tools import sequences
from ..basics.log import no_debugging
from ..simulation.remote import is_finished_status, is_running_status
from ..basics.configuration import prompt_string, prompt_variable
from ..remote.host import find_host_ids
from ..simulation.shower import show_simulation, show_analysis, compare_simulations, compare_analysis
from ..simulation.adapter import adapt_simulation, adapt_analysis
from ..config.show_simulation_settings import definition as show_simulation_definition
from ..config.adapt_simulation_settings import definition as adapt_simulations_definition
from ..config.show_analysis_options import definition as show_analysis_definition
from .batchlauncher import BatchLauncher
from ..basics.table import SmartTable
from ..tools import tables
from ..tools import time
from ..config.analyse_simulation import definition as analyse_simulation_definition
from .analyser import all_steps

# -----------------------------------------------------------------

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

plot_x_limits = (0., None)

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

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Run commands
        if self.do_commands: self.run_commands()

        # Interactive
        if self.do_interactive: self.interactive()

        # Move simulations
        if self.do_moving: self.move()

        # Launch
        if self.do_launch: self.launch()

        # Show
        if self.do_showing: self.show()

        # 3. Write
        if self.do_writing: self.write()

        # Plotting
        if self.do_plotting: self.plot()

        # Re-analyse
        if self.do_reanalysis: self.reanalyse()

        # Analyse
        if self.do_analysis: self.analyse()

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

        # Backup
        #if self.config.backup_simulations:
        #if self.config.backup_assignment:

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
        elif self.config.backup_dir_path is not None: return fs.create_directory_in(self.config.backup_dir_path, time.unique_name("manage"))
        else: return self.output_path_directory(time.unique_name("manage"), create=True)

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

        # Get the simulation
        return self.simulations[host_id][simulation_name]

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

                    # Get cluster name
                    cluster_name = simulation.cluster_name

                    # Change the host for this simulation
                    self.assignment.set_host_for_simulation(simulation_name, actual_host_id, cluster_name=cluster_name)
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

            # Not yet retrieved
            else:
                host_id = simulation.host_id
                screen_states = self.screen_states[host_id] if host_id in self.screen_states else None
                jobs_status = self.jobs_status[host_id] if host_id in self.jobs_status else None
                with no_debugging(): simulation_status = self.get_remote(host_id).get_simulation_status(simulation, screen_states=screen_states, jobs_status=jobs_status)

            # Add the status
            status_list.append(simulation_status)

        # Create the table and return
        return SimulationStatusTable.from_columns(self.simulation_names, status_list)

    # -----------------------------------------------------------------

    def reset_status(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Restting the simulation status ...")

        # Reset the jobs status
        del self.jobs_status

        # Reset the screen states
        del self.screen_states

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

    @lazyproperty
    def jobs_status(self):

        """
        This function gets the status of jobs
        :return:
        """

        # Initialize dictionary
        jobs_status_hosts = dict()

        # Loop over the hosts
        for host_id in self.simulations:
            if not self.get_remote(host_id).scheduler: continue

            # Get the status of the jobs
            jobs_status_hosts[host_id] = self.get_remote(host_id).get_jobs_status()

        # Return
        return jobs_status_hosts

    # -----------------------------------------------------------------

    @lazyproperty
    def screen_states(self):

        """
        This function gets the status of screen sessions
        :return:
        """

        # Initialize dictionary
        screen_states_hosts = dict()

        # Loop over the hosts
        for host_id in self.simulations:
            if self.get_remote(host_id).scheduler: continue

            # Get the status of the screen sessions
            screen_states_hosts[host_id] = self.get_remote(host_id).screen_states()

        # Return
        return screen_states_hosts

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

            # Process command
            success = True
            try: self.process_command(command)
            except Exception as e:
                log.warning(str(e))
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

    def process_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Show simulation status
        if command == "status": self.show_status()

        # Show things
        elif command.startswith("show"): self.show_command(command)

        # Plot things
        elif command.startswith("plot"): self.plot_command(command)

        # Move simulation
        elif command.startswith("move"): self.move_simulation_command(command)

        # Stop simulation?
        elif command.startswith("stop"): self.stop_simulation_command(command)

        # Relaunch simulation
        elif command.startswith("relaunch"): self.relaunch_simulation_command(command)

        # Show log of simulation
        elif command.startswith("log"): self.show_simulation_log_command(command)

        # Show simulation settings
        elif command.startswith("settings"): self.show_simulation_settings_command(command)

        # Show analysis options
        elif command.startswith("analysis"): self.show_analysis_options_command(command)

        # Adapt simulation settings or analysis options
        elif command.startswith("adapt"): self.adapt_simulation_command(command)

        # Compare simulation settings or analysis options
        elif command.startswith("compare"): self.compare_simulations_command(command)

        # Analyse simulation
        elif command.startswith("analyse"): self.analyse_simulation_command(command)

        # Re-analyse simulation
        elif command.startswith("reanalyse"): self.reanalyse_simulation_command(command)

        # Invalid command
        else: raise ValueError("Invalid command: '" + command + "'")

    # -----------------------------------------------------------------

    def show_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)

        # Assignment
        if splitted[1] == "assignment": self.show_assignment()

        # Status
        elif splitted[1] == "status": self.show_status()

        # Runtimes
        elif splitted[1] == "runtimes": self.show_runtimes_command(command)

        # Memory
        elif splitted[1] == "memory": self.show_memory_command(command)

        # Invalid
        else: raise ValueError("Invalid command")

    # -----------------------------------------------------------------

    def plot_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)

        # Runtimes
        if splitted[1] == "runtimes": self.plot_runtimes_command(command)

        # Memory
        elif splitted[1] == "memory": self.plot_memory_command(command)

        # Invalid
        else: raise ValueError("Invalid command")

    # -----------------------------------------------------------------

    def get_host_and_parallelization_from_command(self, command, name=None, index=1):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_positional_optional("host", "host", "remote host")
        definition.add_positional_optional("parallelization", "parallelization", "parallelization scheme")
        parse_command = splitted[index:]

        # Parse arguments
        config = parse_arguments(name, definition, command=parse_command, error="exception")

        # Return
        return config.host, config.parallelization

    # -----------------------------------------------------------------

    def parse_simulation_command(self, command, command_definition=None, name=None, index=1):

        """
        This function ...
        :param command:
        :param name:
        :param index:
        :param command_definition:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("simulation", "integer_or_string", "simulation index or name")
        if command_definition is not None:
            definition.import_settings(command_definition, required_to="optional")
            parse_command = splitted[index:]
        else: parse_command = splitted[index:index+1] # only simulation name

        # Parse arguments
        config = parse_arguments(name, definition, command=parse_command, error="exception")

        # Get simulation name
        if types.is_integer_type(config.simulation): simulation_name = self.simulation_names[config.pop("simulation")]
        else: simulation_name = config.pop("simulation")

        # Return
        return splitted, simulation_name, config

    # -----------------------------------------------------------------

    def parse_two_simulations_command(self, command, command_definition=None, name=None, index=1):

        """
        This function ....
        :param command:
        :param command_definition:
        :param name:
        :param index:
        :return:
        """

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("simulation_a", "integer_or_string", "simulation index or name")
        definition.add_required("simulation_b", "integer_or_string", "simulation index or name")
        definition.import_settings(command_definition, required_to="optional")

        # Parse arguments
        config = parse_arguments(name, definition, command=splitted[index:], error="exception")

        # Get simulation_a name
        if types.is_integer_type(config.simulation_a): simulation_a_name = self.simulation_names[config.pop("simulation_a")]
        else: simulation_a_name = config.pop("simulation_a")

        # Get simulation_b name
        if types.is_integer_type(config.simulation_b): simulation_b_name = self.simulation_names[config.pop("simulation_b")]
        else: simulation_b_name = config.pop("simulation_b")

        # Return
        return splitted, simulation_a_name, simulation_b_name, config

    # -----------------------------------------------------------------

    def get_simulation_name_from_command(self, command, name=None):

        """
        This function ...
        :param command:
        :param name:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, name=name)

        # Return the simulation name
        return simulation_name

    # -----------------------------------------------------------------

    def get_simulation_name_and_config_from_command(self, command, command_definition, name=None):

        """
        This function ...
        :param command:
        :param command_definition:
        :param name:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, command_definition=command_definition, name=name)

        # Return simulation name and config
        return simulation_name, config

    # -----------------------------------------------------------------

    def show_simulation_log_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command)

        # Check
        if not self.is_running_or_finished(simulation_name): raise ValueError("Simulation is not finished or running")

        # Show
        self.show_simulation_log(simulation_name)

    # -----------------------------------------------------------------

    def show_simulation_log(self, simulation_name):

        """
        This function ...
        :param simulation_name:
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
        for line in lines: print(line)
        print("")

    # -----------------------------------------------------------------

    def show_simulation_settings_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, show_simulation_definition, "show_simulation_settings")

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

    def show_analysis_options_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, show_analysis_definition, "show_analysis_options")

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

    def adapt_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Parse the command
        splitted, simulation_name, config = self.parse_simulation_command(command, command_definition=adapt_simulations_definition, name="adapt_simulation", index=2)

        # Simulation
        if splitted[1] == "simulation": self.adapt_simulation_settings(simulation_name, config)

        # Analysis
        elif splitted[1] == "analysis": self.adapt_analysis_options(simulation_name, config)

        # Invalid
        else: raise ValueError("Invalid argument")

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

    def compare_simulations_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Parse
        splitted, simulation_a_name, simulation_b_name, config = self.parse_two_simulations_command(command, show_simulation_definition, name="compare_simulations", index=2)

        # Compare simulation
        if splitted[1] == "simulation": self.compare_simulation_settings(simulation_a_name, simulation_b_name, config)

        # Compare analysis
        elif splitted[1] == "analysis": self.compare_analysis_options(simulation_a_name, simulation_b_name, config)

        # Invalid
        else: raise ValueError("Invalid argument")

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

    def move_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get the simulation name
        splitted, simulation_name, config = self.parse_simulation_command(command, name="move_simulation")

        # Check status
        if self.is_finished(simulation_name): raise ValueError("Simulation is already finished")
        if self.is_running(simulation_name): log.warning("Simulation is already running")

        # Get host and check
        host = load_host(splitted[2])
        if host == self.get_simulation(simulation_name).host: raise ValueError("Simulation '" + simulation_name + "' is already queued/running on host '" + tostr(host) + "'")

        # Move simulation
        self.move_simulation(simulation_name, host)

    # -----------------------------------------------------------------

    def move_simulation(self, simulation_name, host):

        """
        This function ...
        :param simulation_name:
        :param host:
        :return:
        """

        # Debugging
        log.debug("Moving simulation '" + simulation_name + "' to host '" + tostr(host) + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Add entry to moved table
        self.moved.add_simulation(simulation, host)

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

    def stop_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command)

        # Check whether queued or running
        if not self.is_queued_or_running(simulation_name): raise ValueError("The simulation '" + simulation_name + "' is not queued or running")

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

    # -----------------------------------------------------------------

    def relaunch_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name
        simulation_name = self.get_simulation_name_from_command(command)

        # Check simulation
        if not self.is_failed(simulation_name): raise ValueError("Simulation '" + simulation_name + "' is running, finished or still queued")

        # # Relaunch finished simulations?
        # if config.finished:
        #
        #     # Remove local output, if present
        #     if simulation.retrieved:
        #         fs.clear_directory(simulation.output_path)
        #         simulation.set_retrieved(False)
        #         simulation.save()
        #     if simulation.analysed_any_extraction:
        #         fs.clear_directory(simulation.analysis.extraction.path)
        #         simulation.unset_analysed_extraction()
        #         simulation.save()
        #     if simulation.analysed_any_plotting:
        #         fs.clear_directory(simulation.analysis.plotting.path)
        #         simulation.unset_analysed_plotting()
        #         simulation.save()
        #     if simulation.analysed_any_misc:
        #         fs.clear_directory(simulation.analysis.misc.path)
        #         simulation.unset_analysed_misc()
        #         simulation.save()

        # Relaunch the simulation
        self.relaunch_simulation(simulation_name)

    # -----------------------------------------------------------------

    def relaunch_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Debugging
        log.debug("Relaunching simulation '" + simulation_name + "' ...")

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

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching simulations ...")

        # Run the launcher
        self.launcher.run()

        # Set moved IDs
        for simulation in self.launcher.launched_simulations:
            if not self.is_moved(simulation.name): continue
            self.moved.set_new_id(simulation.name, simulation.id)

        # Reset status?
        if self.has_moved or self.has_relaunched: self.reset_status()

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

    def show_assignment(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the assignment scheme ...")

        # Show
        print(self.assignment)

    # -----------------------------------------------------------------

    @memoize_method
    def get_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return SKIRTRemote(host_id=host_id)

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

    def show_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the simulation status ...")

        # Show
        print("")
        print(fmt.bold + "Total number of simulations: " + fmt.reset + str(self.nsimulations))
        print(fmt.bold + "Number of finished simulations: " + fmt.reset + str(self.nfinished) + " (" + tostr(self.percentage_nfinished, round=True, ndigits=2) + "%)")
        print(fmt.bold + "Number of retrieved simulations: " + fmt.reset + str(self.nretrieved) + " (" + tostr(self.percentage_nretrieved, round=True, ndigits=2) + "%)")
        print(fmt.bold + "Number of analysed simulations: " + fmt.reset + str(self.nanalysed) + " (" + tostr(self.percentage_nanalysed, round=True, ndigits=2) + "%)")
        print("")

        # Loop over the simulations
        for index, simulation_name in enumerate(self.simulation_names):

            # Get the simulation
            simulation = self.get_simulation(simulation_name)

            # Get display name
            display_name = self.get_display_name(simulation, id_size=3) #host_id_size=10)

            # Get the status
            status = self.get_status(simulation_name)

            # Get index string
            index_string = "[" + strings.integer(index, 3, fill=" ") + "] "

            # Show status
            if status == "analysed": print(" - " + index_string + fmt.green + display_name + ": analysed" + fmt.reset)
            elif status == "retrieved": print(" - "  + index_string + fmt.yellow + display_name + ": retrieved" + fmt.reset)
            elif status == "finished": print(" - " + index_string + fmt.yellow + display_name + ": finished" + fmt.reset)
            elif "running" in status: print(" - " + index_string + display_name + ": " + status + fmt.reset)
            else: print(" - " + index_string + fmt.red + display_name + ": " + status + fmt.reset)

        # End with empty line
        print("")

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

    def show_runtimes_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, name="show_runtimes")

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

    def show_memory_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, name="show_memory")

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

        # Write the assignment
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

    def plot_runtimes_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, name="plot_runtimes")

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

    def plot_memory_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get host and parallelization
        host, parallelization = self.get_host_and_parallelization_from_command(command, name="plot_memory")

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

    def reanalyse_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
        definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")
        definition.add_optional("not_steps", "string_list", "don't analyse these steps", choices=all_steps)
        definition.add_optional("not_features", "string_list", "don't analyse these features (if a single not_step is defined)")

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, command_definition=definition)
        steps = config.steps
        features = config.features
        not_steps = config.not_steps
        not_features = config.not_features

        # Reanalyse the simulation
        self.reanalyse_simulation(simulation_name, steps=steps, features=features, not_steps=not_steps, not_features=not_features)

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

    def analyse_simulation_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Get simulation name and config
        simulation_name, config = self.get_simulation_name_and_config_from_command(command, command_definition=analyse_simulation_definition)

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
