#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.batchlauncher Contains the BatchLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import traceback
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.remote import SKIRTRemote
from ..remote.remote import Remote
from .options import LoggingOptions
from ..tools import introspection, time
from ..tools import filesystem as fs
from ..basics.log import log
from ..remote.host import load_host
from .analyser import SimulationAnalyser
from .options import AnalysisOptions
from ..simulation.definition import create_definitions
from ..basics.handle import ExecutionHandle
from ..simulation.execute import SkirtExec
from ..simulation.input import SimulationInput
from ..simulation.skifile import SkiFile
from ..simulation.arguments import SkirtArguments
from ..tools import formatting as fmt
from ..simulation.parallelization import Parallelization, get_possible_nprocesses_in_memory
from ..tools import monitoring
from ..advanced.memoryestimator import estimate_memory
from ..tools import sequences
from ..tools.utils import lazyproperty
from ..basics.map import Map
from ..advanced.parallelizationtool import determine_parallelization
from ..tools import parallelization
from .timing import TimingTable
from .memory import MemoryTable
from ..tools.stringify import tostr
from ..tools import numbers
from ..tools import types
from ..simulation.simulation import SkirtSimulation
from .tables import SimulationAssignmentTable, QueuedSimulationsTable

# -----------------------------------------------------------------

class DifferentMemoryUsages(Exception):

    """
    This class ...
    """

    def __init__(self, message):

        """
        Thisf unction ...
        :param message:
        """

        # Call the base class constructor with the parameters it needs
        super(DifferentMemoryUsages, self).__init__(message)

# -----------------------------------------------------------------

class MissingSimulation(Exception):

    """
    This class ...
    """

    def __init__(self, name, missing_from=None):

        """
        Thisf unction ...
        :param name:
        :param missing_from:
        :param:
        """

        # Set message
        message = "The simulation '" + name + "' is missing"
        if missing_from is not None: message += " from " + missing_from

        # Call the base class constructor with the parameters it needs
        super(MissingSimulation, self).__init__(message)

# -----------------------------------------------------------------

class DifferentNwavelengths(Exception):

    """
    This class ...
    """

    def __init__(self, message):

        """
        This function ...
        :param message:
        """

        # Call the base class constructor
        super(DifferentNwavelengths, self).__init__(message)

# -----------------------------------------------------------------

class DifferentDustLibDimensions(Exception):

    """
    This function ...
    """

    def __init__(self, message):

        """
        This function ...
        :param message:
        """

        # Call the base class constructor
        super(DifferentDustLibDimensions, self).__init__(message)

# -----------------------------------------------------------------

class DifferentNDustCells(Exception):

    """
    Thisf unction ...
    """

    def __init__(self, message):

        """
        This function ...
        :param message:
        """

        # Call the base class constructor
        super(DifferentNDustCells, self).__init__(message)

# -----------------------------------------------------------------

default_extraction_dirname = "extr"
default_plotting_dirname = "plot"
default_misc_dirname = "misc"

# -----------------------------------------------------------------

class BatchLauncher(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(BatchLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The local SKIRT execution context
        self.skirt = SkirtExec()

        # Initialize a list to contain different SKIRTRemote instances for the different remote hosts
        self.remotes = []

        # For some remote hosts, if they use a scheduling system, defines the name of the cluster to be used
        # (the default cluster will be used if not defined)
        self.cluster_names = dict()

        # The queue for locally executed simulations
        self.local_queue = []

        # The queues for the different remote hosts
        self.queues = defaultdict(list)

        # The scheduling options for (some of) the simulations and (some of) the remote hosts. This is a nested
        # dictionary where the first key represents the remote host ID and the next key represents the name of the
        # simulation (see 'name' parameter of 'add_to_queue')
        self.scheduling_options = defaultdict(dict)

        # The assignment from items in the queue to the different remote hosts
        self.assignment = None

        # The desired number of processes
        self.nprocesses_local = None

        # Flag for data parallelization
        self.data_parallel_local = None
        self.data_parallel_hosts = dict()

        # The memory usage (if defined, assumed to be the same for each simulation of the batch)
        self.memory = None

        # The number of dust cells (if defined, assumed to be the same for each simulation of the batch)
        self.ncells = None

        # The number of wavelengths (if defined, assumed to be the same for each simulation of the batch)
        self.nwavelengths = None

        # The parallelization scheme for local execution
        self.parallelization_local = None

        # The parallelization scheme for the different remote hosts
        self.parallelization_hosts = dict()

        # The number of processes for the different remote hosts
        self.nprocesses_hosts = dict()
        self.nprocesses_per_node_hosts = dict()

        # The parallelization scheme for the different simulations
        self.parallelization_simulations = dict()

        # The logging options for the different simulations
        self.logging_simulations = dict()

        # The paths to the directories for placing the (screen/job) scripts (for manual inspection) for the different remote hosts
        self.script_paths = dict()

        # The list of remote hosts for which the screen output should be saved remotely (for debugging)
        self.save_screen_output = []

        # The simulations that have been launched
        self.launched_simulations = []

        # The simulations that have been retrieved
        self.simulations = []

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # The logging options
        self.logging_options = None

        # Create temporary local directory
        self.temp_path = introspection.create_temp_dir(time.unique_name("BatchLauncher"))

        # Original definitions of local simulations (when definitions are changed because of shared input)
        self.original_local_definitions = dict()

        # The simulation assignment table
        self.assignment = None

        # Groups of simulations with shared input
        self.shared_input_groups = []

        # Remote input paths
        self.remote_input_paths_hosts = dict()

        # The timing and memory table
        self.timing_table = None
        self.memory_table = None

    # -----------------------------------------------------------------

    @property
    def uses_remotes(self):

        """
        This function ...
        :return:
        """

        return len(self.config.remotes) > 0 or len(self.remotes) > 0

    # -----------------------------------------------------------------

    @property
    def only_local(self):

        """
        This function ...
        :return:
        """

        return len(self.config.remotes) == 0 and len(self.remotes) == 0

    # -----------------------------------------------------------------

    def in_queue_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return len(self.queues[host_id]) if host_id in self.queues else 0

    # -----------------------------------------------------------------

    @property
    def shortest_queue_host_id(self):

        """
        This function ...
        :return:
        """

        size = float('inf')
        host_id = None

        for host in self.hosts:
            queue_size = self.in_queue_for_host(host.id)
            if queue_size < size:
                size = queue_size
                host_id = host.id

        # Return the ID of the shortest
        return host_id

    # -----------------------------------------------------------------

    def add_to_queue(self, definition, name, host_id=None, parallelization=None, logging_options=None,
                     scheduling_options=None, analysis_options=None, local=False):

        """
        This function ...
        :param definition:
        :param name: a name that is given to the simulation
        :param host_id: the host on which this simulation should be executed
        :param parallelization: individual parallelization scheme for this particular simulation (only allowed if host_id is also specified)
        :param logging_options:
        :param scheduling_options:
        :param analysis_options: analysis options (if None, analysis options will be created from batch launcher configuration)
        :param local:
        :return:
        """

        # Debugging
        log.debug("Adding a simulation '" + name + "' to the queue with:")
        log.debug("")
        log.debug(" - input: " + tostr(definition.input_path))
        log.debug(" - ski path: " + definition.ski_path)
        log.debug(" - output path: " + definition.output_path)
        log.debug("")

        # Check whether the simulation name doesn't contain spaces
        if " " in name: raise ValueError("The simulation name cannot contain spaces")

        # Checks
        if (local or self.has_no_remotes) and host_id is not None: raise ValueError("Host ID cannot be specified for local execution")
        if (local or self.has_no_remotes) and scheduling_options is not None: raise ValueError("Scheduling options cannot be specified for local execution")

        # Local execution
        if local or self.has_no_remotes: self.add_to_local_queue(definition, name, parallelization, logging_options, analysis_options)

        # Remote execution
        else: self.add_to_remote_queue(definition, name, host_id, parallelization, logging_options, scheduling_options, analysis_options)

    # -----------------------------------------------------------------

    def add_to_local_queue(self, definition, name, parallelization=None, logging_options=None, analysis_options=None):

        """
        This function ...
        :param definition:
        :param name:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :param analysis_options:
        :return:
        """

        # Add to the local queue
        self.local_queue.append((definition, name, analysis_options))

        # Set parallelization, logging options
        if parallelization is not None: self.set_parallelization_for_simulation(name, parallelization)
        if logging_options is not None: self.set_logging_options_for_simulation(name, logging_options)

    # -----------------------------------------------------------------

    def add_to_remote_queue(self, definition, name, host_id=None, parallelization=None, logging_options=None,
                            scheduling_options=None, analysis_options=None):

        """
        This function ...
        :param definition:
        :param name:
        :param host_id:
        :param parallelization:
        :param logging_options:
        :param scheduling_options:
        :param analysis_options:
        :return:
        """

        # Checks
        # Check whether, if parallelization is specified, that host ID is also specified
        if host_id is None and parallelization is not None: raise ValueError("If parallelization is specified, host ID must also be specified")
        if host_id is None and scheduling_options is not None: raise ValueError("If scheduling options are specified, host ID must also be specified")

        # Check if this is a valid host ID
        if host_id is not None and host_id not in self.host_ids: raise ValueError("Invalid host ID")

        # Set host ID if undefined
        # Determine the ID of the hosts with the shortest queue (or the first if the queues are equally long)
        if host_id is None: host_id = self.shortest_queue_host_id

        # Add to the queue of the specified host
        self.queues[host_id].append((definition, name, analysis_options))

        # Set parallelization, logging_options, scheduling options
        if parallelization is not None: self.set_parallelization_for_simulation(name, parallelization)
        if logging_options is not None: self.set_logging_options_for_simulation(name, logging_options)
        if scheduling_options is not None: self.set_scheduling_options(host_id, name, scheduling_options)

    # -----------------------------------------------------------------

    def get_remote(self, host_id):

        """
        This function ..
        :param host_id:
        :return:
        """

        if isinstance(host_id, Remote): return host_id

        # Search in remotes
        for remote in self.remotes:
            if remote.host_id == host_id: return remote

        # Create new remote (shouldn't happen if the setup has been called)
        return Remote(host_id=host_id)

    # -----------------------------------------------------------------

    def get_skirt_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        if isinstance(host_id, SKIRTRemote): return host_id

        # Search in remotes
        for remote in self.remotes:
            if remote.host_id == host_id:
                assert isinstance(remote, SKIRTRemote)
                return remote

        # Create new remote (shouldn't happen if the setup has been called)
        return SKIRTRemote(host_id=host_id)

    # -----------------------------------------------------------------

    def set_cluster_for_host(self, host_id, cluster_name):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :return:
        """

        self.cluster_names[host_id] = cluster_name

    # -----------------------------------------------------------------

    def get_clustername_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Cluster name is specifically defined
        if host_id in self.cluster_names and self.cluster_names[host_id] is not None: return self.cluster_names[host_id]

        # The setup has not been called yet
        elif len(self.remotes) == 0:

            # Get the default cluster name
            host = load_host(host_id)

            # Return the default cluster name
            return host.clusters.default

        # The setup has been called
        else: return self.get_remote(host_id).cluster_name

    # -----------------------------------------------------------------

    def set_scheduling_options(self, host_id, name, options):

        """
        This function ...
        :param name:
        :param host_id:
        :param options:
        :return:
        """

        # Set the scheduling options for the specified host and for the specified simulation
        self.scheduling_options[host_id][name] = options

    # -----------------------------------------------------------------

    @property
    def in_queue(self):

        """
        This function ...
        :return:
        """

        # Return the total number of simulations in the queues
        return self.in_local_queue + self.in_remote_queues

    # -----------------------------------------------------------------

    @property
    def in_local_queue(self):

        """
        This function ...
        :return:
        """

        return len(self.local_queue)

    # -----------------------------------------------------------------

    @property
    def in_remote_queues(self):

        """
        This function ...
        :return:
        """

        total = 0
        for host_id in self.host_ids: total += self.in_remote_queue(host_id)
        return total

    # -----------------------------------------------------------------

    @property
    def has_queued_local(self):

        """
        This function ...
        :return:
        """

        return self.in_local_queue > 0

    # -----------------------------------------------------------------

    @property
    def has_queued_remotes(self):

        """
        This function ...
        :return:
        """

        return self.in_remote_queues > 0

    # -----------------------------------------------------------------

    @property
    def has_queued(self):

        """
        This function ...
        :return:
        """

        return self.has_queued_local or self.has_queued_remotes

    # -----------------------------------------------------------------

    def has_queued_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.in_remote_queue(host_id) > 0

    # -----------------------------------------------------------------

    def in_remote_queue(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.in_queue_for_host(host_id)

    # -----------------------------------------------------------------

    @property
    def single_remote(self):

        """
        This function ...
        :return:
        """

        # The setup has not been called yet
        if len(self.remotes) == 0:

            # Warning
            log.warning("The setup has not been called yet, creating temporary remote execution environment (limited to basic capabilities)")

            # Get the single host ID
            host_id = self.single_host_id

            # Create a remote execution context
            remote = Remote()

            # Check whether a cluster name is defined
            cluster_name = self.cluster_names[host_id] if host_id in self.cluster_names else None

            # Setup the remote for the specified host
            remote.setup(host_id, cluster_name=cluster_name)

            # Return the remote environment
            return remote

        elif len(self.remotes) == 1: return self.remotes[0]
        else: raise RuntimeError("Multiple remotes have been configured for this batch launcher")

    # -----------------------------------------------------------------

    @property
    def single_host(self):

        """
        This function ...
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            if len(self.host_ids) != 1: raise RuntimeError("Multiple remotes have been configured for this batch launcher")
            else: return self.hosts[0]

        else: return self.single_remote.host

    # -----------------------------------------------------------------

    @property
    def single_host_id(self):

        """
        This function ...
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            if len(self.host_ids) != 1: raise RuntimeError("Multiple remotes have been configured for this batch launcher")
            else: return self.host_ids[0]

        else: return self.single_remote.host_id

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function returns the IDs of the hosts that will be used by the batch launcher
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # If a list of remotes is defined
            if self.config.remotes is not None: host_ids = self.config.remotes

            # If a list of remotes is not defined, create a remote for all of the hosts that have a configuration file
            else: host_ids = introspection.remote_host_ids()

            # Return the list of host IDs
            return host_ids

        # If the setup has already been called
        else: return [remote.host_id for remote in self.remotes]

    # -----------------------------------------------------------------

    @property
    def queue_host_ids(self):

        """
        This function ...
        :return:
        """

        return [host_id for host_id in self.queues.keys() if self.has_queued_remote(host_id)]

    # -----------------------------------------------------------------

    @property
    def hosts(self):

        """
        This function returns the Host objects for all the hosts that will be used by the batch launcher
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Create the list of hosts
            hosts = []
            #self.get_clustername_for_host(host_id)
            for host_id in self.host_ids: hosts.append(load_host(host_id))

            # Return the list of hosts
            return hosts

        # If the setup has already been called
        else: return [remote.host for remote in self.remotes]

    # -----------------------------------------------------------------

    @property
    def scheduler_host_ids(self):

        """
        This function returns just the IDS of the hosts which use a scheduling system. This is useful for when example
        the user wants to specify the parallelization scheme only for these hosts, while the parallelization strategy
        for the other hosts can be left up to the BatchLauncher class based on the current load of that system.
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Initialize a list to contain the host IDs
            host_ids = []

            # Loop over the IDs of all the hosts used by the BatchLauncher
            for host_id in self.host_ids:

                # Create Host instance
                host = load_host(host_id)

                # If it's a scheduler, add it to the list
                if host.scheduler: host_ids.append(host_id)

            # Return the list of hosts which use a scheduling system
            return host_ids

        # If the setup has already been called
        else: return [remote.host_id for remote in self.remotes if remote.scheduler]

    # -----------------------------------------------------------------

    @property
    def scheduler_hosts(self):

        """
        This function is similar to 'scheduler_host_ids' but returns the Host objects instead of just the IDs of these
        hosts.
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Initialize a list to contain the hosts
            hosts = []

            # Loop over the IDs of all the hosts used by the BatchLauncher
            for host_id in self.host_ids:

                # Create a Host instance
                host = load_host(host_id)

                # If it's a scheulder, add it to the list
                if host.scheduler: hosts.append(host)

            # Return the list of hosts
            return hosts

        # If the setup has already been called
        else: return [remote.host for remote in self.remotes if remote.scheduler]

    # -----------------------------------------------------------------

    @property
    def no_scheduler_host_ids(self):

        """
        This function returns just the IDS of the hosts which DON'T use a scheduling system.
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Initialize a list to contain the host IDs
            host_ids = []

            # Loop over the IDs of all the hosts used by the BatchLauncher
            for host_id in self.host_ids:

                # Create Host instance
                host = load_host(host_id)

                # If it's not a scheduler, add it to the list
                if not host.scheduler: host_ids.append(host_id)

            # Return the list of hosts which don't use a scheduling system
            return host_ids

        # If the setup has already been called
        else: return [remote.host_id for remote in self.remotes if not remote.scheduler]

    # -----------------------------------------------------------------

    @property
    def no_scheduler_hosts(self):

        """
        This function is similar to 'no_scheduler_host_ids' but returns the Host objects instead of just the IDs of these
        hosts.
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Initialize a list to contain the hosts
            hosts = []

            # Loop over the IDs of all the hosts used by the BatchLauncher
            for host_id in self.host_ids:

                # Create a Host instance
                host = load_host(host_id)

                # If it's a not scheulder, add it to the list
                if not host.scheduler: hosts.append(host)

            # Return the list of hosts
            return hosts

        # If the setup has already been called
        else: return [remote.host for remote in self.remotes if not remote.scheduler]

    # -----------------------------------------------------------------

    @property
    def uses_schedulers(self):

        """
        This function ...
        :return:
        """

        return len(self.scheduler_host_ids) > 0

    # -----------------------------------------------------------------

    @property
    def nremotes(self):

        """
        This function ...
        :return:
        """

        #return len(self.remotes)
        if len(self.remotes) == 0: # perhaps setup not been called
            return len(self.config.remotes)
        else: return len(self.remotes)

    # -----------------------------------------------------------------

    @property
    def has_no_remotes(self):

        """
        Thisn function ...
        :return:
        """

        return self.nremotes == 0

    # -----------------------------------------------------------------

    @property
    def has_single_remote(self):

        """
        This function ...
        :return:
        """

        return self.nremotes == 1

    # -----------------------------------------------------------------

    def set_remote_input_path_for_host(self, host_id, path):

        """
        This function ...
        :param host_id:
        :param path:
        :return:
        """

        self.remote_input_paths_hosts[host_id] = path

    # -----------------------------------------------------------------

    def has_remote_input_path_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return host_id in self.remote_input_paths_hosts and self.remote_input_paths_hosts[host_id] is not None

    # -----------------------------------------------------------------

    def get_remote_input_path_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.remote_input_paths_hosts[host_id]

    # -----------------------------------------------------------------

    def set_logging_options_for_simulation(self, name, logging_options):

        """
        This function ...
        :param name:
        :param logging_options:
        :return:
        """

        # Set the logging options for this simulation
        self.logging_simulations[name] = logging_options

    # -----------------------------------------------------------------

    def has_logging_options_for_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.logging_simulations and self.logging_simulations[name] is not None

    # -----------------------------------------------------------------

    @property
    def has_logging_options_for_any_simulation(self):

        """
        This function ...
        :return:
        """

        for simulation_name in self.all_simulation_names:
            if self.has_logging_options_for_simulation(simulation_name): return True
        return False

    # -----------------------------------------------------------------

    def get_logging_options_for_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.logging_simulations[name] if name in self.logging_simulations else None

    # -----------------------------------------------------------------

    def set_parallelization_for_local(self, parallelization):

        """
        This function ...
        :param parallelization:
        :return:
        """

        # Set
        self.parallelization_local = parallelization

    # -----------------------------------------------------------------

    def set_parallelization_for_host(self, host_id, parallelization):

        """
        This function ...
        :param host_id:
        :param parallelization:
        :return:
        """

        # Set the parallelization properties for the specified host
        self.parallelization_hosts[host_id] = parallelization

    # -----------------------------------------------------------------

    def get_parallelization_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        if not self.has_parallelization_for_host(host_id): raise ValueError("Parallelization is not defined for remote host '" + host_id + "'")
        return self.parallelization_hosts[host_id]

    # -----------------------------------------------------------------

    def set_parallelization_for_simulation(self, name, parallelization):

        """
        This function ...
        :param name:
        :param parallelization:
        :return:
        """

        # Set the parallelization for this simulation
        self.parallelization_simulations[name] = parallelization

    # -----------------------------------------------------------------

    def set_nprocesses_local(self, nprocesses):

        """
        This function ...
        :param nprocesses:
        :return:
        """

        self.nprocesses_local = nprocesses

    # -----------------------------------------------------------------

    def set_nprocesses_for_host(self, host_id, nprocesses):

        """
        This function ...
        :param host_id:
        :param nprocesses:
        :return:
        """

        # Set the number of processes for the specified host
        self.nprocesses_hosts[host_id] = nprocesses

    # -----------------------------------------------------------------

    def set_nprocesses_per_node_for_host(self, host_id, nprocesses):

        """
        Thisn function ...
        :param host_id:
        :param nprocesses:
        :return:
        """

        # Set the number of processes per node for the specified host
        self.nprocesses_per_node_hosts[host_id] = nprocesses

    # -----------------------------------------------------------------

    def set_data_parallel_local(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.data_parallel_local = value

    # -----------------------------------------------------------------

    def set_data_parallel_for_host(self, host_id, value):

        """
        This function ...
        :param host_id:
        :param value:
        :return:
        """

        self.data_parallel_hosts[host_id] = value

    # -----------------------------------------------------------------

    def parallelization_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.parallelization_hosts[host_id] if host_id in self.parallelization_hosts else None

    # -----------------------------------------------------------------

    def has_parallelization_for_host(self, host_id):

        """
        This funtion ...
        :param host_id:
        :return:
        """

        return host_id in self.parallelization_hosts and self.parallelization_hosts[host_id] is not None

    # -----------------------------------------------------------------

    @property
    def has_parallelization_for_any_host(self):

        """
        This function ...
        :return:
        """

        for host_id in self.host_ids:
            if self.has_parallelization_for_host(host_id): return True
        return False

    # -----------------------------------------------------------------

    def has_parallelization_for_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.parallelization_simulations and self.parallelization_simulations[name] is not None

    # -----------------------------------------------------------------

    @property
    def has_parallelization_for_any_simulation(self):

        """
        This function ...
        :return:
        """

        for simulation_name in self.all_simulation_names:
            if self.has_parallelization_for_simulation(simulation_name): return True
        return False

    # -----------------------------------------------------------------

    def has_parallelization_for_all_simulations_of_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        for simulation_name in self.get_simulation_names_for_host(host_id):
            if not self.has_parallelization_for_simulation(simulation_name): return False
        return True

    # -----------------------------------------------------------------

    def nprocesses_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.nprocesses_hosts[host_id] if host_id in self.nprocesses_hosts else None

    # -----------------------------------------------------------------

    def has_nprocesses_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return host_id in self.nprocesses_hosts and self.nprocesses_hosts[host_id] is not None

    # -----------------------------------------------------------------

    def get_nprocesses_for_host(self, host_id):

        """
        Thisnfunction ...
        :param host_id:
        :return:
        """

        return self.nprocesses_hosts[host_id]

    # -----------------------------------------------------------------

    def nprocesses_per_node_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.nprocesses_per_node_hosts[host_id] if host_id in self.nprocesses_per_node_hosts else None

    # -----------------------------------------------------------------

    def has_nprocesses_per_node_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return host_id in self.nprocesses_per_node_hosts and self.nprocesses_per_node_hosts[host_id] is not None

    # -----------------------------------------------------------------

    def get_nprocesses_per_node_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.nprocesses_per_node_hosts[host_id]

    # -----------------------------------------------------------------

    def get_parallelization_for_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        #if not self.has_parallelization_for_simulation(name): raise ValueError("Parallelization is not defined for simulation '" + name + "'")
        return self.parallelization_simulations[name] if name in self.parallelization_simulations else None

    # -----------------------------------------------------------------

    @lazyproperty
    def runtime_estimator(self):

        """
        This function ...
        :return:
        """

        from ..advanced.runtimeestimator import RuntimeEstimator

        # Create a RuntimeEstimator instance
        return RuntimeEstimator(self.timing_table)

    # -----------------------------------------------------------------

    @property
    def has_timing_table(self):

        """
        This function ...
        :return:
        """

        return self.timing_table is not None

    # -----------------------------------------------------------------

    @property
    def has_memory_table(self):

        """
        This function ...
        :return:
        """

        return self.memory_table is not None

    # -----------------------------------------------------------------

    @property
    def testing(self):

        """
        This function ...
        :return:
        """

        return self.config.test

    # -----------------------------------------------------------------

    @property
    def do_check_input(self):

        """
        This function ...
        :return:
        """

        return self.config.shared_input is None

    # -----------------------------------------------------------------

    @property
    def do_estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        return self.uses_schedulers and self.has_timing_table

    # -----------------------------------------------------------------

    @property
    def do_show(self):

        """
        This function ...
        :return:
        """

        return self.config.show

    # -----------------------------------------------------------------

    @property
    def do_launch(self):

        """
        This function ...
        :return:
        """

        return not self.testing

    # -----------------------------------------------------------------

    @property
    def do_write(self):

        """
        This function ...
        :return:
        """

        return self.config.write

    # -----------------------------------------------------------------

    @property
    def do_retrieve(self):

        """
        This function ...
        :return:
        """

        return self.config.retrieve and not self.testing

    # -----------------------------------------------------------------

    @property
    def do_show_finished(self):

        """
        This function ...
        :return:
        """

        return self.has_simulations and self.config.show_finished

    # -----------------------------------------------------------------

    @property
    def do_analyse(self):

        """
        This function ...
        :return:
        """

        return self.config.analyse and self.has_simulations

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Check the input files for all simulations
        if self.do_check_input: self.check_input()

        # 3. Set the parallelization scheme for the remote hosts for which this was not specified by the user
        self.set_parallelization()

        # 4. Estimate the runtimes if necessary
        if self.do_estimate_runtimes: self.estimate_runtimes()

        # 5. Show
        if self.do_show: self.show()

        # 6. Launch the simulations
        if self.do_launch: self.launch()

        # 7. Write (before trying to receive and analyse)
        if self.do_write: self.write()

        # 8. Retrieve the simulations that are finished
        if self.do_retrieve: self.try_retrieving()

        # 9. Show the simulations that are finished
        if self.do_show_finished: self.show_finished()

        # 10. Analyse the output of the retrieved simulations
        if self.do_analyse: self.try_analysing()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Log out from the remotes
        for remote in self.remotes: remote.logout()

        # Clear the list of remotes
        self.remotes = []

        # Clear the cluster names
        self.cluster_names = dict()

        # Clear the queues
        self.local_queue = []
        self.queues = defaultdict(list)

        # Clear the scheduling options
        self.scheduling_options = defaultdict(dict())

        # Clear the assignment
        self.assignment = None

        # Clear the launched simulations
        self.launched_simulations = []

        # Clear the simulations
        self.simulations = []

        # Clear the script path dictionary
        self.script_paths = dict()

    # -----------------------------------------------------------------

    @property
    def has_remotes(self):

        """
        This function ...
        :return: 
        """

        return len(self.remotes) > 0

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchLauncher, self).setup(**kwargs)

        # Restrict remotes only to the remotes for which there are queued simulations
        self.config.remotes = self.queue_host_ids

        # Remotes are passed as instances?
        if "remotes" in kwargs:
            remotes = kwargs.pop("remotes")
            if types.is_dictionary(remotes): self.remotes = remotes.values()
            elif types.is_sequence(remotes): self.remotes = remotes
            else: raise ValueError("Invalid type for 'remotes'")

        # Setup the remote instances
        if not self.has_remotes: self.setup_remotes()

        # Create the logging options
        self.logging_options = LoggingOptions()
        self.logging_options.set_options(self.config.logging)

        # Add simulation definitions for ski file present in the current working directory (if nothing is added to the queue manually)
        if self.in_queue == 0: self.load_queue()

        # If attached mode is enabled, check whether all remotes are non-schedulers
        if self.config.attached and self.uses_schedulers: raise ValueError("Cannot use attached mode when remotes with scheduling system are used")

        # If dry mode is enabled, check whether all remotes are non-schedulers
        if self.config.dry and self.uses_schedulers: log.warning("Using dry mode when remotes with scheduling system are used can be impractical, since it requires uploading and submitting each seperate job script")

        # Get input

        ## Get the memory information passed to this instance
        self.memory = kwargs.pop("memory", None)

        ## Get the number of dust cells if given
        if "ncells" in kwargs: self.ncells = kwargs.pop("ncells")

        ## Get the number of wavelengths if given
        if "nwavelengths" in kwargs: self.nwavelengths = kwargs.pop("nwavelengths")

        ## Get the parallelization
        if "parallelization_local" in kwargs: self.parallelization_local = kwargs.pop("parallelization_local")

        ## Get the number of processes
        if "nprocesses_local" in kwargs: self.nprocesses_local = kwargs.pop("nprocesses_local")

        # Initialize the simulation assignment table
        self.assignment = SimulationAssignmentTable()

        # Get timing table
        self.get_timing_table(**kwargs)

        # Get memory table
        self.get_memory_table(**kwargs)

    # -----------------------------------------------------------------

    def get_timing_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "timing" in kwargs: self.timing_table = kwargs.pop("timing")
        elif self.config.timing_table_path is not None: self.timing_table = TimingTable.from_file(self.config.timing_table_path)

    # -----------------------------------------------------------------

    def get_memory_table(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        if "memory" in kwargs: self.memory_table = kwargs.pop("memory")
        elif self.config.memory_table_path is not None: self.memory_table = MemoryTable.from_file(self.config.memory_table_path)

    # -----------------------------------------------------------------

    def setup_remotes(self):

        """
        This function ...
        :return:
        """

        # If a list of remotes is defined
        if self.config.remotes is not None: host_ids = self.config.remotes
        else: host_ids = []

        # If a list of remotes is not defined, execute locally
        if len(host_ids) == 0: log.warning("No remote set: simulations will be performed locally")

        # Loop over all the remote host ids
        for host_id in host_ids:

            # Create a remote SKIRT execution context
            remote = SKIRTRemote()

            # Check whether a cluster name is defined
            cluster_name = self.cluster_names[host_id] if host_id in self.cluster_names else None

            # Setup the remote for the specified host, check if succesful
            if not remote.setup(host_id, cluster_name=cluster_name):
                log.warning("Remote host '" + host_id + "' is down")
                continue

            # Add the remote to the list of remote objects
            self.remotes.append(remote)

        # Check if any remotes are initialized
        if len(host_ids) > 0 and len(self.remotes) == 0: raise RuntimeError("None of the remotes are available")

    # -----------------------------------------------------------------

    def load_queue(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading ski files in the current working directory into the queue ...")

        # Create simulation definitions from the working directory and add them to the queue
        for definition in create_definitions(self.config.path, self.config.simulation_output, self.config.simulation_input, recursive=self.config.recursive):

            # Add the definition to the queue
            self.add_to_queue(definition, definition.prefix)

    # -----------------------------------------------------------------

    def get_definition_for_local_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        for definition, name, _ in self.local_queue:
            if name == simulation_name: return definition

        raise ValueError("No simulation found in the local queue with the name '" + simulation_name + "'")

    # -----------------------------------------------------------------

    def check_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the input of the simulations ...")

        # Local
        if self.in_local_queue: self.check_input_local()

        # Remote
        if self.in_remote_queues: self.check_input_remotes()

    # -----------------------------------------------------------------

    def check_input_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the input of the simulations that are run locally ...")

        # A dictionary with the simulation names for each input file that is shared
        shared_input = defaultdict(list)

        # Check which definitions in the local queue use the same input
        for definition, simulation_name, analysis_options in self.local_queue:

            # Skip simulations without input
            if definition.input_path is None: continue

            # Convert definition input to actual SimulationInput object
            simulation_input = SimulationInput.from_any(definition.input_path)

            # Loop over the input files
            for filename, filepath in simulation_input:

                # Set
                shared_input[filepath].append((simulation_name, filename))

            # Set the simulation input object
            definition.input_path = simulation_input

        # New paths for each shared input file
        new_paths = dict()

        # List of participating simulations
        participating_simulations = set()

        # Loop over the filepaths that are shared
        for filepath in shared_input:

            # Not shared
            if len(shared_input[filepath]) < 2: continue

            # Add all corresponding simulations to the list of participating simulations
            for simulation_name, _ in shared_input[filepath]: participating_simulations.add(simulation_name)

            # Generate a unique filename for this file
            original_filename = fs.name(filepath)
            extension = fs.get_extension(original_filename)
            bare_original_filename = fs.strip_extension(original_filename)
            new_filename = time.unique_name(bare_original_filename, precision="micro") + "." + extension

            # Copy the file to the temporary directory
            fs.copy_file(filepath, self.temp_path, new_name=new_filename)
            new_filepath = fs.join(self.temp_path, new_filename)

            # Set the new path
            new_paths[filepath] = new_filepath

        # Copy the rest of the input files of the participating simulations
        for simulation_name in participating_simulations:

            # Get the definition
            definition = self.get_definition_for_local_simulation(simulation_name)

            # Set the original definition
            self.original_local_definitions[simulation_name] = definition.copy()

            # Get the original simulation input specification
            original_simulation_input = definition.input_path

            # Set the input path to the temporary directory path
            definition.input_path = self.temp_path

            # Open the ski file
            ski = SkiFile(definition.ski_path)

            # List of shared filenames for this simulation
            shared_filenames_for_simulation = []

            # Loop over the shared input files, set the new file name
            for filepath in shared_input:

                ski_filename = None
                for sim_name, filename in shared_input[filepath]:
                    if sim_name == simulation_name:
                        ski_filename = filename
                        break
                else: continue

                # Get the new filepath
                new_filepath = new_paths[filepath]
                new_filename = fs.name(new_filepath)

                # Change in the ski file
                ski.change_input_filename(ski_filename, new_filename)

                # Add to list
                shared_filenames_for_simulation.append(new_filename)

            # Loop over other input files for this simulation, also copy them to the temporary directory with a unique name
            for filename in ski.input_files:

                # Already changed
                if filename in shared_filenames_for_simulation: continue

                # Get path for this filename
                filepath = original_simulation_input[filename]

                # Generate unique filename
                extension = fs.get_extension(filename)
                bare_filename = fs.strip_extension(filename)
                new_filename = time.unique_name(bare_filename, precision="micro") + "." + extension

                # Save to temporary directory
                fs.copy_file(filepath, self.temp_path, new_name=new_filename)

                # Change filename in ski
                ski.change_input_filename(filename, new_filename)

            # Save the ski file
            new_ski_path = fs.join(self.temp_path, time.unique_name(simulation_name, precision="micro") + ".ski")
            ski.saveto(new_ski_path)

            # Set the new ski path to the definition
            definition.ski_path = new_ski_path

    # -----------------------------------------------------------------

    def check_input_remotes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the input of the simulations that are run remotely ...")

        # Loop over the remotes
        for host_id in self.host_ids: self.check_input_for_remote(host_id)

    # -----------------------------------------------------------------

    def check_input_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # A dictionary with the simulation names for each input file that is shared
        shared_input = defaultdict(list)

        # Check which definitions in the queue of the remote host use the same input
        for definition, simulation_name, analysis_options_item in self.queues[host_id]:

            # Skip simulations without input
            if definition.input_path is None: continue

            # Convert definition input to actual SimulationInput object
            simulation_input = SimulationInput.from_any(definition.input_path)

            # Loop over the input files
            for filename, filepath in simulation_input:

                # Set
                shared_input[filepath].append((simulation_name, filename))

            # Set the simulation input object
            definition.input_path = simulation_input

        # List of participating simulations (simulations with a file that is shared with at least one other simulation)
        participating_simulations = set()

        # List of truly shared (at least 2 simulations) input files
        truly_shared_filepaths = []

        # Loop over the filepaths that are shared
        for filepath in shared_input:

            # Not shared
            if len(shared_input[filepath]) < 2: continue

            # Add all corresponding simulations to the list of participating simulations
            for simulation_name, _ in shared_input[filepath]: participating_simulations.add(simulation_name)

            # Add filepath
            truly_shared_filepaths.append(filepath)

        # Loop over the participating simulations
        for definition, simulation_name, analysis_options_item in self.queues[host_id]:
            if simulation_name not in participating_simulations: continue

            # If at least one input file is not shared with others, skip the simulation
            simulation_input = definition.input_path
            if not all_input_files_are_shared(simulation_input, truly_shared_filepaths): continue

            # For each input file, get the simulation names that share the file
            shared_simulations = defaultdict(list)
            for filename, filepath in simulation_input:
                for sim_name, sim_filename in shared_input[filepath]:
                    if sim_name == simulation_name: continue
                    shared_simulations[filepath].append(sim_name)

            # Check if there are other simulations that share all of the input files
            other_simulation_names = sequences.get_values_in_all(*shared_simulations.values())
            if len(other_simulation_names) == 0: continue

            # Set list of shared simulations
            shared_simulation_names = [simulation_name] + other_simulation_names

            # Set group of simulations with shared input
            self.set_shared_input_group(shared_simulation_names)

    # -----------------------------------------------------------------

    def set_shared_input_group(self, simulation_names):

        """
        This function ...
        :param simulation_names:
        :return:
        """

        self.shared_input_groups.append(simulation_names)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting parallelization schemes ...")

        # Set parallelization scheme for local execution
        if self.parallelization_local is None and self.has_queued_local: self.set_parallelization_local()
        elif self.has_queued_local and self.config.check_paralleliation: self.check_parallelization_local()

        # Set parallelization schemes for remote execution
        if self.has_remotes: self.set_parallelization_remote()

    # -----------------------------------------------------------------

    def check_parallelization_local(self):

        """
        Thisj function ...
        :return:
        """

        # Inform the user
        log.info("Checking the local parallelization scheme ...")

        # Determine the total number of requested threads
        requested_threads = self.parallelization_local.nthreads

        # Determine the total number of hardware threads that can be used on the remote host
        hardware_threads = parallelization.ncores()
        if parallelization.has_hyperthreading(): hardware_threads *= parallelization.nthreads_per_core()

        # If the number of requested threads is greater than the allowed number of hardware threads, raise
        # an error
        if requested_threads > hardware_threads: raise RuntimeError("The requested number of processes and threads "
                                                                    "exceeds the total number of hardware threads")

    # -----------------------------------------------------------------

    def set_parallelization_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the optimal parallelization scheme for local execution ...")

        # Uniform
        try: self.set_parallelization_local_uniform()

        # Different
        except DifferentMemoryUsages: self.set_parallelization_local_different()

    # -----------------------------------------------------------------

    def set_parallelization_local_uniform(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Determining one parallelization scheme for all local simulations ...")

        # Determine the number of processes
        processes = self.get_nprocesses_local()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        free_cpus = monitoring.free_cpus()
        threads = int(free_cpus / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:
            log.warning("The number of processes was " + str(processes) + " but the number of free CPU's is only " + str(free_cpus))
            processes = max(int(free_cpus), 1)
            log.warning("Adjusting the number of processes to " + str(processes) + " ...")
            threads = 1

        # Determine number of cores
        cores = processes * threads
        threads_per_core = 2

        # Debugging
        log.debug("The number of cores is " + str(cores))
        log.debug("The number of thread per core is " + str(threads_per_core))
        log.debug("The number of processes is " + str(processes))

        # Data parallelization?
        if self.data_parallel_local is not None: data_parallel = self.data_parallel_local
        else: data_parallel = self.config.data_parallel_local

        # Set the parallelization scheme
        self.parallelization_local = Parallelization(cores, threads_per_core, processes, data_parallel=data_parallel)

        # Debugging
        log.debug("The parallelization scheme for local execution is " + str(self.parallelization_local))

    # -----------------------------------------------------------------

    def set_parallelization_local_different(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Determining the best parallelization scheme for each local simulation separately ...")

        # Data parallelization?
        if self.data_parallel_local is not None: data_parallel = self.data_parallel_local
        else: data_parallel = self.config.data_parallel_local

        # Loop over the simulations in the local queue
        for definition, simulation_name, _ in self.local_queue:

            # NO MPI
            if not introspection.has_mpi():

                # Check nprocesses
                if self.has_nprocesses_local and self.nprocesses_local > 1: raise ValueError("The number of processes that is specified is not possible: MPI installation not present")

                # Set number of processes to 1
                processes = 1

            # MPI present and number of processes is defined
            elif self.has_nprocesses_local: processes = self.nprocesses_local

            # MPI present and number of processes not defined
            else:

                # Estimate the memory
                memory = estimate_memory(definition.ski_path, input_path=definition.input_path, ncells=self.ncells)

                # Determine the number of possible nprocesses
                processes = get_possible_nprocesses_in_memory(monitoring.free_memory(), memory.serial, memory.parallel, data_parallel=data_parallel)

            # Calculate the maximum number of threads per process based on the current cpu load of the system
            free_cpus = monitoring.free_cpus()
            threads = int(free_cpus / processes)

            # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
            if threads < 1:
                log.warning("The number of processes was " + str(processes) + " but the number of free CPU's is only " + str(free_cpus))
                processes = max(int(free_cpus), 1)
                log.warning("Adjusting the number of processes to " + str(processes) + " ...")
                threads = 1

            # Determine number of cores
            cores = processes * threads
            threads_per_core = 2

            # Debugging
            log.debug("The number of cores is " + str(cores))
            log.debug("The number of thread per core is " + str(threads_per_core))
            log.debug("The number of processes is " + str(processes))

            # Set the parallelization scheme
            parallelization_simulation = Parallelization(cores, threads_per_core, processes, data_parallel=data_parallel)

            # Debugging
            log.debug("The parallelization scheme for the '" + simulation_name + "' simulation in the local queue is " + str(parallelization_simulation))

            # Set the parallelization scheme for this simulation
            self.set_parallelization_for_simulation(simulation_name, parallelization_simulation)

    # -----------------------------------------------------------------

    @property
    def has_nprocesses_local(self):

        """
        This function ...
        :return:
        """

        return self.nprocesses_local is not None

    # -----------------------------------------------------------------

    @property
    def local_simulation_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for definition, simulation_name, _ in self.local_queue: names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @lazyproperty
    def memory_for_local_simulations(self):

        """
        Thisn function ...
        :return:
        """

        #simulation_names = self.simulations

        memories = dict()

        # If memory requirement is not set
        if self.memory is None:

            # Estimate the memory usage for each simulation in the local queue
            for definition, simulation_name, _ in self.local_queue:

                # Estimate the memory
                memory = estimate_memory(definition.ski_path, input_path=definition.input_path, ncells=self.ncells)

                # Add to list
                #memories.append(memory)
                memories[simulation_name] = memory

        # Get the memory usage
        else:
            for name in self.local_simulation_names:
                # TODO: what is supposed to be implemented here?
                # Get the memory required estimated on previous local simulations?
                # self.memory is a memory table for multiple simulations
                #memories[name] = self.memory
                pass

        # Return the memory usages
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def serial_memory_for_local_simulations(self):

        """
        This function ...
        :return:
        """

        memories = dict()
        for name in self.memory_for_local_simulations: memories[name] = self.memory_for_local_simulations[name].serial
        return memories

    # -----------------------------------------------------------------

    @lazyproperty
    def parallel_memory_for_local_simulations(self):

        """
        This function ...
        :return:
        """

        memories = dict()
        for name in self.memory_for_local_simulations: memories[name] = self.memory_for_local_simulations[name].parallel
        return memories

    # -----------------------------------------------------------------

    def get_memory_for_local_simulations(self):

        """
        This function ...
        :return:
        """

        # Memory usages between simulations in local queue are very similar
        if sequences.all_close(self.serial_memory_for_local_simulations.values(), rtol=0.1) and sequences.all_close(self.parallel_memory_for_local_simulations.values(), rtol=0.1): return self.memory_for_local_simulations[self.local_simulation_names[0]]

        # Not similar
        else: raise DifferentMemoryUsages("Different memory usage for simulations in local queue")

    # -----------------------------------------------------------------

    def get_nprocesses_local(self):

        """
        This function ...
        :return:
        """

        # Check whether MPI is available
        if not introspection.has_mpi():

            # Check nprocesses
            if self.has_nprocesses_local and self.nprocesses_local > 1: raise ValueError("The number of processes that is specified is not possible: MPI installation not present")

            # Set number of processes to 1
            processes = 1

        # MPI present and number of processes is defined
        elif self.has_nprocesses_local: processes = self.nprocesses_local

        # MPI present and number of processes not defined
        else:

            # Get the memory
            memory = self.get_memory_for_local_simulations()

            # Determine the number of possible nprocesses
            processes = get_possible_nprocesses_in_memory(monitoring.free_memory(), memory.serial, memory.parallel, data_parallel=self.config.data_parallel_local)

        # Return
        return processes

    # -----------------------------------------------------------------

    def get_simulation_names_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        names = []
        for definition, simulation_name, _ in self.queues[host_id]: names.append(simulation_name)
        return names

    # -----------------------------------------------------------------

    @property
    def all_remote_simulation_names(self):

        """
        This function ...
        :return:
        """

        names = []
        for host_id in self.host_ids: names.extend(self.get_simulation_names_for_host(host_id))
        return names

    # -----------------------------------------------------------------

    @property
    def all_simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.local_simulation_names + self.all_remote_simulation_names

    # -----------------------------------------------------------------

    def get_first_simulation_name_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_names_for_host(host_id)[0]

    # -----------------------------------------------------------------

    def get_last_simulation_name_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_names_for_host(host_id)[-1]

    # -----------------------------------------------------------------

    def get_simulation_definition(self, host_id, simulation_name):

        """
        Thisf unction ....
        :param host_id:
        :param simulation_name:
        :return:
        """

        for definition, name, _ in self.queues[host_id]:
            if name == simulation_name: return definition
        raise ValueError("No simulation with the name '" + simulation_name + "' in the queue for remote host '" + host_id + "'")

    # -----------------------------------------------------------------

    def get_simulation_ski_path(self, host_id, simulation_name):

        """
        This function ...
        :param host_id:
        :param simulation_name:
        :return:
        """

        return self.get_simulation_definition(host_id, simulation_name).ski_path

    # -----------------------------------------------------------------

    def get_simulation_skifile(self, host_id, simulation_name):

        """
        This function ...
        :param host_id:
        :param simulation_name:
        :return:
        """

        return self.get_simulation_definition(host_id, simulation_name).ski

    # -----------------------------------------------------------------

    def get_first_simulation_definition_for_host(self, host_id):

        """
        Thisf untion ...
        :param host_id:
        :return:
        """

        return self.get_simulation_definition(host_id, self.get_first_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def get_last_simulation_definition_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_definition(host_id, self.get_last_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def get_first_simulation_ski_path_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_ski_path(host_id, self.get_first_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def get_last_simulation_ski_path_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_ski_path(host_id, self.get_last_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def get_first_simulation_skifile_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_skifile(host_id, self.get_first_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def get_last_simulation_skifile_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.get_simulation_skifile(host_id, self.get_last_simulation_name_for_host(host_id))

    # -----------------------------------------------------------------

    def nwavelengths_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        nwavelengths = dict()

        # If number of wavelengths is not set
        if self.nwavelengths is None:

            # Get the number of wavelengths for each simulation in the queue
            for definition, simulation_name, _ in self.queues[host_id]:

                # Get the number of wavelengths
                nwave = definition.nwavelengths

                # Add to dictionary
                nwavelengths[simulation_name] = nwave

        # Get the number of wavelengths for each simulation
        else:
            for name in self.local_simulation_names: nwavelengths[name] = self.nwavelengths

        # Return the dictionary
        return nwavelengths

    # -----------------------------------------------------------------

    def get_nwavelengths_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Check whether the number of wavelengths is the same for each simulation in this host's queue
        nwavelengths = self.nwavelengths_for_host(host_id)
        simulation_names = nwavelengths.keys()
        if sequences.all_equal(nwavelengths.values()): return nwavelengths[simulation_names[0]]

        # Not the same
        else: raise DifferentNwavelengths("Different number of wavelengths for simulations in queue of remote host '" + host_id + "'")

    # -----------------------------------------------------------------

    def get_first_nwavelengths_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        nwavelengths = self.nwavelengths_for_host(host_id)
        return nwavelengths[self.get_first_simulation_name_for_host(host_id)]

    # -----------------------------------------------------------------

    def dustlib_dimensions_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        dimensions = dict()

        # Get the dustlib dimension for each simulation in the queue
        for definition, simulation_name, _ in self.queues[host_id]:

            # Get the dimension
            dimension = definition.dustlib_dimension

            # Add to dictionary
            dimensions[simulation_name] = dimension

        # Return the dictionary
        return dimensions

    # -----------------------------------------------------------------

    def get_dustlib_dimension_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Check whether dimensions are the same for each simulation in this host's queue
        dimensions = self.dustlib_dimensions_for_host(host_id)
        simulation_names = dimensions.keys()
        if sequences.all_equal(dimensions.values()): return dimensions[simulation_names[0]]

        # Not the same
        else: raise DifferentDustLibDimensions("Different dustlib dimensions for simulations in queue of remote host '" + host_id + "'")

    # -----------------------------------------------------------------

    def get_first_dustlib_dimension_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        dimensions = self.dustlib_dimensions_for_host(host_id)
        return dimensions[self.get_first_simulation_name_for_host(host_id)]

    # -----------------------------------------------------------------

    def ncells_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        ncells = dict()

        # Get the number of dust cells for each simulation in the queue
        for definition, simulation_name, _ in self.queues[host_id]:

            # Get the number of dust cells
            ndustcells = definition.ndust_cells

            # Add to dictionary
            ncells[simulation_name] = ndustcells

        # Return the dictionary
        return ncells

    # -----------------------------------------------------------------

    def get_ncells_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Check whether the number of dust cells are the same for each simulation in this host's queue
        ncells = self.ncells_for_host(host_id)
        simulation_names = ncells.keys()
        if sequences.all_equal(ncells.values()): return ncells[simulation_names[0]]

        # Not the same
        else: raise DifferentNDustCells("Different number of dust cells for simulations in queue of remote host '" + host_id + "'")

    # -----------------------------------------------------------------

    def get_first_ncells_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        ncells = self.ncells_for_host(host_id)
        return ncells[self.get_first_simulation_name_for_host(host_id)]

    # -----------------------------------------------------------------

    def get_nnodes_for_host(self, host_id, return_none=False):

        """
        THis function ...
        :param host_id:
        :param return_none:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler:

            # Get the number of nodes from configuration
            if self.config.nnodes_per_host is not None and host_id in self.config.nnodes_per_host: return self.config.nnodes_per_host[host_id]
            elif self.config.nnodes is not None: return self.config.nnodes
            elif return_none: return None
            else: raise ValueError("Number of nodes must be specified for host with scheduling system")

        # Remote does not use a scheduling system
        else: return 1

    # -----------------------------------------------------------------

    def get_nsockets_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return remote.host.cluster.sockets_per_node

        # Remote does not use a scheduling system
        else:
            if self.config.nsockets is not None: return self.config.nsockets
            elif self.config.all_sockets: return int(remote.sockets_per_node) # should be int already
            else: return int(math.floor(remote.free_sockets))

    # -----------------------------------------------------------------

    def get_ncores_for_host(self, host_id):

        """
        This function returns the number of cores per socket for the specified remote host
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return remote.host.cluster.cores_per_socket

        # Remote does not use a scheduling system
        else: return remote.cores_per_socket

    # -----------------------------------------------------------------

    def get_memory_for_host(self, host_id):

        """
        Thisn function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return remote.host.cluster.memory

        # Remote does not use a scheduling system
        else: return remote.free_memory

    # -----------------------------------------------------------------

    def get_mpi_for_host(self, host_id):

        """
        Thins function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return True

        # Remote does not use a scheduling system
        else: return remote.has_mpi

    # -----------------------------------------------------------------

    def get_hyperthreading_for_host(self, host_id):

        """
        Thisn function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        #remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        #if remote.scheduler: return remote.host.use_hyperthreading

        # Remote does not use a scheduling system
        #else: return remote.host.use_hyperthreading

        # Get the remote instance
        remote = self.get_skirt_remote(host_id)

        # Return
        return remote.use_hyperthreading_skirt

    # -----------------------------------------------------------------

    def get_threads_per_core_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return remote.host.cluster.threads_per_core

        # Remote does not use a scheduling system
        else: return remote.threads_per_core

    # -----------------------------------------------------------------

    def get_properties_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # Set properties
        properties = Map()
        properties.nnodes = self.get_nnodes_for_host(remote, return_none=True)
        properties.nsockets = self.get_nsockets_for_host(remote)
        properties.ncores = self.get_ncores_for_host(remote)
        properties.memory = self.get_memory_for_host(remote)
        properties.mpi = self.get_mpi_for_host(remote)
        properties.hyperthreading = self.get_hyperthreading_for_host(remote)
        properties.threads_per_core = self.get_threads_per_core_for_host(remote)

        # Return the properties
        return properties

    # -----------------------------------------------------------------

    def set_parallelization_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization schemes for the different remote hosts (for hosts and simulation for which"
                 " the parallelization scheme has not been defined yet) ...")

        # Loop over the different remote hosts
        for host_id in self.host_ids:

            # Check if any simulation needs a parallelization
            if self.has_parallelization_for_all_simulations_of_host(host_id): continue

            # Debugging
            log.debug("Setting parallelization scheme(s) for remote host '" + host_id + "' ...")

            # Check whether the parallelization has already been defined by the user for this remote host
            if host_id in self.parallelization_hosts:

                # Check the parallelization?
                if self.config.check_parallelization: self.check_parallelization_for_host(host_id)
                else: continue

            # Number of processes is defined for this host
            if self.has_nprocesses_for_host(host_id):

                # Try to set one parallelization for the host
                try: self.set_parallelization_for_remote_from_nprocesses(host_id)
                except (DifferentNwavelengths, DifferentDustLibDimensions) as e: self.set_parallelization_different_for_remote(host_id)

            # The number of processes per node is defined for this host
            elif self.has_nprocesses_per_node_for_host(host_id):

                # Try to set one parallelization for the host
                try: self.set_parallelization_for_remote_from_nprocesses_per_node(host_id)
                except (DifferentNwavelengths, DifferentDustLibDimensions) as e: self.set_parallelization_different_for_remote(host_id)

            # Determine the parallelization scheme for each remote host
            elif self.config.same_requirements: self.set_parallelization_same_for_remote(host_id)

            # Determine the parallelization scheme with the parallelization tool separately for each simulation
            else: self.set_parallelization_different_for_remote(host_id)

    # -----------------------------------------------------------------

    def check_parallelization_for_host(self, host_id):

        """
        Thisn function ...
        :param host_id:
        :return:
        """

        # Inform the user
        log.info("Checking the parallelization scheme for remote host '" + host_id + "' ...")

        # Get the remote instance
        remote = self.get_skirt_remote(host_id)

        # Get the specified parallelization
        parallelization = self.parallelization_hosts[host_id]

        # If the remote host uses a scheduling system, check whether the parallelization options are possible
        # based on the cluster properties defined in the configuration
        if remote.scheduler:

            # Determine the total number of hardware threads that can be used on the remote cluster
            hardware_threads_per_node = remote.cores_per_node
            if remote.use_hyperthreading_skirt: hardware_threads_per_node *= remote.threads_per_core

            # Raise an error if the number of requested threads per process exceeds the number of hardware threads
            # per node
            if parallelization.threads > hardware_threads_per_node: raise RuntimeError("The number of requested threads per process exceeds the number of allowed threads per node")

            # Determine the number of processes per node (this same calculation is also done in JobScript)
            processes_per_node = remote.cores_per_node // parallelization.threads

            # Determine the amount of requested nodes based on the total number of processes and the number of processes per node
            requested_nodes = math.ceil(parallelization.processes / processes_per_node)

            # Raise an error if the number of requested nodes exceeds the number of nodes of the system
            if requested_nodes > remote.nodes: raise RuntimeError("The required number of computing nodes for"
                                                                   "the requested number of processes and threads "
                                                                   "exceeds the existing number of nodes")

        # No scheduling system
        else:

            # Determine the total number of requested threads
            requested_threads = parallelization.processes * parallelization.threads

            # Determine the total number of hardware threads that can be used on the remote host
            hardware_threads = remote.cores_per_node
            if remote.use_hyperthreading_skirt: hardware_threads *= remote.threads_per_core

            # If the number of requested threads is greater than the allowed number of hardware threads, raise
            # an error
            if requested_threads > hardware_threads: raise RuntimeError("The requested number of processes and threads exceeds the total number of hardware threads")

    # -----------------------------------------------------------------

    def get_data_parallel(self, host_id, nprocesses):

        """
        This function ...
        :param host_id:
        :param nprocesses:
        :return:
        """

        # Set data parallelization flag
        if host_id in self.data_parallel_hosts: data_parallel = self.data_parallel_hosts[host_id]
        elif self.config.data_parallel_remote is not None: data_parallel = self.config.data_parallel
        else:

            # User thinks all simulations have the same requirements
            if self.config.same_requirements: data_parallel = self.get_first_nwavelengths_for_host(host_id) >= 10 * nprocesses and self.get_first_dustlib_dimension_for_host(host_id) == 3  # assured to be safe even if simulations differ

            # Check whether we can find that all simulations have the same number of wavelengths and dustlib dimension
            elif self.get_nwavelengths_for_host(host_id) >= 10 * nprocesses and self.get_dustlib_dimension_for_host(host_id) == 3: data_parallel = True  # can return errors for non-uniformity

            # Not data-parallel
            else: data_parallel = False

        # Return the flag
        return data_parallel

    # -----------------------------------------------------------------

    def set_parallelization_for_remote_from_nprocesses(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the number of processes
        nprocesses = self.get_nprocesses_for_host(host_id)

        # Debugging
        log.debug("Determining one parallelization scheme for all simulations of remote host '" + host_id + "' with " + str(nprocesses) + " processes ...")

        # Get host properties
        prop = self.get_properties_for_host(host_id)

        # Get number of sockets to use and the number of cores per socket
        nnodes = prop.nnodes
        nsockets = prop.nsockets # nsockets per node
        ncores = prop.ncores # ncores per socket

        # Check the number of processes and the number of nodes (the total number of processes must be a multiple of the number of nodes)
        if not numbers.is_multiple_of(nprocesses, nnodes): raise ValueError("The total number of processes (" + str(nprocesses) + ") must be a multiple of the number of nodes (" + str(nnodes) + ")")
        nprocesses_per_node = numbers.as_integer_check_division(nprocesses, nnodes)

        # Check the number of processes per node and the number of sockets
        if self.config.allow_multisocket_processes:
            if not numbers.is_multiple_or_divisor_of(nprocesses_per_node, nsockets): raise ValueError("The number of processes per node (" + str(nprocesses_per_node) + ") must be a multiple or divisor of (or equal to) the number of sockets per node (" + str(nsockets) + ")")
        else:
            if not numbers.is_multiple_of(nprocesses_per_node, nsockets): raise ValueError("The number of processes per node (" + str(nprocesses_per_node) + ") must be a multiple of (or equal to) the number of sockets per node (" + str(nsockets) + ")")

        # Determine cores per node and total number of cores
        cores_per_node = nsockets * ncores
        total_ncores = prop.nnodes * cores_per_node

        # Check number of processes
        if nprocesses > cores_per_node: raise ValueError("The number of processes (" + str(nprocesses) + ") cannot be larger than the number of cores per node (" + str(cores_per_node) + ")")

        # Determine other parameters
        ppn = nsockets * ncores
        nprocesses = nprocesses_per_node * nnodes
        ncores_per_process = ppn / nprocesses_per_node
        threads_per_core = prop.threads_per_core if prop.hyperthreading else 1
        threads_per_process = threads_per_core * ncores_per_process

        # Get data parallelization flag
        data_parallel = self.get_data_parallel(host_id, nprocesses)

        # Create the parallelization object
        para = Parallelization.from_mode("hybrid", total_ncores, threads_per_core, threads_per_process=threads_per_process, data_parallel=data_parallel)

        # Set the parallelization for this remote
        self.set_parallelization_for_host(host_id, para)

    # -----------------------------------------------------------------

    def set_parallelization_for_remote_from_nprocesses_per_node(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the number of processes per node
        nprocesses_per_node = self.get_nprocesses_per_node_for_host(host_id)

        # Debugging
        log.debug("Determining one parallelization scheme for all simulations of remote host '" + host_id + "' with " + str(nprocesses_per_node) + " processes on each node ...")

        # Get host properties
        prop = self.get_properties_for_host(host_id)

        # Get the number of sockets to use and the number of cores per socket
        nsockets = prop.nsockets
        ncores = prop.ncores

        # Check the number of processes per node with the number of sockets per node
        if self.config.allow_multisocket_processes:
            if not numbers.is_multiple_or_divisor_of(nprocesses_per_node, nsockets): raise ValueError("The number of processes per node (" + str(nprocesses_per_node) + ") must be a multiple or divisor of (or equal to) the number of sockets per node (" + str(nsockets) + ")")
        else:
            if not numbers.is_multiple_of(nprocesses_per_node, nsockets): raise ValueError("The number of processes per node (" + str(nprocesses_per_node) + ") must be a multiple of (or equal to) the number of sockets per node (" + str(nsockets) + ")")

        # Determine other parameters
        ppn = nsockets * ncores
        nprocesses = nprocesses_per_node * self.config.nnodes
        ncores_per_process = ppn / nprocesses_per_node
        threads_per_core = prop.threads_per_core if prop.hyperthreading else 1
        threads_per_process = threads_per_core * ncores_per_process
        total_ncores = prop.nnodes * prop.nsockets * prop.ncores

        # Get data parallelization flag
        data_parallel = self.get_data_parallel(host_id, nprocesses)

        # Create the parallelization object
        para = Parallelization.from_mode("hybrid", total_ncores, threads_per_core, threads_per_process=threads_per_process, data_parallel=data_parallel)

        # Set
        self.set_parallelization_for_host(host_id, para)

    # -----------------------------------------------------------------

    def set_parallelization_same_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # TODO: use data parallelization flags?

        # Debugging
        log.debug("Determining one parallelization scheme for all simulations of remote host '" + host_id + "' ...")

        # Determine the number of wavelengths
        if self.nwavelengths is not None: nwavelengths = self.nwavelengths
        else: nwavelengths = self.get_first_nwavelengths_for_host(host_id)

        # Determine the number of dust cells
        if self.ncells is not None: ncells = self.ncells
        else: ncells = self.get_first_ncells_for_host(host_id)

        # Get host properties
        prop = self.get_properties_for_host(host_id)

        # Get first simulation definition for this remote host
        definition = self.get_first_simulation_definition_for_host(host_id)

        # Estimate the memory requirement
        # TODO: use the self.memory table??
        memory = estimate_memory(definition.ski_path, input_path=definition.input_path, ncells=self.ncells)

        # Determine parallelization
        parallelization_remote = determine_parallelization(definition.ski_path, definition.input_path,
                                                    memory, prop.nnodes, prop.nsockets, prop.ncores, prop.memory, prop.mpi,
                                                    prop.hyperthreading, prop.threads_per_core, ncells=ncells,
                                                    nwavelengths=nwavelengths)

        # Show
        log.debug("The parallelization scheme for the '" + host_id + "' host is a " + str(parallelization_remote))

        # Set the parallelization scheme for this host
        self.set_parallelization_for_host(host_id, parallelization_remote)

    # -----------------------------------------------------------------

    def set_parallelization_different_for_remote(self, host_id):

        """
        This function ...
        :param host_id: 
        :return: 
        """

        # TODO: use data parallelization flags?

        # Debugging
        log.debug("Determining the parallelization schemes for host '" + host_id + "' for each simulation separately ...")

        # Loop over the simulations in the queue for the current remote host
        for definition, simulation_name, _ in self.queues[host_id]:

            # Debugging
            log.debug("Setting the parallelization scheme for the '" + simulation_name + "' simulation in the queue of remote host '" + host_id + "' ...")

            # Determine the number of wavelengths
            if self.nwavelengths is not None: nwavelengths = self.nwavelengths
            else: nwavelengths = definition.nwavelengths

            # Determine the number of dust cells
            if self.ncells is not None: ncells = self.ncells
            else: ncells = definition.ndust_cells
            #if ncells is None: estimate_ncells() # TODO

            # Get host properties
            prop = self.get_properties_for_host(host_id)

            # Estimate the memory
            # TODO: use the self.memory table??
            memory = estimate_memory(definition.ski_path, input_path=definition.input_path, ncells=self.ncells)

            # ski_path, input_path, memory, nnodes, nsockets, ncores, host_memory, mpi, hyperthreading, threads_per_core, ncells=None
            parallelization_simulation = determine_parallelization(definition.ski_path, definition.input_path,
                                                                memory, prop.nnodes, prop.nsockets, prop.ncores, prop.memory, prop.mpi,
                                                                prop.hyperthreading, prop.threads_per_core, ncells=ncells,
                                                                nwavelengths=nwavelengths)

            # Show
            log.debug("The parallelization scheme for the '" + simulation_name + "' simulation is " + str(parallelization_simulation))

            # Set the parallelization scheme for this simulation
            self.set_parallelization_for_simulation(simulation_name, parallelization_simulation)

    # -----------------------------------------------------------------

    def set_script_path(self, host_id, path):

        """
        This function ...
        :param host_id:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Setting the local script path for remote host '" + host_id + "' to '" + path + "' ...")

        # Set the path to the script paths dictionary for the remote host
        self.script_paths[host_id] = path

    # -----------------------------------------------------------------

    def enable_screen_output(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Debugging
        log.debug("Saving the screen output remotely will be enabled for remote host '" + host_id + "' ...")

        # Add the ID of the host to the list
        self.save_screen_output.append(host_id)

    # -----------------------------------------------------------------

    def has_scheduling_options(self, host_id, name):

        """
        Thisn function ...
        :param host_id:
        :param name:
        :return:
        """

        return host_id in self.scheduling_options and name in self.scheduling_options[host_id]

    # -----------------------------------------------------------------

    def get_scheduling_options(self, host_id, name):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        return self.scheduling_options[host_id][name]

    # -----------------------------------------------------------------

    def has_walltime(self, host_id, name):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        return self.has_scheduling_options(host_id, name) and self.get_scheduling_options(host_id, name).walltime is not None

    # -----------------------------------------------------------------

    def get_walltime(self, host_id, name):

        """
        This function ...
        :param host_id:
        :param name:
        :return:
        """

        return self.scheduling_options[host_id][name].walltime

    # -----------------------------------------------------------------

    def set_walltime(self, host_id, name, walltime):

        """
        Thisf unction ...
        :param host_id:
        :param name:
        :param walltime:
        :return:
        """

        self.scheduling_options[host_id][name].walltime = walltime

    # -----------------------------------------------------------------

    def all_have_walltime(self, host_id):

        """
        This function ....
        :param host_id:
        :return:
        """

        # Loop over the simulation names for this host
        for simulation_name in self.get_simulation_names_for_host(host_id):
            if not self.has_walltime(host_id, simulation_name): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtimes based on the timings of previous simulations ...")

        # Loop over the hosts which use a scheduling system
        for host_id in self.scheduler_host_ids:

            # Debugging
            log.debug("Estimating the runtimes for host '" + host_id + "' ...")

            # The simulations have the same requirements; estimate the runtime for the host
            if self.config.same_requirements: self.estimate_runtimes_same_for_remote(host_id)

            # The simulations can have different requirements: estimate the runtime for each simulation separately
            else: self.estimate_runtimes_different_for_remote(host_id)

    # -----------------------------------------------------------------

    def estimate_runtimes_same_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Debugging
        log.debug("Estimating one runtime for all simulations on remote host '" + host_id + "' ...")

        # Check if all simulations for this host already have a walltime defined
        if self.all_have_walltime(host_id): return

        # Get cluster name
        cluster_name = self.get_clustername_for_host(host_id)

        # Get the parallelization scheme that we have defined for this remote host
        parallelization_host = self.parallelization_for_host(host_id)

        # Get the first simulation definition
        ski = self.get_first_simulation_skifile_for_host(host_id)

        # Visualisation of the distribution of estimated runtimes
        if self.config.runtimes_plot_path is not None: plot_path = fs.join(self.config.runtimes_plot_path, time.unique_name("explorer_runtimes_" + host_id) + ".pdf")
        else: plot_path = None

        # Estimate the runtime
        runtime = self.runtime_estimator.runtime_for(ski, parallelization_host, host_id, cluster_name, nwavelengths=self.nwavelengths, plot_path=plot_path)

        # Debugging
        log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

        # Loop over the simulations, set the walltime
        for simulation_name in self.get_simulation_names_for_host(host_id):

            # Check if walltime is already defined: shouldn't be?
            if self.has_walltime(host_id, simulation_name): continue

            # Set walltime
            self.set_walltime(host_id, simulation_name, runtime)

    # -----------------------------------------------------------------

    def estimate_runtimes_different_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Debugging
        log.debug("Estimating the runtime for each simulation in the queue of remote host '" + host_id + "' separately ...")

        # Get cluster name
        cluster_name = self.get_clustername_for_host(host_id)

        # Loop over the simulations
        for simulation_name in self.get_simulation_names_for_host(host_id):

            # Get the parallelization scheme
            parallelization_item = self.get_parallelization(host_id, simulation_name)

            # Check if walltime is already defined
            if self.has_walltime(host_id, simulation_name): continue

            # Visualisation of the distribution of estimated runtimes
            if self.config.runtimes_plot_path is not None: plot_path = fs.join(self.config.runtimes_plot_path, time.unique_name("explorer_runtimes_" + host_id) + ".pdf")
            else: plot_path = None

            # Get the ski file
            ski = self.get_simulation_skifile(host_id, simulation_name)

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = self.runtime_estimator.runtime_for(ski, parallelization_item, host_id, cluster_name, nwavelengths=self.nwavelengths, plot_path=plot_path)

            # Debugging
            log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

            # Set the walltime
            self.set_walltime(host_id, simulation_name, runtime)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # The complete list of simulations
        simulations = []

        # Launch locally
        if self.has_queued_local: simulations += self.launch_local()

        # Launch remotely
        if self.has_queued_remote: simulations += self.launch_remote()

        # Set the launched simulations
        self.launched_simulations = simulations

    # -----------------------------------------------------------------

    def launch_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching simulations locally ...")

        # Initialize a list of simulations
        simulations = []

        # Get the number of simulations in the local queue
        nqueued = self.in_local_queue

        # Loop over the simulation in the queue for this remote host
        for index in range(nqueued):

            # Get the last item from the queue (it is removed)
            definition, name, analysis_options_item = self.local_queue.pop()

            # Get the parallelization scheme that has been defined for this simulation
            parallelization_item = self.get_parallelization_for_simulation(name)
            if parallelization_item is not None: pass # OK
            elif self.parallelization_local is not None: parallelization_item = self.parallelization_local
            else: raise RuntimeError("Parallelization has not been defined for local simulation '" + name + "' and no general parallelization scheme has been set for local execution")

            # Get original definition if applicable
            if name in self.original_local_definitions: original_definition = self.original_local_definitions[name]
            else: original_definition = None

            # Generate the analysis options: THIS DOES NOT MODIFY THE DEFINITION
            options_definition = original_definition if original_definition is not None else definition
            logging_options, analysis_options = self.generate_options(name, options_definition, analysis_options_item, local=True)

            # Try launching the simulation
            simulation = self._try_launching_simulation(index, nqueued, name, definition, parallelization_item, logging_options,
                                                        analysis_options, original_definition=original_definition)

            # Add the simulation if succesful
            success = simulation is not None
            if success:

                # Add the simulation to the list
                simulations.append(simulation)

            # Something went wrong: cancel queue?
            elif self.config.cancel_launching_after_fail:

                # Show error message
                log.error("Cancelling following simulations in the local queue ...")

                # Set the assignment for all failed simulations
                self.add_assignment_failed_local()

                # Exit the loop
                break

            # Don't cancel queue
            else:

                # Show warning
                log.warning("Tagging the simulation as unsuccesful but continuing trying to launch other simulations in the local queue ...")

                # Set unsuccesful
                self.assignment.add_local_simulation(name, success=False)

        # Return the list of simulations
        return simulations

    # -----------------------------------------------------------------

    def _try_launching_simulation(self, index, nqueued, name, definition, parallelization, logging_options, analysis_options,
                                  original_definition=None):

        """
        This function ...
        :param index:
        :param name:
        :param definition:
        :param parallelization:
        :param logging_options:
        :param analysis_options:
        :return:
        """

        # Perform the simulation locally
        try: simulation = self._launch_simulation(index, nqueued, name, definition, parallelization, logging_options,
                                                  analysis_options, original_definition=original_definition)

        # Error occured during simulation
        except Exception:

            # Show error messages and traceback; cancel all following simulations
            log.error("Launching simulation '" + name + "' failed:")
            traceback.print_exc()

            # Add the failed simulation back to the queue
            failed_item = (definition, name, analysis_options) # was (definition, name, analysis_options_item)
            self.local_queue.append(failed_item)

            # Set the simulation to None
            simulation = None

        # Return the simulation
        return simulation

    # -----------------------------------------------------------------

    def _launch_simulation(self, index, nqueued, name, definition, parallelization, logging_options, analysis_options,
                           original_definition=None):

        """
        This function ...
        :param index:
        :param name:
        :param definition:
        :param parallelization:
        :param logging_options:
        :param analysis_options:
        :param original_definition:
        :return:
        """

        # Inform the user
        log.info("Launching simulation " + str(index + 1) + " out of " + str(nqueued) + " in the local queue ...")

        # Dry: don't actually run the simulation
        if self.config.dry:

            # Show warning
            log.warning("[DRY] Not actually launching simulation '" + name + "' in dry mode")

            # Create the arguments
            arguments = SkirtArguments.from_definition(definition, logging_options, parallelization)

            # Get the command string
            command = arguments.to_command(self.skirt.scheduler, skirt_path=self.skirt.path, mpirun_path=self.skirt.mpi_command)

            # Show the command
            log.warning("[DRY] To launch the simulation, the command is: '" + command + "'")

            # Create the simulation object
            simulation = SkirtSimulation(definition.prefix, inpath=arguments.input_path, outpath=arguments.output_path, ski_path=definition.ski_path, name=name)

        # Run the simulation
        else: simulation = self.skirt.run(definition, logging_options=logging_options,
                                        parallelization=parallelization, silent=(not log.is_debug),
                                        show_progress=self.config.show_progress)

        # Overwrite the simulation object when the definition had been altered by this class
        if original_definition is not None:

            # Get modified prefix and original prefix
            prefix = simulation.prefix()
            original_prefix = original_definition.prefix

            # Change the names of the output files so that they start with the right prefix (and not the timestamped prefix of the temporarily created ski file)
            for filename in fs.files_in_path(simulation.output_path, returns="name", extensions=True):
                if not filename.startswith(prefix): continue
                original_filename = filename.replace(prefix, original_prefix)
                fs.rename_file(simulation.output_path, filename, original_filename)

            # Create new simulation object
            arguments = SkirtArguments.from_definition(original_definition, logging_options, parallelization)
            simulation = arguments.simulations(simulation_name=name)

        # Success
        log.success("Finished simulation " + str(index + 1) + " out of " + str(nqueued) + " in the local queue ...")

        # Set the parallelization scheme
        simulation.parallelization = parallelization

        # Set the analysis options
        simulation.set_analysis_options(analysis_options)

        # Add analyser classes
        if self.config.analysers is not None:
            for class_path in self.config.analysers: simulation.add_analyser(class_path)

        # Also add the simulation directly to the list of simulations to be analysed
        self.simulations.append(simulation)

        # Add the simulation to the assignment table
        self.assignment.add_local_simulation(name, success=True)

        # Return the simulation
        return simulation

    # -----------------------------------------------------------------

    def add_assignment_failed_local(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting local simulations that failed to be assigned ...")

        # Loop over the leftover simulations
        for definition, name, analysis_options_item in self.local_queue:

            # Add entry to the assignment table
            self.assignment.add_local_simulation(name, success=False)

    # -----------------------------------------------------------------

    def launch_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching simulations on the remote hosts ...")

        # Initialize a list of simulations
        simulations = []

        # Loop over the different remotes
        for remote in self.remotes:

            # Schedule simulations for the current remote host
            simulations_remote = self.schedule_simulations(remote)

            # Start the simulations
            self.start_simulations(remote, simulations_remote)

            # Add the simulations of this remote to the total list of simulations
            simulations += simulations_remote

        # Return the list of simulations
        return simulations

    # -----------------------------------------------------------------

    def get_parallelization(self, host_id, simulation_name):

        """
        This function ...
        :param host_id:
        :param simulation_name:
        :return:
        """

        # Get the parallelization scheme that has been defined for this simulation
        parallelization_item = self.get_parallelization_for_simulation(simulation_name)

        # If no parallelization scheme has been defined for this simulation, use the parallelization scheme
        # defined for the current remote host
        if parallelization_item is None: return self.parallelization_for_host(host_id) # Get the parallelization scheme for this remote host

        # Return individual parallelization
        else: return parallelization_item

    # -----------------------------------------------------------------

    def schedule_simulations(self, remote):

        """
        This function ...
        :param remote
        :return: 
        """

        # Inform the user
        log.info("Scheduling simulations on remote host '" + remote.host_id + "' ...")

        # The remote input path
        if self.has_remote_input_path_for_host(remote.host_id): remote_input_path = self.get_remote_input_path_for_host(remote.host_id)
        else: remote_input_path = None

        # Cache the simulation objects scheduled to the current remote
        simulations_remote = []

        # Get the number of queued simulations for this host
        nqueued = self.in_queue_for_host(remote.host_id)

        # Loop over the simulation in the queue for this remote host
        for index in range(nqueued):

            # Get the last item from the queue (it is removed)
            definition, name, analysis_options_item = self.queues[remote.host_id].pop()

            # Get the parallelization
            parallelization_item = self.get_parallelization(remote.host_id, name)

            # Check whether scheduling options are defined for this simulation and for this remote host
            if self.has_scheduling_options(remote.host_id, name): scheduling_options = self.scheduling_options[remote.host_id][name]
            else: scheduling_options = None

            # Generate the analysis options
            logging_options, analysis_options = self.generate_options(name, definition, analysis_options_item)

            # Queue the simulation
            simulation = self._try_scheduling_simulation(index, nqueued, name, definition, remote, parallelization_item,
                                                         logging_options, analysis_options, scheduling_options,
                                                         remote_input_path=remote_input_path)

            # Add the simulation if succesful
            success = simulation is not None
            if success:

                # Add the simulation to list
                simulations_remote.append(simulation)

                # If the input directory is shared between the different simulations
                if self.config.shared_input and remote_input_path is None: remote_input_path = simulation.remote_input_path

            # Something went wrong: cancel queue?
            elif self.config.cancel_scheduling_after_fail:

                # Show error message
                log.error("Cancelling following simulations in the queue for remote host '" + remote.host_id + "' ...")

                # Set the assignment for all failed simulations
                self.add_assignment_failed_for_host(remote.host_id)

                # Exit the loop
                break

            # Don't cancel queue
            else:

                # Show warning
                log.warning("Tagging the simulation scheduling as unsuccesful but continuing trying to schedule other simulations in the queue of remote host '" + remote.host_id + "' ...")

                # Add entry to assignment table
                cluster_name = self.get_clustername_for_host(remote.host_id)
                self.assignment.add_remote_simulation(name, remote.host_id, cluster_name=cluster_name, success=False)

        # Return the simulations scheduled for the host
        return simulations_remote

    # -----------------------------------------------------------------

    def _try_scheduling_simulation(self, index, nqueued, name, definition, remote, parallelization, logging_options,
                                   analysis_options, scheduling_options=None, remote_input_path=None):

        """
        This function ...
        :param index:
        :param nqueued:
        :param name:
        :param definition:
        :param remote:
        :param parallelization:
        :param logging_options:
        :param analysis_options:
        :param scheduling_options:
        :param remote_input_path:
        :return:
        """

        # Try scheduling
        try: simulation = self._schedule_simulation(index, nqueued, name, definition, remote, parallelization, logging_options,
                                                    analysis_options, scheduling_options, remote_input_path)

        # Exception was raised
        except Exception as e:

            # Print error
            log.error("Adding simulation '" + name + "' to the queue failed:")
            traceback.print_exc()
            log.error(str(e))

            # Add the failed simulation back to the queue
            failed_item = (definition, name, analysis_options) # was (definition, name, analysis_options_item)
            self.queues[remote.host_id].append(failed_item)

            # Set simulation to None
            simulation = None

        # Return the simulation
        return simulation

    # -----------------------------------------------------------------

    def _schedule_simulation(self, index, nqueued, name, definition, remote, parallelization, logging_options,
                             analysis_options, scheduling_options=None, remote_input_path=None):

        """
        This function ...
        :param index:
        :param nqueued:
        :param name:
        :param definition:
        :param remote:
        :param parallelization:
        :param logging_options:
        :param analysis_options:
        :param scheduling_options:
        :param remote_input_path:
        :return:
        """

        # Inform the user
        log.info("Adding simulation " + str(index + 1) + " out of " + str(nqueued) + " to the queue of remote host " + remote.host_id + " ...")

        # Add to the queue
        simulation = remote.add_to_queue(definition, logging_options, parallelization, name=name,
                                         scheduling_options=scheduling_options, remote_input_path=remote_input_path,
                                         analysis_options=analysis_options, emulate=self.config.emulate, dry=self.config.dry,
                                         clear_existing=self.config.clear_existing)

        # Success
        log.success("Added simulation " + str(index + 1) + " out of " + str(nqueued) + " to the queue of remote host " + remote.host_id)

        # Set the parallelization scheme of the simulation (important since SKIRTRemote does not know whether
        # hyperthreading would be enabled if the user provided the parallelization_item when adding the
        # simulation to the queue
        simulation.parallelization = parallelization

        # SET OPTIONS

        # Remove remote files
        simulation.remove_remote_input = not self.config.keep and not self.config.shared_input
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep and not self.config.shared_input

        # Retrieval
        simulation.retrieve_types = self.config.retrieve_types

        # Add analyser classes
        if self.config.analysers is not None:
            for class_path in self.config.analysers: simulation.add_analyser(class_path)

        # Save the simulation object
        if not self.config.dry: simulation.save()
        else: log.warning("[DRY] Not saving the simulation object ...")

        # Get cluster name
        cluster_name = self.get_clustername_for_host(remote.host_id)

        # Add to assignment
        self.assignment.add_remote_simulation(name, remote.host_id, cluster_name=cluster_name, simulation_id=simulation.id, success=None)  # We don't know yet whether launching will actually be succesful!

        # Return the simulation
        return simulation

    # -----------------------------------------------------------------

    def add_assignment_failed_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Debugging
        log.debug("Setting simulations that failed to be assigned to remote host '" + host_id + "' ...")

        # Get cluster name
        cluster_name = self.get_clustername_for_host(host_id)

        # Loop over the leftover simulations
        for definition, name, analysis_options_item in self.queues[host_id]:

            # Add entry to assignment table
            self.assignment.add_remote_simulation(name, host_id, cluster_name=cluster_name, success=False)

    # -----------------------------------------------------------------

    def start_simulations(self, remote, simulations_remote):

        """
        This function ...
        :param remote:
        :param simulations_remote:
        :return: 
        """

        # Inform the user
        log.info("Starting the simulations on remote host '" + remote.host_id + "' ...")

        # Determine queue name (name of the screen session or the remote simulation queue)
        queue_name = time.unique_name("batch_launcher")

        # Set a path for the script file to be saved to locally (for manual inspection)
        if remote.host_id in self.script_paths: local_script_path = fs.join(self.script_paths[remote.host_id], queue_name + "_" + remote.host_id  + ".sh")
        else: local_script_path = None

        # Set a path for the screen output to be saved remotely (for debugging)
        if remote.host_id in self.save_screen_output:
            remote_skirt_dir_path = remote.skirt_dir
            remote_skirt_screen_output_path = fs.join(remote_skirt_dir_path, "screen")
            if not remote.is_directory(remote_skirt_screen_output_path): remote.create_directory(remote_skirt_screen_output_path)
            this_screen_output_path = fs.join(remote_skirt_screen_output_path, queue_name)
            remote.create_directory(this_screen_output_path)
            screen_output_path = this_screen_output_path
        else: screen_output_path = None

        # Try starting the queue
        try: self._start_queue(remote, simulations_remote, queue_name, screen_output_path, local_script_path)

        # An error occured
        except Exception as e:

            # Print error
            log.error("Launching queue of simulations for remote host '" + remote.host_id + "' failed")
            traceback.print_exc()
            log.error(str(e))

            # Set the assignment for all failed simulations
            self.set_assignment_failed_for_host(remote.host_id)

    # -----------------------------------------------------------------

    def _start_queue(self, remote, simulations, queue_name, screen_output_path, local_script_path):

        """
        This function ...
        :param remote:
        :param simulations:
        :param queue_name:
        :param screen_output_path:
        :param local_script_path:
        :return:
        """

        # Set scheduling method
        if remote.scheduler:
            if self.config.group_simulations: scheduling_method = "group"
            else: scheduling_method = "separate"
        else: scheduling_method = None

        # Get jobscripts directory path
        jobscripts_path = self.script_paths[remote.host_id] if remote.host_id in self.script_paths else None

        # Start the queue on the remote host
        handles = remote.start_queue(queue_name=queue_name, schedule_method=scheduling_method,
                                     group_walltime=self.config.group_walltime, local_script_path=local_script_path,
                                     screen_output_path=screen_output_path, jobscripts_path=jobscripts_path,
                                     attached=self.config.attached, dry=self.config.dry)

        # SET THE EXECUTION HANDLES
        # Loop over the simulation for this remote
        for simulation in simulations:

            # If all simulations should have the same handle
            if isinstance(handles, ExecutionHandle): simulation.handle = handles
            else: simulation.handle = handles[simulation.name]  # get the handle for this particular simulation

            # Save the simulation with the execution handle
            if not self.config.dry: simulation.save()
            else: log.warning("[DRY] Not saving the simulation object ...")

    # -----------------------------------------------------------------

    @property
    def nassigned(self):

        """
        This function ...
        :return:
        """

        return len(self.assignment)

    # -----------------------------------------------------------------

    @property
    def assigned_local(self):

        """
        This function ...
        :return:
        """

        return self.assignment.local

    # -----------------------------------------------------------------

    def set_assignment_failed_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Debugging
        log.debug("Indicating that launching the simulations for remote host '" + host_id + "' failed ...")

        # Loop over all the simulations in the assignment table
        for index in range(self.nassigned):

            # Set success flag to False
            self.assignment.set_success(index, False)

    # -----------------------------------------------------------------

    def generate_options(self, name, definition, analysis_options_item, local=False):

        """
        This function ...
        :param name:
        :param definition:
        :param analysis_options_item:
        :param local:
        :return:
        """

        # Set add_timing and add_memory settings
        if local:
            add_timing = self.config.add_timing and self.config.add_timing_local
            add_memory = self.config.add_memory and self.config.add_memory_local
        else:
            add_timing = self.config.add_timing
            add_memory = self.config.add_memory

        # Check whether analysis options are specified
        if analysis_options_item is None:

            # Get a copy of the default logging options, create the default analysis options and adjust logging options according to the analysis options
            if self.has_logging_options_for_simulation(name): logging_options = self.get_logging_options_for_simulation(name)
            else: logging_options = self.logging_options.copy()

            # Create analysis options
            analysis_options_item = self.create_analysis_options(definition, name, logging_options, add_timing=add_timing, add_memory=add_memory)

        # AnalysisOptions object
        elif isinstance(analysis_options_item, AnalysisOptions):

            # Get the default logging options
            if self.has_logging_options_for_simulation(name): logging_options = self.get_logging_options_for_simulation(name)
            else: logging_options = self.logging_options

        # Dict-like analysis options
        elif isinstance(analysis_options_item, dict):

            # Create the default analysis options from the configuration of the batch launcher
            default_analysis_options = self.create_analysis_options(definition, name, add_timing=add_timing, add_memory=add_memory)

            # Set the options specified in the analysis_options_item dictionary
            default_analysis_options.set_options(analysis_options_item)

            # Set the modified default analysis options object as the analysis options for this item
            analysis_options_item = default_analysis_options

            # Get logging options
            if self.has_logging_options_for_simulation(name): logging_options = self.get_logging_options_for_simulation(name)
            else: logging_options = self.logging_options.copy()

            # Check the analysis options
            analysis_options_item.check(logging_options)

        # Invalid type for the analysis options
        else: raise ValueError("Invalid type for the analysis options of simulation " + name + ": " + str(type(analysis_options_item)) + " (must be dictionary-like or AnalyisOptions instance)")

        # Return the logging and analysis options
        return logging_options, analysis_options_item

    # -----------------------------------------------------------------

    def try_retrieving(self):

        """
        This function ...
        :return:
        """

        try: self.retrieve()
        except Exception:
            log.error("Retrieving simulations failed:")
            traceback.print_exc()

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving the output of finished simulations ...")

        # Loop over the different remotes
        for remote in self.remotes:

            # Debugging info
            log.debug("Retrieving finished simulations from remote " + remote.system_name + " ...")

            # Get a list of the simulations that have been succesfully retrieved and add the corresponding
            # simulation objects to the list
            self.simulations += remote.retrieve()

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulations)

    # -----------------------------------------------------------------

    @property
    def has_simulations(self):

        """
        Thisf unction ...
        :return:
        """

        return self.nsimulations > 0

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Showing basic info
        if self.config.show_info: self.show_info()

        # Show remote status
        if self.config.show_status and self.uses_remotes: self.show_status()

        # Show parallelizations
        if self.config.show_parallelizations: self.show_parallelization()

    # -----------------------------------------------------------------

    def show_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing basic info ...")

        print("")
        print(fmt.green + fmt.underlined + "Basic info:" + fmt.reset)
        print("")

        # Show local queue
        if self.has_queued_local: print(" - number of simulations in local queue: " + tostr(self.in_local_queue))

        # Show remote queues
        if self.has_queued_remote:
            print(" - number of simulations in remote queues: " + tostr(self.in_remote_queues))
            for host_id in self.host_ids: print("    * " + fmt.bold + host_id + fmt.reset + ": " + tostr(self.in_remote_queue(host_id)))

        print("")

    # -----------------------------------------------------------------

    def show_status(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the status of the execution platforms ...")

        # Remotes
        if self.uses_remotes: self.show_status_remotes()

        # Local
        else: self.show_status_local()

    # -----------------------------------------------------------------

    def show_status_remotes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the status of the remote hosts ...")

        # Loop over the remote hosts
        for host_id in self.host_ids:

            # Get host properties
            prop = self.get_properties_for_host(host_id)

            print("")
            print(fmt.yellow + fmt.bold + host_id.upper() + fmt.reset)
            print("")

            # Show the properties
            for name in prop:
                if prop[name] is None: continue
                print(" - " + fmt.bold + name + fmt.reset + ": " + tostr(prop[name]))
            print("")

    # -----------------------------------------------------------------

    def show_status_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the status of the local host ...")

        # Get the current memory and cpu load of the system
        memory = monitoring.free_memory()
        free_cpus = monitoring.free_cpus()

        print("")
        print(fmt.yellow + fmt.bold + "LOCAL" + fmt.reset)
        print("")

        print(" - number of free cores: " + tostr(free_cpus))
        print(" - free virtual memory: " + tostr(memory))
        print("")

    # -----------------------------------------------------------------

    @property
    def has_local_parallelization(self):

        """
        This function ...
        :return:
        """

        return self.parallelization_local is not None

    # -----------------------------------------------------------------

    def show_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parallelization schemes ...")

        # Local
        if self.has_local_parallelization: self.show_local_parallelization()

        # Remotes
        if self.has_parallelization_for_any_host: self.show_parallelization_hosts()

        # Simulations
        if self.has_parallelization_for_any_simulation: self.show_parallelization_simulations()

    # -----------------------------------------------------------------

    def show_local_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parallelization scheme for local execution ...")

        print("")
        print(self.parallelization_local)
        print("")

    # -----------------------------------------------------------------

    def show_parallelization_hosts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parallelization schemes for the remote hosts ...")

        print("")
        print(fmt.green + fmt.underlined + "Host parallelization schemes:" + fmt.reset)
        print("")

        # Loop over the hosts
        for host_id in self.host_ids:

            # Check whether has parallelization
            if not self.has_parallelization_for_host(host_id): continue

            # Show host name
            print(fmt.yellow + fmt.bold + host_id.upper() + fmt.reset)
            print("")

            # Show parallelization
            parallelization = self.get_parallelization_for_host(host_id)
            print(parallelization)
            print("")

    # -----------------------------------------------------------------

    def show_parallelization_simulations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parallelization schemes for the simulations ...")

        print("")
        print(fmt.green + fmt.underlined + "Simulation parallelization schemes:" + fmt.reset)
        print("")

        # Loca:
        if self.has_queued_local: self.show_parallelization_simulations_local()

        # Remote
        if self.has_queued_remotes: self.show_parallelization_simulations_remote()

    # -----------------------------------------------------------------

    def show_parallelization_simulations_local(self):

        """
        This function ...
        :return:
        """

        # Show 'LOCAL'
        print(fmt.yellow + fmt.bold + "LOCAL" + fmt.reset)
        print("")

        # Loop over the local simulations
        for simulation_name in self.local_simulation_names:

            # Check whether has parallelization
            if not self.has_parallelization_for_simulation(simulation_name): continue

            # Show simulation name
            print(fmt.blue + fmt.underlined + simulation_name + fmt.reset)

            # Show parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)
            print(parallelization)
            print("")

    # -----------------------------------------------------------------

    def show_parallelization_simulations_remote(self):

        """
        This function ...
        :return:
        """

        # Loop over the hosts
        for host_id in self.host_ids:
            if not self.has_queued_remote(host_id): continue

            # Show
            self.show_parallelization_simulations_for_remote(host_id)

    # -----------------------------------------------------------------

    def show_parallelization_simulations_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Show host name
        print(fmt.yellow + fmt.bold + host_id.upper() + fmt.reset)
        print("")

        # Loop over the simulations in the queue
        for simulation_name in self.get_simulation_names_for_host(host_id):

            # Check whether has parallelization
            if not self.has_parallelization_for_simulation(simulation_name): continue

            # Show simulation name
            print(fmt.blue + fmt.underlined + simulation_name + fmt.reset)

            # Show parallelization
            parallelization = self.get_parallelization_for_simulation(simulation_name)
            print(parallelization)
            print("")

    # -----------------------------------------------------------------

    def show_finished(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the output of retrieved simulations ...")

        # Loop over the simulations
        print("")
        for simulation in self.simulations:

            # Print the simulation name
            print(fmt.blue + simulation.prefix() + fmt.reset + ":")

            # Show the output
            simulation.output.show(line_prefix="  ")

    # -----------------------------------------------------------------

    def try_analysing(self):

        """
        This function ...
        :return:
        """

        try: self.analyse()
        except Exception, err:
            log.error("Analysing finished simulations failed:")
            traceback.print_exc()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the output of retrieved simulations ...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation=simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the simulation assignment
        if self.config.write_assignment: self.write_assignment()

        # Write queues
        if self.config.write_queues: self.write_queues()

    # -----------------------------------------------------------------

    def write_assignment(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulation assignment ...")

        # Determine the path
        path = self.output_path_file("assignment.dat")

        # Write
        self.assignment.saveto(path)

    # -----------------------------------------------------------------

    def write_queues(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the queues ...")

        # Write local queue, if any
        if self.has_queued_local: self.write_local_queue()

        # Write remote queues
        if self.has_queued_remotes: self.write_remote_queues()

    # -----------------------------------------------------------------

    def write_local_queue(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing local queue ...")

        # Determine path
        path = self.output_path_file("queue_local.dat")

        # Write
        #write_list(self.local_queue, path)

        # Create table
        table = QueuedSimulationsTable.from_queue(self.local_queue)

        # Save
        table.saveto(path)

    # -----------------------------------------------------------------

    def write_remote_queues(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing remote queues ...")

        # Loop over the remotes
        for host_id in self.queues:

            # Any still queued?
            if not self.has_queued_remote(host_id): continue

            # Determine path
            path = self.output_path_file("queue_" + host_id + ".dat")

            # Write queue
            #write_list(self.queues[host_id], path)

            # Create table
            table = QueuedSimulationsTable.from_queue(self.queues[host_id])

            # Save
            table.saveto(path)

    # -----------------------------------------------------------------

    def create_analysis_options(self, definition, simulation_name, logging_options=None, add_timing=True, add_memory=True):

        """
        This function ...
        :param definition:
        :param simulation_name:
        :param logging_options:
        :param add_timing:
        :param add_memory:
        :return:
        """

        # Create the analysis options object
        analysis_options = AnalysisOptions()
        analysis_options.set_options(self.config.analysis)

        # Determine the path to the extraction directory for this simulation
        if self.config.analysis.extraction.path is None:
            extraction_path = fs.join(definition.base_path, default_extraction_dirname)
            # the ski files of multiple simulations can be int he same directory
            if fs.is_directory(extraction_path): extraction_path = fs.join(definition.base_path, definition.name, default_extraction_dirname)
        elif self.config.relative_analysis_paths: extraction_path = fs.join(definition.base_path, self.config.analysis.extraction.path)
        else: extraction_path = fs.join(self.config.analysis.extraction.path, simulation_name)

        # Create the extraction directory only if it is necessary
        if analysis_options.any_extraction:
            if not fs.is_directory(extraction_path): fs.create_directory(extraction_path, recursive=True)
            analysis_options.extraction.path = extraction_path

        # Set the extraction path to None otherwise
        else: analysis_options.extraction.path = None

        # Determine the path to the plotting directory for this simulation
        if self.config.analysis.plotting.path is None:
            plotting_path = fs.join(definition.base_path, default_plotting_dirname)
            # the ski files of multiple simulations can be in the same directory
            if fs.is_directory(plotting_path): plotting_path = fs.join(definition.base_path, simulation_name, default_plotting_dirname)
        elif self.config.relative_analysis_paths: plotting_path = fs.join(definition.base_path, self.config.analysis.plotting.path)
        else: plotting_path = fs.join(self.config.analysis.plotting.path, simulation_name)

        # Create the plotting directory only if it is necessary
        if analysis_options.any_plotting:
            if not fs.is_directory(plotting_path): fs.create_directory(plotting_path, recursive=True)
            analysis_options.plotting.path = plotting_path

        # Set the plotting path to None otherwise
        else: analysis_options.plotting.path = None

        # Determine the 'misc' directory for this simulation (and create it if necessary)
        if self.config.analysis.misc.path is None:
            misc_path = fs.join(definition.base_path, default_misc_dirname)
            # the ski files of multiple simulations can be int he same directory
            if fs.is_directory(misc_path): misc_path = fs.join(definition.base_path, simulation_name, default_misc_dirname)
        elif self.config.relative_analysis_paths: misc_path = fs.join(definition.base_path, self.config.analysis.misc.path)
        else: misc_path = fs.join(self.config.analysis.misc.path, simulation_name)

        # Create the misc directory only if it is necessary
        if analysis_options.any_misc:
            if not fs.is_directory(misc_path): fs.create_directory(misc_path, recursive=True)
            analysis_options.misc.path = misc_path

        # Set the misc path to None otherwise
        else: analysis_options.misc.path = None

        # Set timing and memory table paths (if specified for this batch launcher)
        if self.has_timing_table and add_timing:

            # Check path
            if self.timing_table.path is None: raise ValueError("Path of the timing table is undefined")

            # Set the table path
            analysis_options.timing_table_path = self.timing_table.path

            # Check whether the extract timing option has been enabled
            if not analysis_options.extraction.timeline:
                log.warning("Timeline extraction will be enabled for writing to the timing table ...")
                analysis_options.extraction.timeline = True

        # Set memory table path
        if self.has_memory_table and add_memory:

            # Check path
            if self.memory_table.path is None: raise ValueError("Path of the memory table is undefined")

            # Set
            analysis_options.memory_table_path = self.memory_table.path

            # Check whether the extract memory option has been enabled
            if not analysis_options.extraction.memory:
                log.warning("Memory extraction will be enabled for writing to the mrmory table ...")
                analysis_options.extraction.memory = True

        # Check the analysis options
        if logging_options is not None: analysis_options.check(logging_options)

        # Return the analysis options
        return analysis_options

# -----------------------------------------------------------------

def all_input_files_are_shared(simulation_input, shared_filepaths):

    """
    This function ...
    :param simulation_input:
    :param shared_filepaths:
    :return:
    """

    for filename, filepath in simulation_input:
        if filepath not in shared_filepaths: return False
    return True

# -----------------------------------------------------------------
