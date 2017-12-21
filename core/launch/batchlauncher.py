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
from ..remote.host import Host, load_host
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

        # THe parallelization scheme for the different simulations
        self.parallelization_simulations = dict()

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

        return len(self.queues[host_id])

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

    def add_to_queue(self, definition, name, host_id=None, parallelization=None, analysis_options=None, local=False):

        """
        This function ...
        :param definition:
        :param name: a name that is given to the simulation
        :param host_id: the host on which this simulation should be executed
        :param parallelization: individual parallelization scheme for this particular simulation (only allowed if host_id is also specified)
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

        # Local execution
        if local or self.has_no_remotes:

            # Add to the local queue
            self.local_queue.append((definition, name, analysis_options))

            # If parallelization is specified, set it
            if parallelization is not None: self.set_parallelization_for_simulation(name, parallelization)

        # If a host ID is specified
        elif host_id is not None:

            # Check if this is a valid host ID
            if not host_id in self.host_ids: raise ValueError("Invalid host ID")

            # Add to the queue of the specified host
            self.queues[host_id].append((definition, name, analysis_options))

            # If parallelization is specified, set it
            if parallelization is not None: self.set_parallelization_for_simulation(name, parallelization)

        else:

            # Check whether, if parallelization is specified, that host ID is also specified
            if parallelization is not None: raise ValueError("If parallelization is specified, host ID must also be specified")

            # Determine the ID of the hosts with the shortest queue (or the first if the queues are equally long)
            host_id = self.shortest_queue_host_id

            # Add to the queue of the host
            self.queues[host_id].append((definition, name, analysis_options))

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
    def hosts(self):

        """
        This function returns the Host objects for all the hosts that will be used by the batch launcher
        :return:
        """

        # If the setup has not been called yet
        if len(self.remotes) == 0:

            # Create the list of hosts
            hosts = []
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

    def set_parallelization_for_simulation(self, name, parallelization):

        """
        This function ...
        :param name:
        :param parallelization:
        :return:
        """

        # Set the parallelization for this host and this simulation
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

    def has_paralleliation_for_host(self, host_id):

        """
        This funtion ...
        :param host_id:
        :return:
        """

        return host_id in self.parallelization_hosts[host_id] and self.parallelization_hosts[host_id] is not None

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

    def parallelization_for_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.parallelization_simulations[name] if name in self.parallelization_simulations else None

    # -----------------------------------------------------------------

    @lazyproperty
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return TimingTable.from_file(self.config.timing_table_path) if self.config.timing_table_path is not None else None

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

        return self.config.timing_table_path is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return MemoryTable.from_file(self.config.memory_table_path) if self.config.memory_table_path is not None else None

    # -----------------------------------------------------------------

    @property
    def has_memory_table(self):

        """
        This function ...
        :return:
        """

        return self.config.memory_table_path is not None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Check the input files for all simulations
        self.check_input()

        # 3. Set the parallelization scheme for the remote hosts for which this was not specified by the user
        self.set_parallelization()

        # 4. Estimate the runtimes if necessary
        if self.uses_schedulers and self.has_timing_table: self.estimate_runtimes()

        # 5. Launch the simulations
        self.launch()

        # 6. Retrieve the simulations that are finished
        self.try_retrieving()

        # 7. Show the simulations that are finished
        if self.has_simulations and self.config.show: self.show()

        # 8. Analyse the output of the retrieved simulations
        if self.has_simulations: self.try_analysing()

        # 9. Write
        self.write()

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
        for definition in create_definitions(self.config.path, self.config.output, self.config.input, recursive=self.config.recursive):

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

        if self.data_parallel_local is not None: data_parallel = self.data_parallel_local
        else: data_parallel = self.config.data_parallel_local

        # Set the parallelization scheme
        self.parallelization_local = Parallelization(cores, threads_per_core, processes, data_parallel=data_parallel)

        # Debugging
        log.debug("The parallelization scheme for local execution is " + str(self.parallelization_local))

    # -----------------------------------------------------------------

    def set_parallelization_local_different(self):

        """
        Thisj function ...
        :return:
        """

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
            for name in self.local_simulation_names: memories[name] = self.memory

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

    def get_nnodes_for_host(self, host_id):

        """
        THis function ...
        :param host_id:
        :return:
        """

        # Get the remote instance
        remote = self.get_remote(host_id)

        # If the remote uses a scheduling system
        if remote.scheduler: return self.config.nnodes

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
        else: return int(math.floor(remote.free_sockets))

    # -----------------------------------------------------------------

    def get_ncores_for_host(self, host_id):

        """
        This function ...
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
        properties.nnodes = self.get_nnodes_for_host(remote)
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
        for remote in self.remotes:

            # Check whether the parallelization has already been defined by the user for this remote host
            if remote.host_id in self.parallelization_hosts:
                # CHECK THE PARALLELIZATION?
                if self.config.check_parallelization: self.check_parallelization_for_host(remote.host_id)
                else: continue

            # Number of processes is defined for this host
            if self.has_nprocesses_for_host(remote.host_id):

                # Try to set one parallelization for the host
                try: self.set_parallelization_for_remote_from_nprocesses(remote.host_id)
                except (DifferentNwavelengths, DifferentDustLibDimensions) as e: self.set_parallelization_different_for_remote(remote.host_id)

            # The number of processes per node is defined for this host
            elif self.has_nprocesses_per_node_for_host(remote.host_id):

                # Try to set one parallelization for the host
                try: self.set_parallelization_for_remote_from_nprocesses_per_node(remote.host_id)
                except (DifferentNwavelengths, DifferentDustLibDimensions) as e: self.set_parallelization_different_for_remote(remote.host_id)

            # Determine the parallelization scheme for each remote host
            elif self.config.same_requirements: self.set_parallelization_same_for_remote(remote.host_id)

            # Determine the parallelization scheme with the parallelization tool separately for each simulation
            else: self.set_parallelization_different_for_remote(remote.host_id)

    # -----------------------------------------------------------------

    def check_parallelization_for_host(self, host_id):

        """
        Thisn function ...
        :param host_id:
        :return:
        """

        # Inform the user
        log.info("Checking the parallelization scheme ...")

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
            # self.remote.cores = cores per node
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

    def set_parallelization_for_remote_from_nprocesses(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the number of processes
        nprocesses = self.get_nprocesses_for_host(host_id)

        # Get host properties
        prop = self.get_properties_for_host(host_id)

        # Determine cores per node and total number of cores
        cores_per_node = prop.nsockets * prop.ncores
        total_ncores = prop.nnodes * cores_per_node

        # Check number of processes
        if nprocesses > cores_per_node: raise ValueError("The number of processes cannot be larger than the number of cores per node (" + str(cores_per_node) + ")")

        # Determine other parameters
        ppn = prop.nsockets * prop.ncores
        nprocesses_per_node = int(nprocesses / prop.nnodes)
        nprocesses = nprocesses_per_node * prop.nnodes
        ncores_per_process = ppn / nprocesses_per_node
        threads_per_core = prop.threads_per_core if prop.hyperthreading else 1
        threads_per_process = threads_per_core * ncores_per_process

        # Set data parallelization flag
        if host_id in self.data_parallel_hosts: data_parallel = self.data_parallel_hosts[host_id]
        elif self.config.data_parallel_remote is not None: data_parallel = self.config.data_parallel
        else:

            # Determine data-parallel flag

            # User thinks all simulations have the same requirements
            if self.config.same_requirements: data_parallel = self.get_first_nwavelengths_for_host(host_id) >= 10 * nprocesses and self.get_first_dustlib_dimension_for_host(host_id) == 3 # assured to be safe even if simulations differ

            # Check whether we can find that all simulations have the same number of wavelengths and dustlib dimension
            elif self.get_nwavelengths_for_host(host_id) >= 10 * nprocesses and self.get_dustlib_dimension_for_host(host_id) == 3: data_parallel = True # can return errors for non-uniformity

            # Not data-parallel
            else: data_parallel = False

        # Create the parallelization object
        parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
                                                    threads_per_process=threads_per_process,
                                                    data_parallel=data_parallel)

        # Set
        self.set_parallelization_for_host(host_id, parallelization)

    # -----------------------------------------------------------------

    def set_parallelization_for_remote_from_nprocesses_per_node(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # Get the number of processes per node
        nprocesses_per_node = self.get_nprocesses_per_node_for_host(host_id)

        # Get host properties
        prop = self.get_properties_for_host(host_id)

        # Determine other parameters
        ppn = prop.nsockets * prop.ncores
        nprocesses = nprocesses_per_node * self.config.nnodes
        ncores_per_process = ppn / nprocesses_per_node
        threads_per_core = prop.threads_per_core if prop.hyperthreading else 1
        threads_per_process = threads_per_core * ncores_per_process
        total_ncores = prop.nnodes * prop.nsockets * prop.ncores

        # Set data parallelization flag
        if host_id in self.data_parallel_hosts: data_parallel = self.data_parallel_hosts[host_id]
        elif self.config.data_parallel_remote is not None: data_parallel = self.config.data_parallel
        else:

            # Determine data-parallel flag
            #if self.config.data_parallel_remote is None:

            if self.get_nwavelengths_for_host(host_id) >= 10 * nprocesses and self.get_dustlib_dimension_for_host(host_id) == 3: data_parallel = True
            else: data_parallel = False

        # Up to the user
        #else: data_parallel = self.config.data_parallel_remote

        # Create the parallelization object
        parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
                                                         threads_per_process=threads_per_process,
                                                         data_parallel=data_parallel)

        # Set
        self.set_parallelization_for_host(host_id, parallelization)

    # -----------------------------------------------------------------

    def set_parallelization_same_for_remote(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        # TODO: use data parallelization flags?

        # Debugging
        log.debug("Setting the parallelization schemes for host '" + host_id + "' ...")

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

        # Determine parallelization
        parallelization = determine_parallelization(definition.ski_path, definition.input_path,
                                                    self.memory, prop.nnodes, prop.nsockets, prop.ncores, prop.memory, prop.mpi,
                                                    prop.hyperthreading, prop.threads_per_core, ncells=ncells,
                                                    nwavelengths=nwavelengths)

        # Show
        log.debug("The parallelization scheme for the '" + host_id + "' host is a " + str(parallelization))

        # Set the parallelization scheme for this host
        self.set_parallelization_for_host(host_id, parallelization)

    # -----------------------------------------------------------------

    def set_parallelization_different_for_remote(self, host_id):

        """
        This function ...
        :param host_id: 
        :return: 
        """

        # TODO: use data parallelization flags?

        # Debugging
        log.debug("Setting the parallelization schemes for host '" + host_id + "' ...")

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

            # ski_path, input_path, memory, nnodes, nsockets, ncores, host_memory, mpi, hyperthreading, threads_per_core, ncells=None
            parallelization_simulation = determine_parallelization(definition.ski_path, definition.input_path,
                                                                self.memory, prop.nnodes, prop.nsockets, prop.ncores, prop.memory, prop.mpi,
                                                                prop.hyperthreading, prop.threads_per_core, ncells=ncells,
                                                                nwavelengths=nwavelengths)

            # Show
            log.debug("The parallelization scheme for the '" + simulation_name + "' simulation is " + str(parallelization))

            # Set the parallelization scheme for this simulation
            self.set_parallelization_for_simulation(simulation_name, parallelization)

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

        return host_id in self.scheduling_options and name in self.scheduling_options[name]

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

        # Loop over the simulation in the queue for this remote host
        total_queued = len(self.local_queue)
        for index in range(total_queued):

            # Get the last item from the queue (it is removed)
            definition, name, analysis_options_item = self.local_queue.pop()

            # Get the parallelization scheme that has been defined for this simulation
            parallelization_item = self.parallelization_for_simulation(name)
            if parallelization_item is not None: pass # OK
            elif self.parallelization_local is not None: parallelization_item = self.parallelization_local
            else: raise RuntimeError("Parallelization has not been defined for local simulation '" + name + "' and no general parallelization scheme has been set for local execution")

            # Get original definition if applicable
            if name in self.original_local_definitions: original_definition = self.original_local_definitions[name]
            else: original_definition = None

            # Generate the analysis options: THIS DOES NOT MODIFY THE DEFINITION
            options_definition = original_definition if original_definition is not None else definition
            logging_options, analysis_options = self.generate_options(name, options_definition, analysis_options_item, local=True)

            # Perform the simulation locally
            try:

                # Inform the user
                log.info("Launching simulation " + str(index + 1) + " out of " + str(total_queued) + " in the local queue ...")

                # Run the simulation
                simulation = self.skirt.run(definition, logging_options=logging_options,
                                            parallelization=parallelization_item, silent=(not log.is_debug()),
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
                    arguments = SkirtArguments.from_definition(original_definition, logging_options, parallelization_item)
                    simulation = arguments.simulations(simulation_name=name)

                # Success
                log.success("Finished simulation " + str(index + 1) + " out of " + str(total_queued) + " in the local queue ...")

                # Set the parallelization scheme
                simulation.parallelization = parallelization_item

                # Set the analysis options
                simulation.set_analysis_options(analysis_options)

                # Add analyser classes
                if self.config.analysers is not None:
                    for class_path in self.config.analysers: simulation.add_analyser(class_path)

                # Add the simulation to the list
                simulations.append(simulation)

                # Also add the simulation directly to the list of simulations to be analysed
                self.simulations.append(simulation)

            # Error occured during simulation
            except Exception:

                # Show error messages and traceback; cancel all following simulations
                log.error("Launching simulation '" + name + "' failed:")
                traceback.print_exc()
                log.error("Cancelling following simulations in the queue ...")
                break

        # Return the list of simulations
        return simulations

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
        parallelization_item = self.parallelization_for_simulation(simulation_name)

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
        remote_input_path = None

        # Cache the simulation objects scheduled to the current remote
        simulations_remote = []

        # Loop over the simulation in the queue for this remote host
        total_queued_host = len(self.queues[remote.host_id])
        for index in range(total_queued_host):

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
            try:

                # Inform the user
                log.info("Adding simulation " + str(index + 1) + " out of " + str(total_queued_host) + " to the queue of remote host " + remote.host_id + " ...")

                # Add to the queue
                simulation = remote.add_to_queue(definition, logging_options, parallelization_item, name=name,
                                                 scheduling_options=scheduling_options, remote_input_path=remote_input_path,
                                                 analysis_options=analysis_options, emulate=self.config.emulate)
                simulations_remote.append(simulation)

                # Success
                log.success("Added simulation " + str(index + 1) + " out of " + str(total_queued_host) + " to the queue of remote host " + remote.host_id)

                # Set the parallelization scheme of the simulation (important since SKIRTRemote does not know whether
                # hyperthreading would be enabled if the user provided the parallelization_item when adding the
                # simulation to the queue
                simulation.parallelization = parallelization_item

                # If the input directory is shared between the different simulations
                if self.config.shared_input and remote_input_path is None: remote_input_path = simulation.remote_input_path

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
                simulation.save()

            # Exception was raised
            except Exception:

                log.error("Adding simulation '" + name + "' to the queue failed:")
                traceback.print_exc()
                log.error("Cancelling following simulations in the queue for remote host '" + remote.host_id + "' ...")
                break

        # Return the simulations scheduled for the host
        return simulations_remote

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

        # Start the queue
        jobscripts_path = self.script_paths[remote.host_id] if remote.host_id in self.script_paths else None
        handles = remote.start_queue(queue_name=queue_name, group_simulations=self.config.group_simulations,
                                     group_walltime=self.config.group_walltime, use_pts=self.config.use_pts,
                                     local_script_path=local_script_path, screen_output_path=screen_output_path,
                                     jobscripts_path=jobscripts_path, attached=self.config.attached, dry=self.config.dry)

        # SET THE EXECUTION HANDLES
        # Loop over the simulation for this remote
        for simulation in simulations_remote:

            # If all simulations should have the same handle
            if isinstance(handles, ExecutionHandle): simulation.handle = handles
            else: simulation.handle = handles[simulation.name] # get the handle for this particular simulation
            simulation.save()

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
            logging_options = self.logging_options.copy()

            # Create analysis options
            analysis_options_item = self.create_analysis_options(definition, name, logging_options, add_timing=add_timing, add_memory=add_memory)

        # AnalysisOptions object
        elif isinstance(analysis_options_item, AnalysisOptions):

            # Get the default logging options
            logging_options = self.logging_options

        # Dict-like analysis options
        elif isinstance(analysis_options_item, dict):

            # Create the default analysis options from the configuration of the batch launcher
            default_analysis_options = self.create_analysis_options(definition, name, add_timing=add_timing, add_memory=add_memory)

            # Set the options specified in the analysis_options_item dictionary
            default_analysis_options.set_options(analysis_options_item)

            # Set the modified default analysis options object as the analysis options for this item
            analysis_options_item = default_analysis_options

            # Get logging options
            logging_options = self.logging_options.copy()

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
        except Exception, err:
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
            extraction_path = fs.join(definition.base_path, "extr")
            # the ski files of multiple simulations can be int he same directory
            if fs.is_directory(extraction_path): extraction_path = fs.join(definition.base_path, definition.name, "extr")
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
            plotting_path = fs.join(definition.base_path, "plot")
            # the ski files of multiple simulations can be in the same directory
            if fs.is_directory(plotting_path): plotting_path = fs.join(definition.base_path, simulation_name, "plot")
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
            misc_path = fs.join(definition.base_path, "misc")
            # the ski files of multiple simulations can be int he same directory
            if fs.is_directory(misc_path): misc_path = fs.join(definition.base_path, simulation_name, "misc")
        elif self.config.relative_analysis_paths: misc_path = fs.join(definition.base_path, self.config.analysis.misc.path)
        else: misc_path = fs.join(self.config.analysis.misc.path, simulation_name)

        # Create the misc directory only if it is necessary
        if analysis_options.any_misc:
            if not fs.is_directory(misc_path): fs.create_directory(misc_path, recursive=True)
            analysis_options.misc.path = misc_path

        # Set the misc path to None otherwise
        else: analysis_options.misc.path = None

        # Set timing and memory table paths (if specified for this batch launcher)
        if self.config.timing_table_path is not None and add_timing:

            # Set the table path
            analysis_options.timing_table_path = self.config.timing_table_path

            # Check whether the extract timing option has been enabled
            if not analysis_options.extraction.timeline:
                log.warning("Timeline extraction will be enabled for writing to the timing table ...")
                analysis_options.extraction.timeline = True

        # Set memory table path
        if self.config.memory_table_path is not None and add_memory:

            # Set
            analysis_options.memory_table_path = self.config.memory_table_path

            # Check whether the extract memory option has been enabled
            if not analysis_options.extraction.memory:
                log.warning("Memory extraction will be enabled for writing to the mrmory table ...")
                analysis_options.extraction.memory = True

        # Check the analysis options
        if logging_options is not None: analysis_options.check(logging_options)

        # Return the analysis options
        return analysis_options

# -----------------------------------------------------------------
