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
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.configurable import OldConfigurable
from ..simulation.remote import SkirtRemote
from .options import LoggingOptions
from ..tools import introspection, time
from ..tools import filesystem as fs
from ..tools.logging import log
from ..basics.host import Host
from .parallelization import Parallelization
from .analyser import SimulationAnalyser
from .options import AnalysisOptions

# -----------------------------------------------------------------

class BatchLauncher(OldConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(BatchLauncher, self).__init__(config, "core")

        # -- Attributes --

        # Initialize a list to contain different SkirtRemote instances for the different remote hosts
        self.remotes = []

        # For some remote hosts, if they use a scheduling system, defines the name of the cluster to be used
        # (the default cluster will be used if not defined)
        self.cluster_names = dict()

        # The queue
        self.queue = []

        # The scheduling options for (some of) the simulations and (some of) the remote hosts. This is a nested
        # dictionary where the first key represents the remote host ID and the next key represents the name of the
        # simulation (see 'name' parameter of 'add_to_queue')
        self.scheduling_options = defaultdict(dict)

        # The assignment from items in the queue to the different remote hosts
        self.assignment = None

        # The parallelization scheme for the different remote hosts
        self.parallelization = dict()

        # The paths to the directories for placing the (screen/job) scripts (for manual inspection) for the different remote hosts
        self.script_paths = dict()

        # The list of remote hosts for which the screen output should be saved remotely (for debugging)
        self.save_screen_output = []

        # The simulations that have been retrieved
        self.simulations = []

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # The logging options
        self.logging_options = None

    # -----------------------------------------------------------------

    def add_to_queue(self, definition, name=None, parallelization=None, analysis_options=None):

        """
        This function ...
        :param definition:
        :param name: a name that is given to the simulation
        :param parallelization: individual parallelization scheme for this particular simulation
        :param analysis_options: analysis options (if None, analysis options will be created from batch launcher configuration)
        :return:
        """

        # Add the simulation definition object to the queue
        self.queue.append((definition, name, parallelization, analysis_options))

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

        return len(self.queue)

    # -----------------------------------------------------------------

    @property
    def single_remote(self):

        """
        This function ...
        :return:
        """

        if len(self.remotes) == 0: raise RuntimeError("The setup has not been called yet")
        elif len(self.remotes) == 1: return self.remotes[0]
        else: raise RuntimeError("Multiple remotes have been configured for this batch launcher")

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
            for host_id in self.host_ids: hosts.append(Host(host_id))

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
                host = Host(host_id)

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
            for id in self.host_ids:

                # Create a Host instance
                host = Host(id)

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
                host = Host(host_id)

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
            for id in self.host_ids:

                # Create a Host instance
                host = Host(id)

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

    def set_parallelization_for_host(self, host_id, parallelization):

        """
        This function ...
        :param host_id:
        :param parallelization:
        :return:
        """

        # Set the parallelization properties for the specified host
        self.parallelization[host_id] = parallelization

    # -----------------------------------------------------------------

    def parallelization_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return self.parallelization[host_id] if host_id in self.parallelization else None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Determine how many simulations are assigned to each remote
        self.assign()

        # 3. Set the parallelization scheme for the remote hosts for which this was not specified by the user
        self.set_parallelization()

        # 4. Launch the simulations
        simulations = self.simulate()

        # 5. Retrieve the simulations that are finished
        self.retrieve()

        # 6. Analyse the output of the retrieved simulations
        self.analyse()

        # 7. Return the simulations that are just scheduled
        return simulations

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

        # Clear the queue
        self.queue = []

        # Clear the scheduling options
        self.scheduling_options = defaultdict(dict())

        # Clear the assignment
        self.assignment = None

        # Clear the simulations
        self.simulations = []

        # Clear the script path dictionary
        self.script_paths = dict()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BatchLauncher, self).setup()

        # Setup the remote instances
        if len(self.remotes) == 0: self.setup_remotes()

        # Create the logging options
        self.logging_options = LoggingOptions()
        self.logging_options.set_options(self.config.logging)

    # -----------------------------------------------------------------

    def setup_remotes(self):

        """
        This function ...
        :return:
        """

        # If a list of remotes is defined
        if self.config.remotes is not None: host_ids = self.config.remotes

        # If a list of remotes is not defined, create a remote for all of the hosts that have a configuration file
        else: host_ids = introspection.remote_host_ids()

        # Loop over all the remote host ids
        for host_id in host_ids:

            # Create a remote SKIRT execution context
            remote = SkirtRemote()

            # Check whether a cluster name is defined
            cluster_name = self.cluster_names[host_id] if host_id in self.cluster_names else None

            # Setup the remote for the specified host
            remote.setup(host_id, cluster=cluster_name)

            # Add the remote to the list of remote objects
            self.remotes.append(remote)

    # -----------------------------------------------------------------

    def assign(self):

        """
        This function ...
        :return:
        """

        # Determine how many simulations are assigned to each remote
        queue_length = len(self.queue)
        quotient = queue_length // len(self.remotes)
        remainder = queue_length % len(self.remotes)
        assignment = []
        for i in range(len(self.remotes)):
            nvalues = quotient + 1 if i < remainder else quotient
            assignment.append(nvalues)

        # Set the assignment
        self.assignment = iter(assignment)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for the different remote hosts ...")

        # Loop over the different remote hosts
        for remote in self.remotes:

            # Check whether the parallelization has already been defined by the user for this remote host
            if remote.host_id in self.parallelization: continue

            # Debugging
            log.debug("Setting the parallelization scheme for host '" + remote.host_id + "' ...")

            # Get the number of cores per process as defined in the configuration
            cores_per_process = self.config.cores_per_process

            # Get the amount of (currently) free cores on the remote host
            cores = int(remote.free_cores)

            # Determine the number of thread to be used per core
            threads_per_core = remote.threads_per_core if remote.use_hyperthreading else 1

            # Create the parallelization object
            parallelization = Parallelization.from_free_cores(cores, cores_per_process, threads_per_core)

            # Debugging
            log.debug("Using " + str(parallelization.processes) + " processes and " + str(parallelization.threads) + " threads per process on this remote")
            log.debug("Parallelization scheme: " + str(parallelization))

            # Set the parallelization scheme for this host
            self.parallelization[remote.host_id] = parallelization

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

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # The remote input path
        remote_input_path = None

        # The complete list of simulations
        simulations = []

        # Loop over the different remotes
        for remote in self.remotes:

            # Get the parallelization scheme for this remote host
            parallelization = self.parallelization[remote.host_id]

            # Cache the simulation objects scheduled to the current remote
            simulations_remote = []

            # Repeat for a specified number of times
            for _ in range(next(self.assignment)):

                # Get the last item from the queue (it is removed)
                definition, name, parallelization_item, analysis_options_item = self.queue.pop()
                if parallelization_item is None: parallelization_item = parallelization

                # Check whether scheduling options are defined for this simulation and for this remote host
                if remote.host_id in self.scheduling_options and name in self.scheduling_options[remote.host_id]:
                    scheduling_options = self.scheduling_options[remote.host_id][name]
                else: scheduling_options = None

                # Check whether analysis options are specified
                if analysis_options_item is None:
                    logging_options = self.logging_options.copy()
                    analysis_options_item = self.create_analysis_options(definition, name, logging_options)
                else: logging_options = self.logging_options

                # Queue the simulation
                simulation = remote.add_to_queue(definition, logging_options, parallelization_item, name=name,
                                                 scheduling_options=scheduling_options, remote_input_path=remote_input_path,
                                                 analysis_options=analysis_options_item)
                simulations_remote.append(simulation)

                # Set the parallelization scheme of the simulation (important since SkirtRemote does not know whether
                # hyperthreading would be enabled if the user provided the parallelization_item when adding the
                # simulation to the queue
                simulation.parallelization = parallelization_item

                # If the input directory is shared between the different simulations
                if self.config.shared_input and remote_input_path is None: remote_input_path = simulation.remote_input_path

                ## SET OPTIONS

                # Remove remote files
                simulation.remove_remote_input = not self.config.keep and not self.config.shared_input
                simulation.remove_remote_output = not self.config.keep
                simulation.remove_remote_simulation_directory = not self.config.keep and not self.config.shared_input

                # Retrieval
                simulation.retrieve_types = self.config.retrieve_types

                ##

                # Save the simulation object
                simulation.save()

            # Set a path for the script file to be saved to locally (for manual inspection)
            if remote.host_id in self.script_paths:
                local_script_path = fs.join(self.script_paths[remote.host_id], time.unique_name() + ".sh")
            else: local_script_path = None

            # Determine screen name
            screen_name = time.unique_name("batch_launcher")

            # Set a path for the screen output to be saved remotely (for debugging)
            if remote.host_id in self.save_screen_output:
                remote_skirt_dir_path = remote.skirt_dir
                remote_skirt_screen_output_path = fs.join(remote_skirt_dir_path, "screen output")
                if not remote.is_directory(remote_skirt_screen_output_path): remote.create_directory(remote_skirt_screen_output_path)
                this_screen_output_path = fs.join(remote_skirt_screen_output_path, screen_name)
                remote.create_directory(this_screen_output_path)
                screen_output_path = this_screen_output_path
            else: screen_output_path = None

            # Start the queue
            remote.start_queue(screen_name=screen_name, group_simulations=self.config.group_simulations, local_script_path=local_script_path, screen_output_path=screen_output_path)

            # Set the screen name for all of the simulation objects
            for simulation in simulations_remote:

                simulation.screen_name = screen_name
                simulation.remote_screen_output_path = screen_output_path
                simulation.save()

            # Add the simulations of this remote to the total list of simulations
            simulations += simulations_remote

        # Return the list of simulations
        return simulations

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
            self.analyser.run(simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def create_analysis_options(self, definition, simulation_name, logging_options):

        """
        This function ...
        :param definition:
        :param simulation_name:
        :param logging_options:
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
        elif self.config.analysis.relative: extraction_path = fs.join(definition.base_path, self.config.analysis.extraction.path)
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
        elif self.config.analysis.relative: plotting_path = fs.join(definition.base_path, self.config.analysis.plotting.path)
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
        elif self.config.analysis.relative: misc_path = fs.join(definition.base_path, self.config.analysis.misc.path)
        else: misc_path = fs.join(self.config.analysis.misc.path, simulation_name)

        # Create the misc directory only if it is necessary
        if analysis_options.any_misc:
            if not fs.is_directory(misc_path): fs.create_directory(misc_path, recursive=True)
            analysis_options.misc.path = misc_path

        # Set the misc path to None otherwise
        else: analysis_options.misc.path = None

        # Set timing and memory table paths (if specified for this batch launcher)
        if self.config.timing_table_path is not None: analysis_options.timing_table_path = self.config.timing_table_path
        if self.config.memory_table_path is not None: analysis_options.memory_table_path = self.config.memory_table_path

        # Check the analysis options
        analysis_options.check(logging_options)

        # Return the analysis options
        return analysis_options

# -----------------------------------------------------------------
