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

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..tools import inspection, filesystem
from ..tools.logging import log
from ..basics.host import Host
from .parallelization import Parallelization

# -----------------------------------------------------------------

class BatchLauncher(Configurable):

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

        # The queue
        self.queue = []

        # The scheduling options for (some of) the simulations, accesed by keys that are the names given to the
        # simulations (see 'name' parameter of 'add_to_queue')
        self.scheduling_options = dict()

        # The assignment from items in the queue to the different remote hosts
        self.assignment = None

        # The parallelization scheme for the different remote hosts
        self.parallelization = dict()

        # The simulations that have been retrieved
        self.simulations = []

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new BatchLauncher instance
        launcher = cls()

        # Return the new batch launcher
        return launcher

    # -----------------------------------------------------------------

    def add_to_queue(self, arguments, name=None):

        """
        This function ...
        :param arguments:
        :param name: a name that is given to the simulation
        :return:
        """

        # Add the SkirtArguments object to the queue
        self.queue.append((arguments, name))

    # -----------------------------------------------------------------

    def set_scheduling_options(self, name, options):

        """
        This function ...
        :param name:
        :param options:
        :return:
        """

        self.scheduling_options[name] = options

    # -----------------------------------------------------------------

    @property
    def host_ids(self):

        """
        This function returns the IDs of the hosts that will be used by the batch launcher
        :return:
        """

        # If a list of remotes is defined
        if self.config.remotes is not None: host_ids = self.config.remotes

        # If a list of remotes is not defined, create a remote for all of the hosts that have a configuration file
        else: host_ids = inspection.remote_host_ids()

        # Return the list of host IDs
        return host_ids

    # -----------------------------------------------------------------

    @property
    def hosts(self):

        """
        This function returns the Host objects for all the hosts that will be used by the batch launcher
        :return:
        """

        # Create the list of hosts
        hosts = []
        for host_id in self.host_ids: hosts.append(Host(host_id))

        # Return the list of hosts
        return hosts

    # -----------------------------------------------------------------

    @property
    def scheduler_host_ids(self):

        """
        This function returns just the IDS of the hosts which use a scheduling system. This is useful for when example
        the user wants to specify the parallelization scheme only for these hosts, while the parallelization strategy
        for the other hosts can be left up to the BatchLauncher class based on the current load of that system.
        :return:
        """

        # Initialize a list to contain the host IDs
        host_ids = []

        # Loop over the IDs of all the hosts used by the BatchLauncher
        for id in self.host_ids:

            # Create Host instance
            host = Host(id)

            # If it's a scheduler, add it to the list
            if host.scheduler: host_ids.append(id)

        # Return the list of hosts which use a scheduling system
        return host_ids

    # -----------------------------------------------------------------

    @property
    def scheduler_hosts(self):

        """
        This function is similar to 'scheduler_host_ids' but returns the Host objects instead of just the IDs of these
        hosts.
        :return:
        """

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

        # Clear the queue
        self.queue = []

        # Clear the scheduling options
        self.scheduling_options = dict()

        # Clear the assignment
        self.assignment = None

        # Clear the simulations
        self.simulations = []

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BatchLauncher, self).setup()

        # If a list of remotes is defined
        if self.config.remotes is not None: host_ids = self.config.remotes

        # If a list of remotes is not defined, create a remote for all of the hosts that have a configuration file
        else: host_ids = inspection.remote_host_ids()

        # Loop over all the remote host ids
        for host_id in host_ids:

            # Create a remote SKIRT execution context
            remote = SkirtRemote()

            # Setup the remote for the specified host
            remote.setup(host_id)

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
            cores = remote.free_cores

            # Using the number of cores per process defined in the configuration, determine the number of processes
            processes = int(cores / cores_per_process)

            # Create the parallelization object
            parallelization = Parallelization()

            # Set the parallelization properties
            parallelization.cores = cores
            parallelization.threads_per_core = remote.threads_per_core if remote.use_hyperthreading else 1
            parallelization.processes = processes

            # Debugging
            log.debug("Using " + str(parallelization.processes) + " processes and " + str(parallelization.threads) + " threads per process on this remote")

            # Set the parallelization scheme for this host
            self.parallelization[remote.host_id] = parallelization

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # The remote input path
        remote_input_path = None

        # The complete list of simulations
        simulations = []

        # Loop over the different remotes
        for remote in self.remotes:

            # Get the parallelization scheme for this remote host
            parallellization = self.parallelization[remote.host_id]
            processes = parallellization.processes
            threads = parallellization.threads

            # Cache the simulation objects scheduled to the current remote
            simulations_remote = []

            # Repeat for a specified number of times
            for _ in range(next(self.assignment)):

                # Get the last item from the queue (it is removed)
                arguments, name = self.queue.pop()

                # Set the parallelization
                arguments.parallel.processes = processes
                arguments.parallel.threads = threads

                # Check whether scheduling options are defined for this simulation
                scheduling_options = self.scheduling_options[name] if name in self.scheduling_options else None

                # Queue the simulation
                simulation = remote.add_to_queue(arguments, name=name, scheduling_options=scheduling_options, remote_input_path=remote_input_path)
                simulations_remote.append(simulation)

                # If the input directory is shared between the different simulations,
                if self.config.shared_input and remote_input_path is None: remote_input_path = simulation.remote_input_path

                # Set the parallelization properties for the simulation (this is actually already done by SkirtRemote in the add_to_queue function)
                simulation.parallelization = parallellization

                # Set the analysis options for the simulation
                self.set_analysis_options(simulation)

                # Save the simulation object
                simulation.save()

            # Start the queue if the remote host does not use a scheduling system
            if not remote.scheduler:

                # Start the queue
                screen_name = remote.start_queue()

                # Set the screen name for all of the simulation objects
                for simulation in simulations_remote:
                    simulation.screen_name = screen_name
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

    def set_analysis_options(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Set the options
        simulation.set_analysis_options(self.config.analysis)

        # Determine the extraction directory for this simulation (and create it if necessary)
        if self.config.analysis.extraction.path is not None: extraction_path = filesystem.join(self.config.analysis.extraction.path, simulation.name)
        else: extraction_path = filesystem.join(simulation.output_path, "extr")
        if simulation.analysis.any_extraction:
            if not filesystem.is_directory(extraction_path): filesystem.create_directory(extraction_path, recursive=True)
            simulation.analysis.extraction.path = extraction_path

        # Determine the plotting directory for this simulation (and create it if necessary)
        if self.config.analysis.plotting.path is not None: plotting_path = filesystem.join(self.config.analysis.plotting.path, simulation.name)
        else: plotting_path = filesystem.join(simulation.output_path, "plot")
        if simulation.analysis.any_plotting:
            if not filesystem.is_directory(plotting_path): filesystem.create_directory(plotting_path, recursive=True)
            simulation.analysis.plotting.path = plotting_path

        # Determine the 'misc' directory for this simulation (and create it if necessary)
        if self.config.analysis.misc.path is not None: misc_path = filesystem.join(self.config.analysis.misc.path, simulation.name)
        else: misc_path = filesystem.join(simulation.output_path, "misc")
        if simulation.analysis.any_misc:
            if not filesystem.is_directory(misc_path): filesystem.create_directory(misc_path, recursive=True)
            simulation.analysis.misc.path = misc_path

        # Remove remote files
        simulation.remove_remote_input = not self.config.keep and not self.config.shared_input
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep and not self.config.shared_input

        # Retrieval
        simulation.retrieve_types = self.config.retrieve_types

# -----------------------------------------------------------------
