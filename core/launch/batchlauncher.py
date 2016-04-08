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

        # The assignment from items in the queue to the different remotes
        self.assignment = None

        # The parallelization scheme for the different remotes
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

        ## Adjust the configuration settings according to the command-line arguments

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

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Determine how many simulations are assigned to each remote
        self.assign()

        # 3. Set the parallelization scheme for the different remotes
        self.set_parallelization()

        # 4. Launch the simulations
        simulations = self.simulate()

        # 5. Retrieve the simulations that are finished
        self.retrieve()

        # 6. Analyse the output of the retrieved simulations
        self.analyse()

        # Return the simulations that are just scheduled
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

            # Debugging
            log.debug("Setting the parallelization scheme for host '" + remote.host_id + "' ...")

            # Fix the number of processes to 16
            processes = 16

            # Calculate the maximum number of threads per process based on the current cpu load of the system
            threads = int(remote.free_cores / processes)

            # If hyperthreading should be used for the remote host, we can even use more threads
            if remote.use_hyperthreading: threads *= remote.threads_per_core

            # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
            if threads < 1:

                processes = int(remote.free_cores)
                threads = 1

            # Debugging
            log.debug("Using " + str(processes) + " processes and " + str(threads) + " threads per process on this remote")

            # Set the parallelization scheme for this host
            self.parallelization[remote.host_id] = (processes, threads)

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
            processes, threads = self.parallelization[remote.host_id]

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
