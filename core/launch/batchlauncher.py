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

        # The assignment from items in the queue to the different remotes
        self.assignment = None

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

    def add_to_queue(self, arguments, scheduling_options=None):

        """
        This function ...
        :param arguments:
        :param scheduling_options: can for example specify the walltime ...
        :return:
        """

        # Add the SkirtArguments object to the queue
        self.queue.append((arguments, scheduling_options))

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

        # 2. Launch the simulations
        simulations = self.simulate()

        # 3. Retrieve the simulations that are finished
        self.retrieve()

        # 4. Analyse the output of the retrieved simulations
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

    def simulate(self):

        """
        This function ...
        :return:
        """

        # The complete list of simulations
        simulations = []

        # Loop over the different remotes
        for remote in self.remotes:

            # Cache the simulation objects scheduled to the current remote
            simulations_remote = []

            # Repeat for a specified number of times
            for _ in range(next(self.assignment)):

                # Get the last item from the queue (it is removed)
                arguments, scheduling_options = self.queue.pop()

                # Queue the simulation
                simulation = remote.add_to_queue(arguments, scheduling_options=scheduling_options)
                simulations_remote.append(simulation)

                # Add additional information to the simulation object
                self.add_analysis_info(simulation)

                # Save the simulation object
                simulation.save()

            # Start the queue if the remote host does not use a scheduling system
            if remote.scheduler:

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
        log.info("Analysing the output of retrieved simulations...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def add_analysis_info(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Extraction
        simulation.extract_progress = self.config.extraction.progress
        simulation.extract_timeline = self.config.extraction.timeline
        simulation.extract_memory = self.config.extraction.memory

        # Determine the extraction directory for this simulation
        if self.config.extraction.path is not None: extraction_path = filesystem.join(self.config.extraction.path, simulation.id)
        else: extraction_path = filesystem.join(simulation.output_path, "extract")
        simulation.extraction_path = extraction_path

        # Plotting
        simulation.plot_progress = self.config.plotting.progress
        simulation.plot_timeline = self.config.plotting.timeline
        simulation.plot_memory = self.config.plotting.memory
        simulation.plot_seds = self.config.plotting.seds
        simulation.plot_grids = self.config.plotting.grids

        # Determine the plotting directory for this simulation
        if self.config.plotting.path is not None: plotting_path = filesystem.join(self.config.plotting.path, simulation.id)
        else: plotting_path = filesystem.join(simulation.output_path, "plot")
        simulation.plotting_path = plotting_path

        # Advanced
        simulation.make_rgb = self.config.advanced.rgb
        simulation.make_wave = self.config.advanced.wavemovie

        # Remove remote files
        simulation.remove_remote_input = not self.config.keep and not self.config.shared_input
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep and not self.config.shared_input

        # Retrieval
        simulation.retrieve_types = self.config.retrieve_types

# -----------------------------------------------------------------
