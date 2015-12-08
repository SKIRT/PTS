#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module can be used to synchronize remote SKIRT simulations
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from .analyser import SimulationAnalyser
from ..test.scalinganalyser import ScalingAnalyser
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..tools import inspection, configuration

# -----------------------------------------------------------------

class format:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# -----------------------------------------------------------------

class RemoteSynchronizer(Configurable):

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
        super(RemoteSynchronizer, self).__init__(config)

        ## Attributes

        # Initialize a list to contain different SkirtRemote instances for the different remote hosts
        self.remotes = []

        # The simulation results analyser
        self.analyser = SimulationAnalyser()

        # The scaling results analyser
        self.scalinganalyser = ScalingAnalyser()

        # Initialize a list to contain the retreived simulations
        self.simulations = []

        # Set the delete list to None initially
        self.delete = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new RemoteSynchronizer instance
        synchronizer = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Logging
        if arguments.debug: synchronizer.config.logging.level = "DEBUG"

        # Set the remote name and the delete dictionary
        synchronizer.config.remote = arguments.remote
        synchronizer.config.delete = arguments.delete

        # Return the new synchronizer
        return synchronizer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 1. Extract information from the simulation's log files
        self.retreive()

        # 2. Analyse
        self.analyse()

        # 3. Announce the status of the simulations
        self.announce()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(RemoteSynchronizer, self).setup()

        # Search for files that define remote host configurations
        hosts_directory = os.path.join(inspection.pts_user_dir, "hosts")
        if not os.path.isdir(hosts_directory): os.makedirs(hosts_directory)

        # If the hosts directory is empty, place a template host configuration file there and exit with an error
        if len([item for item in os.listdir(hosts_directory) if os.path.isfile(os.path.join(hosts_directory, item))]) == 0:
            config = configuration.new()

            config.name = "server.institute.com"
            config.user = "user000"
            config.password = None
            config.output_path = "~/DATA/SKIRT"
            config.scheduler = True
            config.mpi_command = "mpirun"
            config.modules = ["examplemodule/2016/version2", "examplemodule2/2016/version0.1.3"]
            config.clusters.default = "cluster_a"
            config.clusters.cluster_a.cores = 16
            config.clusters.cluster_b.cores = 32

            config_file_path = os.path.join(hosts_directory, "template.cfg")
            config_file = open(config_file_path, 'w')
            config.save(config_file)

            self.log.error("No remote configuration files were found. Placing a template into PTS/user/hosts. Adjust it"
                           " for the remote hosts you want to use before using 'pts launch' or 'pts status'.")
            exit()

        # Loop over the configuration files in the hosts directory
        for filename in os.listdir(hosts_directory):

            # Skip the template configuration file
            if filename == "template.cfg": continue

            # Determine the full path to the host file
            file_path = os.path.join(hosts_directory, filename)

            # Ignore directories and hidden files
            if filename.startswith(".") or not os.path.isfile(file_path): continue

            # Get the host id for this line
            host_id = filename.split(".")[0]

            # Check whether there are simulation files corresponding to this host ID
            host_run_dir = os.path.join(inspection.skirt_run_dir, host_id)

            # If there are no simulation files for this host, skip it
            if len([item for item in os.listdir(host_run_dir) if item.endswith(".sim") and not item.startswith(".")]) == 0: continue

            # Create a remote SKIRT execution context
            remote = SkirtRemote()

            # Setup the remote execution context
            remote.setup(host_id)

            # Add the remote to the list of remote objects
            self.remotes.append(remote)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the synchronizer")

        # Set default values for attributes
        self.simulations = []
        self.delete = None

    # -----------------------------------------------------------------

    def retreive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Retreiving the output of finished simulations")

        # Loop over the different remotes
        for remote in self.remotes:

            # Retreive simulations
            self.simulations += remote.retreive()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Analysing the output of retreived simulations")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation)

            # If this simulation is part of a scaling test, run the scalinganalyser
            if simulation.scaling_run_name is not None:

                # Run the scaling analyser and clear it afterwards
                self.scalinganalyser.run(simulation, self.analyser.timeline, self.analyser.memory)
                self.scalinganalyser.clear()

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def announce(self):

        """
        This function ...
        :return:
        """

        # Loop over the different remotes
        for remote in self.remotes:

            # Get the status of the different simulations
            status = remote.status

            # Show the name of the current remote
            if len(status) > 0: self.log.info("Simulations on remote " + remote.host_name + ":")
            print()

            # Get the status of the different simulations
            for entry in status:

                # The path to the simulation file
                path = entry.file_path

                # Get the simulation rank
                simulation_rank = os.path.basename(path).split(".")[0]

                # The simulation name
                simulation_name = entry.name

                # The ski file path
                ski_path = entry.ski_path

                # The remote output path
                remote_output_path = entry.remote_output_path

                # The simulation status
                simulation_status = entry.status

                # If this simulation has already finished, check whether the results of this simulation have been retreived
                if simulation_status == "finished":
                    simulation_file = open(path, 'r')
                    lines = simulation_file.readlines()
                    last = lines[len(lines)-1]
                    if "retreived at" in last: simulation_status = "retreived"

                if self.config.delete is not None and entry.id in self.config.delete[remote.config.host_id]:

                    print(format.FAIL + "  [ X ] " + simulation_name + ": " + simulation_status + format.END)

                else:

                    # Show the status of the current simulation
                    print(format.GREEN + "  [" + str(simulation_rank) + "] " + simulation_name + ": " + simulation_status + format.END)

        print()

        """
        self.log.info(scalingtest + ":")

        for simulation in simulations:

            run = simulation[0]
            sim = simulation[1]
            status = simulation[2]
            rank = simulation[3]
            jobfilepath = simulation[4]

            if delete is not None and rank in delete:

                tag = "    [X] "
                deletethisfile = True

            else:

                tag = "    [" + str(rank) + "] "
                deletethisfile = False

            if status == "queued":

                # TODO: figure out what to do when queued jobs should be deleted

                self.log.info(tag + run + " > " + sim + ": " + status)

            elif status == "running":

                # TODO: figure out what to do when running jobs should be deleted

                self.log.warning(tag + run + " > " + sim + ": " + status)

            elif status == "finished":

                # Delete the job file if requested
                if deletethisfile: os.remove(jobfilepath)

                self.log.success(tag + run + " > " + sim + ": " + status)

            elif status == "crashed":

                # Delete the job file if requested
                if deletethisfile: os.remove(jobfilepath)

                self.log.failure(tag + run + " > " + sim + ": " + status)

            elif status == "unknown":

                # TODO: figure out what to do when jobs with unknown status should be deleted

                self.log.failure(tag + run + " > " + sim + ": " + status)
        """

# -----------------------------------------------------------------
