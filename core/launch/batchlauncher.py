#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.batchlauncher Contains the BatchLauncher class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..tools import inspection, configuration

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

        ## Attributes

        # Initialize a list to contain different SkirtRemote instances for the different remote hosts
        self.remotes = []

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

        # Logging
        if arguments.debug: launcher.config.logging.level = "DEBUG"

        # Return the new batch launcher
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Set the parallelization scheme
        self.set_parallelization()

        # 3. Run the simulation
        self.simulate()

        # 4. Retrieve the simulations that are finished
        self.retrieve()

        # 5. Analyse the output of the retrieved simulations
        self.analyse()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BatchLauncher, self).setup()

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

            # If a list of remotes is defined and this remote is not in it, skip it
            if self.config.remote is not None and host_id not in self.config.remote: continue

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

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the batch launcher...")

# -----------------------------------------------------------------
