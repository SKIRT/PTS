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
import os

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..tools import inspection, configuration
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
        log.info("Clearing the batch launcher...")

# -----------------------------------------------------------------
