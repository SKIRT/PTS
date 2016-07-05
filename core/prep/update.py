#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.update Contains the SkirtUpdater class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import OldConfigurable
from ..simulation.execute import SkirtExec
from ..simulation.remote import SkirtRemote

# -----------------------------------------------------------------

class SkirtUpdater(OldConfigurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """
        
        # Call the constructor of the base class
        super(SkirtUpdater, self).__init__(config, "core")

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()
        
    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :return:
        """

        # Create a new SkirtUpdater instance
        updater = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Remote host
        updater.config.remote = arguments.remote

        # Return the new SkirtUpdater instance
        return updater

    # -----------------------------------------------------------------

    def run(self):
        
        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2. Update
        self.update()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtUpdater, self).setup()

        # Setup the remote execution environment if necessary
        if self.config.remote is not None: self.remote.setup(self.config.remote)

    # -----------------------------------------------------------------

    def update(self):

        """
        This function ...
        :return:
        """

        # Update the remote SKIRT executable
        if self.config.remote is not None: self.remote.update()

        # Update the local SKIRT executable
        else: self.skirt.update()

# -----------------------------------------------------------------
