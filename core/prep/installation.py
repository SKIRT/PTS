#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.installation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import urllib

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.execute import SkirtExec
from ..simulation.remote import SkirtRemote

# -----------------------------------------------------------------

class SkirtInstaller(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SkirtInstaller, self).__init__(config)

        ## Attributes

        # Create the SKIRT execution context
        self.skirt = SkirtExec()

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        #self.has_qt = False
        
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

        # Logging
        if arguments.debug: updater.config.logging.level = "DEBUG"

        # Remote host
        updater.config.remote = arguments.remote

        # Use the private repository
        updater.config.private = arguments.private

        # Return the new SkirtUpdater instance
        return updater

    # -----------------------------------------------------------------

    def run(self):
        
        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2.
        #self.check_qt()

        # 3.
        #if not self.has_qt: self.install_qt()

        # 4.
        #self.get()

        # 5.
        self.install()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtInstaller, self).setup()

        # Setup the remote execution environment if necessary
        if self.config.remote is not None: self.remote.setup(self.config.remote, pre_installation=True)

    # -----------------------------------------------------------------

    def check_qt(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def install_qt(self):

        """
        This function ...
        :return:
        """

        link = "http://download.qt.io/official_releases/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"

        # Determine the path for the Qt source code
        path = os.path.join(None, "qt.tar.gz")

        # Create a
        urllib.urlretrieve(self.config.qt_link, path)

    # -----------------------------------------------------------------

    def get(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def install(self):

        """
        This function ...
        :return:
        """

        # If the installation is to be performed on a remote system
        if self.config.remote is not None: self.remote.install(self.config.private)

        # If SKIRT is to be installed on this system
        else: self.skirt.install(self.config.private)

# -----------------------------------------------------------------
