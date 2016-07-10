#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.installation Contains the SkirtInstaller class.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import urllib

# Import the relevant PTS classes and modules
from ..basics.configurable import OldConfigurable
from ..simulation.execute import SkirtExec
from ..simulation.remote import Remote
from ..tools import introspection
from ..tools import filesystem as fs
from ..tools.logging import log

# -----------------------------------------------------------------

class SkirtInstaller(OldConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SkirtInstaller, self).__init__(config, "core")

        # -- Attributes --

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. If Qt is not present, install it
        if self.qmake_path is None: self.install_qt()

        # 3. Get the SKIRT source code
        self.get()

        # 4. Compile SKIRT
        self.install()

        # 5. Test the installation
        self.test()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Translated from ./makeSKIRT.sh

        # Get a list of qmake paths installed on this system
        qmake_paths = []

        for qt_dir in fs.directories_in_path(fs.home(), startswith="Qt"):
            qmake_paths = fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

        for qt_dir in fs.directories_in_path("/usr/local", startswith="Qt"):
            qmake_paths += fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

        qmake_path = introspection.qmake_path()
        qmake_paths += [qmake_path] if qmake_path is not None else []

        # Get the most recent installation

        # Check whether the Qt version is supported
        # if [[ $VERSION > '5.2.0' ]]



    # -----------------------------------------------------------------

    def install_qt(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the SKIRT source code ...")

    # -----------------------------------------------------------------

    def install(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing SKIRT ...")

    # -----------------------------------------------------------------

    def test(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation ...")

        # Create the SKIRT execution context
        skirt = SkirtExec()

# -----------------------------------------------------------------

class SkirtRemoteInstaller(OldConfigurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SkirtRemoteInstaller, self).__init__(config, "core")

        # -- Attributes --

        # Create the remote execution context
        self.remote = Remote()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SkirtUpdater instance
        updater = cls()

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
        if self.qmake_path is None: self.install_qt()

        # 3.
        self.get()

        # 4.
        self.install()

        # 5.
        self.test()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemoteInstaller, self).setup()

        # Setup the remote execution environment
        self.remote.setup(self.config.remote)

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
        path = fs.join(None, "qt.tar.gz")

        # Create a
        urllib.urlretrieve(self.config.qt_link, path)

        # FOR HPC: module load Qt/5.2.1-intel-2015a

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
