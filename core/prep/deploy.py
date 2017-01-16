#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.deploy Contains the Deployer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools.logging import log
from ..tools import introspection
from .update import SKIRTUpdater, PTSUpdater
from .installation import SKIRTInstaller, PTSInstaller
from ..remote.versionchecker import VersionChecker
from ..remote.configurable import RemotesConfigurable

# -----------------------------------------------------------------

class Deployer(RemotesConfigurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(Deployer, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Install and update SKIRT
        if self.config.skirt: self.install_or_update_skirt()

        # 3. Install and update PTS
        if self.config.pts: self.install_or_update_pts()

        # 4. Check versions
        if self.config.check: self.check_versions()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Deployer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def install_or_update_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing or updating SKIRT ...")

        # Locally
        if self.config.local: self.install_or_update_skirt_locally()

        # Remotely
        self.install_or_update_skirt_remotely()

    # -----------------------------------------------------------------

    def install_or_update_skirt_locally(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing or updating SKIRT locally ...")

        # Check whether installed
        installed = introspection.skirt_is_present()

        # Not installed: install
        if not installed:

            # Create installer and run it
            installer = SKIRTInstaller()
            installer.run()

        # Installed and not SKIRT developer: update
        elif not introspection.is_skirt_developer():

            # Create updater and run it
            updater = SKIRTUpdater()
            updater.run()

    # -----------------------------------------------------------------

    def install_or_update_skirt_remotely(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing or updating SKIRT remotely ...")

        # Loop over the remotes
        for remote in self.remotes:

            # Check if SKIRT is present
            installed = remote.has_skirt

            # If installed, update
            if installed:

                # Debugging
                log.debug("SKIRT is present on host '" + remote.host_id + "', updating ...")

                # Create the updater
                updater = SKIRTUpdater()

                # Run the updater
                updater.run(remote=remote)

            # If not installed, install
            else:

                # Debugging
                log.debug("SKIRT is not present on host '" + remote.host_id + "', installing ...")

                # Create the installer
                installer = SKIRTInstaller()

                # Set the remote host ID
                installer.config.force = True     # SKIRT could not be found as an executable, thus remove whatever partial SKIRT installation there is
                installer.config.repository = "origin"

                # Run the installer
                installer.run(remote=remote)

    # -----------------------------------------------------------------

    def install_or_update_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing or updating PTS ...")

        # Locally
        if self.config.local and not introspection.is_pts_developer(): self.update_pts_locally()

        # Remotely
        self.install_or_update_pts_remotely()

    # -----------------------------------------------------------------

    def update_pts_locally(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating PTS locally ...")

        # Create updater and run it
        updater = PTSUpdater()
        updater.run()

    # -----------------------------------------------------------------

    def install_or_update_pts_remotely(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing or updating PTS remotely ...")

        # Loop over the remotes
        for remote in self.remotes:

            # Do not install PTS on some hosts
            if self.config.pts_on is not None and remote.host_id not in self.config.pts_on: continue

            # Debugging
            log.debug("Checking the presence of PTS on remote host '" + remote.host_id + "' ...")

            # Check whether PTS is present
            installed = remote.has_pts

            # If installed, update
            if installed:

                # Debugging
                log.debug("PTS is present on host '" + remote.host_id + "', updating ...")

                # Create the updater
                updater = PTSUpdater()

                # Run the updater
                updater.run(remote=remote)

            # If not installed, install
            else:

                # Debugging
                log.debug("PTS is not present on host '" + remote.host_id + "', installing ...")

                # Create installer
                installer = PTSInstaller()

                # Install
                installer.run(remote=remote)

    # -----------------------------------------------------------------

    def check_versions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking versions after deployment ...")

        # Initialize the version checker
        checker = VersionChecker()
        checker.config.show = False

        # Run the checker
        checker.run(remotes=self.remotes)

        # C++ compiler
        log.info("Local C++ compiler version: " + introspection.cpp_compiler_version())
        log.info("Remote C++ compiler versions:")
        for host_id in checker.cpp_versions: log.info(" - " + host_id + ": " + checker.cpp_versions[host_id])

        # MPI compiler
        if introspection.has_mpi_compiler():
            log.info("Local MPI compiler version: " + introspection.mpi_compiler_version())
            log.info("Remote MPI compiler versions:")
            for host_id in checker.mpi_versions: log.info(" - " + host_id + ": " + checker.mpi_versions[host_id])

        # Qt version
        log.info("Local Qt version: " + introspection.qmake_version())
        log.info("Remote Qt versions:")
        for host_id in checker.qt_versions: log.info(" - " + host_id + ": " + checker.qt_versions[host_id])

        # Python
        log.info("Local Python version: " + introspection.python_version_long())
        log.info("Remote Python versions:")
        for host_id in checker.python_versions: log.info(" - " + host_id + ": " + checker.python_versions[host_id])

        # SKIRT
        log.info("Local SKIRT version: " + introspection.skirt_version())
        log.info("Remote SKIRT versions:")
        for host_id in checker.skirt_versions: log.info(" - " + host_id + ": " + checker.skirt_versions[host_id])

        # PTS
        log.info("Local PTS version: " + introspection.pts_version())
        log.info("Remote PTS versions:")
        for host_id in checker.pts_versions: log.info(" - " + host_id + ": " + checker.pts_versions[host_id])

# -----------------------------------------------------------------
