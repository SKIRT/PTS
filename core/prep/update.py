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

# Import standard modules
import subprocess
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..basics.remote import Remote
from ..tools.logging import log
from ..tools import introspection

# -----------------------------------------------------------------

class Updater(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(Updater, self).__init__(config)

        # The remote execution environment
        self.remote = None

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
        super(Updater, self).setup()

        # Setup the remote execution environment if necessary
        if self.config.remote is not None:

            # Create and setup the remote execution environment
            self.remote = Remote()
            self.remote.setup(self.config.remote)

    # -----------------------------------------------------------------

    def update(self):

        """
        This function ...
        :return:
        """

        # Update locally or remotely
        if self.remote is None: self.update_local()
        else: self.update_remote()

    # -----------------------------------------------------------------

    @abstractmethod
    def update_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def update_remote(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

class SKIRTUpdater(Updater):
    
    """
    This class ...
    """

    def update_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating SKIRT locally ...")

        # Get the SKIRT repo directory
        skirt_repo_path = introspection.skirt_repo_dir

        # Debugging
        log.debug("Getting latest version ...")

        # Call the appropriate git command at the SKIRT repository directory
        subprocess.call(["git", "pull", "origin", "master"], cwd=skirt_repo_path)

        # Debugging
        log.debug("Compiling latest version ...")

        # Recompile SKIRT
        subprocess.call(["sh", "makeSKIRT.sh"], cwd=skirt_repo_path)

    # -----------------------------------------------------------------

    def update_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating SKIRT remotely on host '" + self.config.remote + "' ...")

        # Debugging
        log.debug("Getting latest version ...")

        # Change working directory to repository directory
        skirt_git_path = self.remote.skirt_repo_path
        self.remote.change_cwd(skirt_git_path)

        # Git pull
        self.remote.ssh.sendline("git pull origin master")
        self.remote.ssh.expect(":")

        if "Enter passphrase for key" in self.remote.ssh.before:

            self.remote.execute(self.config.pubkey_password, show_output=True)

        else: self.remote.prompt()

        # Debugging
        log.debug("Compiling latest version ...")

        # Build SKIRT
        self.remote.execute("./makeSKIRT.sh", show_output=True, output=False)

# -----------------------------------------------------------------

class PTSUpdater(Updater):

    """
    This class ...
    """

    def update_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating PTS locally ...")

        # Get the PTS repo directory
        pts_git_path = introspection.pts_package_dir

        # Debugging
        log.debug("Getting latest version ...")

        # Call the appropriate git command at the PTS repository directory
        subprocess.call(["git", "pull", "origin", "master"], cwd=pts_git_path)

    # -----------------------------------------------------------------

    def update_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating PTS remotely on host '" + self.config.remote + "' ...")

        # Debugging
        log.debug("Getting latest version ...")

        # Change working directory to repository directory
        pts_git_path = self.remote.pts_package_path
        self.remote.change_cwd(pts_git_path)

        # Git pull
        self.remote.ssh.sendline("git pull origin master")
        self.remote.ssh.expect(":")

        if "Enter passphrase for key" in self.remote.ssh.before:

            self.remote.execute(self.config.pubkey_password, show_output=True)

        else: self.remote.prompt()

# -----------------------------------------------------------------
