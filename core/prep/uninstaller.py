#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.uninstall Uninstall PTS and/or SKIRT on remote hosts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from pts.core.tools.logging import log
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

class Uninstaller(Configurable):
        
    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(Uninstaller, self).__init__(config)

        # The remotes
        self.remotes = []

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Locate SKIRT directory
        if "skirt" in self.config.skirt_and_or_pts: self.uninstall_skirt()

        # Locate PTS directory
        if "pts" in self.config.skirt_and_or_pts: self.uninstall_pts()

        # Uninstall conda
        if self.config.conda: self.uninstall_conda()

        # Uninstall Qt
        if self.config.qt: self.uninstall_qt()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Uninstaller, self).setup(**kwargs)

        # Loop over the remotes
        for host_id in self.config.remotes:

            # Create remote
            remote = Remote()
            if not remote.setup(host_id):
                log.warning("Remote host '" + host_id + "' is offline")
                continue
            self.remotes.append(remote)

    # -----------------------------------------------------------------

    def uninstall_skirt(self):

        """
        This function ...
        :return:
        """

        # Loop over the remotes
        for remote in self.remotes:

            skirt_root_path = remote.skirt_root_path
            if not remote.is_directory(skirt_root_path): return

            # Inform the user
            log.info("Removing SKIRT from the remote host '" + remote.host_id + "' ...")

            # Remove the entire directory
            remote.remove_directory(skirt_root_path)

            # Remove lines from shell configuration file
            comment = "For SKIRT and FitSKIRT, added by PTS (Python Toolkit for SKIRT)"
            # terminal.remove_aliases_and_variables_with_comment(comment)
            remote.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_pts(self):

        """
        This function ...
        :return:
        """

        for remote in self.remotes:

            pts_root_path = remote.pts_root_path
            if not remote.is_directory(pts_root_path): return

            # Inform the user
            log.info("Removing PTS from the remote host '" + remote.host_id + "' ...")

            # Remove the entire directory
            remote.remove_directory(pts_root_path)

            # Remove lines from shell configuration file
            comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
            # terminal.remove_aliases_and_variables_with_comment(comment)
            remote.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_conda(self):

        """
        This function ...
        :return:
        """

        for remote in self.remotes:

            # Determine path of miniconda installation
            installation_path = fs.join(remote.home_directory, "miniconda")
            if not remote.is_directory(installation_path): return

            # Inform the user
            log.info("Removing Conda from the remote host '" + remote.host_id + "' ...")

            # Remove the directory
            remote.remove_directory(installation_path)

            # Remove lines from shell configuration file
            comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
            remote.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_qt(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
