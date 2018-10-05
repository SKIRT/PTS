#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.uninstall Uninstall PTS and/or SKIRT on a remote host.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from pts.core.basics.og import log
from pts.core.remote.remote import Remote
from pts.core.tools import filesystem as fs
from pts.core.tools import terminal
from pts.core.tools import introspection

# -----------------------------------------------------------------

class Uninstaller(Configurable):
        
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Uninstaller, self).__init__(*args, **kwargs)

        # The remote
        self.remote = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Uninstall SKIRT
        if "skirt" in self.config.skirt_and_or_pts: self.uninstall_skirt()

        # 3. Uninstall PTS
        if "pts" in self.config.skirt_and_or_pts: self.uninstall_pts()

        # 4. Uninstall conda
        if self.config.conda: self.uninstall_conda()

        # 5. Uninstall Qt
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

        # Check if remote instance is passed
        if "remote" in kwargs: self.remote = kwargs.pop("remote")
        elif self.config.remote is not None:

            remote = Remote()
            if not remote.setup(self.config.remote): raise RuntimeError("The remote host is not available")
            self.remote = remote

        else: log.warning('No remote host is specified, will be uninstalling locally')

    # -----------------------------------------------------------------

    def uninstall_skirt(self):

        """
        This function ...
        :return:
        """

        if self.remote is None: self.uninstall_skirt_local()
        else: self.uninstall_skirt_remote()

    # -----------------------------------------------------------------

    def uninstall_skirt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling SKIRT locally ...")

        # Check installation
        skirt_root_path = introspection.skirt_root_dir
        if not fs.is_directory(skirt_root_path):
            log.warning("SKIRT was not found locally")
            return

        # Debugging
        log.debug("Removing the SKIRT directory ...")

        # Remove the entire directory
        fs.remove_directory(skirt_root_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        terminal.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling SKIRT on the remote host ...")

        # Check installation
        skirt_root_path = self.remote.skirt_root_path
        if not self.remote.is_directory(skirt_root_path):
            log.warning("SKIRT was not found on the remote host")
            return

        # Debugging
        log.debug("Removing the SKIRT directory ...")

        # Remove the entire directory
        self.remote.remove_directory(skirt_root_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        self.remote.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_pts(self):

        """
        This function ...
        :return:
        """

        if self.remote is None: self.uninstall_pts_local()
        else: self.uninstall_pts_remote()

    # -----------------------------------------------------------------

    def uninstall_pts_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling PTS locally ...")

        # Check installation
        pts_root_path = introspection.pts_root_dir
        if not fs.is_directory(pts_root_path):
            log.warning("PTS could not be found locally (which is certainly weird) ...")
            return

        # Debugging
        log.debug("Removing the PTS directory ...")

        # Remove the entire directory
        fs.remove_directory(pts_root_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
        terminal.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_pts_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling PTS on the remote host ...")

        # Check installation
        pts_root_path = self.remote.pts_root_path
        if not self.remote.is_directory(pts_root_path):
            log.warning("PTS is not present on the remote host")
            return

        # Debugging
        log.debug("Removing the PTS directory ...")

        # Remove the entire directory
        self.remote.remove_directory(pts_root_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
        self.remote.remove_aliases_and_variables_with_comment(comment)

    # -----------------------------------------------------------------

    def uninstall_conda(self):

        """
        This function ...
        :return:
        """

        if self.remote is None: self.uninstall_conda_local()
        else: self.uninstall_conda_remote()

    # -----------------------------------------------------------------

    def uninstall_conda_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling the Conda python distribution locally ...")

        # Check installation
        installation_path = fs.join(fs.home, "miniconda")
        if not fs.is_directory(installation_path):
            log.warning("Conda was not found locally")
            return

        # Debugging
        log.debug("Removing the Conda directory ...")

        # Remove the directory
        fs.remove_directory(installation_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
        terminal.remove_aliases_and_variables_with_comment(comment)
        terminal.remove_from_path_variable_containing("miniconda/bin")

    # -----------------------------------------------------------------

    def uninstall_conda_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uninstalling the Conda python distribution on the remote host ...")

        # Determine path of miniconda installation
        installation_path = fs.join(self.remote.home_directory, "miniconda")
        if not self.remote.is_directory(installation_path):
            log.warning("Conda was not found on the remote host")
            return

        # Debugging
        log.debug("Removing the Conda directory ...")

        # Remove the directory
        self.remote.remove_directory(installation_path)

        # Debugging
        log.debug("Removing lines from shell configuration ...")

        # Remove lines from shell configuration file
        comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
        self.remote.remove_aliases_and_variables_with_comment(comment)
        self.remote.remove_from_path_variable_containing("miniconda/bin")

    # -----------------------------------------------------------------

    def uninstall_qt(self):

        """
        This function ...
        :return:
        """

        if self.remote is None: self.uninstall_qt_local()
        else: self.uninstall_qt_remote()

    # -----------------------------------------------------------------

    def uninstall_qt_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def uninstall_qt_remote(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
