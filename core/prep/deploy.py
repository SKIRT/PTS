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
from ..basics.log import log
from ..tools import introspection
from .update import SKIRTUpdater, PTSUpdater
from .installation import SKIRTInstaller, PTSInstaller
from ..remote.versionchecker import VersionChecker
from ..remote.configurable import RemotesConfigurable
from ..tools import git

# -----------------------------------------------------------------

class Deployer(RemotesConfigurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Deployer, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Install and update SKIRT
        if self.config.skirt: self.install_or_update_skirt()

        # 3. Install and update PTS
        if self.config.pts: self.install_or_update_pts()

        # 4. Check versions
        if self.config.check and self.has_remotes: self.check_versions()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set log_conda flag
        kwargs["log_conda"] = True

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
        if self.config.local and not introspection.is_skirt_developer(): self.install_or_update_skirt_locally()

        # Remotely
        if self.has_remotes: self.install_or_update_skirt_remotely()

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

        # ALready installed
        if installed:

            # Clean install
            if self.config.clean:

                # Debugging
                log.debug("Peforming a clean install of SKIRT ...")

                installer = SKIRTInstaller()
                installer.config.force = True
                installer.config.repository = "origin"
                installer.run()

            # Update, if not SKIRT developer
            else:

                # Debugging
                log.debug("SKIRT is present, updating ...")

                updater = SKIRTUpdater()
                updater.config.dependencies = self.config.update_dependencies
                updater.run()

        # Not installed or not found
        else:

            # Debugging
            log.debug("SKIRT is not found, installing ...")

            # Create installer and run it
            installer = SKIRTInstaller()
            installer.config.force = True
            installer.config.repository = "origin"
            installer.run()

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

                # Check origins
                local_origin_url = introspection.skirt_git_remote_url(self.config.skirt_repo_name, pubkey_password=self.config.pubkey_password)
                remote_origin_url = remote.skirt_git_remote_url("origin")

                # Decompose
                host_local, user_or_organization_local, repo_name_local, username_local, password_local = git.decompose_repo_url(local_origin_url)
                host_remote, user_or_organization_remote, repo_name_remote, username_remote, password_remote = git.decompose_repo_url(remote_origin_url)

                # Compose in simple https to compare whether equal
                local_simple = git.compose_https(host_local, user_or_organization_local, repo_name_local)
                remote_simple = git.compose_https(host_remote, user_or_organization_remote, repo_name_remote)

                # Reinstall because of wrong origin
                if local_simple != remote_simple:

                    # Debugging
                    log.debug("Local and remote SKIRT repositories do not have the same source")
                    log.debug("Local original repository for SKIRT: '" + local_simple)
                    log.debug("Remote original repository for SKIRT: '" + remote_simple)
                    log.debug("Reinstalling SKIRT on the remote host '" + remote.host_id + "' ...")

                    # Install
                    installer = SKIRTInstaller()
                    installer.config.force = True
                    installer.config.repository = self.config.skirt_repo_name
                    installer.run(remote=remote)

                # Clean install
                if self.config.clean:

                    # Debugging
                    log.debug("Performing clean install of SKIRT on remote host '" + remote.host_id + "' ...")

                    # Installer
                    installer = SKIRTInstaller()
                    installer.config.force = True
                    installer.config.repository = self.config.skirt_repo_name
                    installer.run(remote=remote)

                else:

                    # Debugging
                    log.debug("SKIRT is present on host '" + remote.host_id + "', updating ...")

                    # Create the updater
                    updater = SKIRTUpdater()
                    updater.config.dependencies = self.config.update_dependencies

                    # Run the updater
                    updater.run(remote=remote)

            # If not installed, install
            else:

                # Debugging
                log.debug("SKIRT is not found on host '" + remote.host_id + "', installing ...")

                # Create the installer
                installer = SKIRTInstaller()

                # Set the remote host ID
                installer.config.force = True   # SKIRT could not be found as an executable, thus remove whatever partial SKIRT installation there is
                installer.config.repository = self.config.skirt_repo_name

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
        if self.has_remotes: self.install_or_update_pts_remotely()

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
        updater.config.dependencies = self.config.update_dependencies
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
            if installed: self.update_pts_on_remote(remote)

            # If not installed, install
            else: self.install_pts_on_remote(remote)

    # -----------------------------------------------------------------

    def update_pts_on_remote(self, remote):

        """
        This function ...
        :param remote:
        :return:
        """

        # Debugging
        #log.debug("Updating PTS on remote host '" + remote.host_id + "' ...")

        # Check origins
        local_origin_url = introspection.pts_git_remote_url(self.config.pts_repo_name, pubkey_password=self.config.pubkey_password)
        remote_origin_url = remote.pts_git_remote_url("origin")

        # Decompose
        host_local, user_or_organization_local, repo_name_local, username_local, password_local = git.decompose_repo_url(local_origin_url)
        host_remote, user_or_organization_remote, repo_name_remote, username_remote, password_remote = git.decompose_repo_url(remote_origin_url)

        # Compose in simple https to compare whether equal
        local_simple = git.compose_https(host_local, user_or_organization_local, repo_name_local)
        remote_simple = git.compose_https(host_remote, user_or_organization_remote, repo_name_remote)

        # Reinstall because of wrong origin
        if local_simple != remote_simple:

            # Debugging
            log.debug("Local and remote PTS repositories do not have the same source")
            log.debug("Local original repository for PTS: '" + local_simple + "'")
            log.debug("Remote original repository for PTS: '" + remote_simple + "'")
            log.debug("Reinstalling PTS on the remote host '" + remote.host_id + "' ...")

            # Install
            installer = PTSInstaller()
            installer.config.force = True
            installer.config.repository = self.config.pts_repo_name
            installer.run(remote=remote)

        # Clean install
        elif self.config.clean:

            # Debugging
            log.debug("Performing a clean install of PTS on remote host '" + remote.host_id + "' ...")

            # Install
            installer = PTSInstaller()
            installer.config.force = True
            installer.config.repository = self.config.pts_repo_name
            installer.run(remote=remote)

        else:

            # Debugging
            log.debug("PTS is present on host '" + remote.host_id + "', updating ...")

            # Create the updater
            updater = PTSUpdater()
            updater.config.dependencies = self.config.update_dependencies

            # Run the updater
            updater.run(remote=remote)

    # -----------------------------------------------------------------

    def install_pts_on_remote(self, remote):

        """
        This function ...
        :param remote:
        :return:
        """

        # Debugging
        log.debug("PTS is not present on host '" + remote.host_id + "', installing ...")

        # Create installer
        installer = PTSInstaller()
        installer.config.force = True
        installer.config.repository = self.config.pts_repo_name

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
        if self.nremotes == 1 and self.single_host_id in checker.cpp_versions: log.info("Remote C++ compiler version: " + checker.cpp_versions[self.single_host_id])
        else:
            log.info("Remote C++ compiler versions:")
            for host_id in checker.cpp_versions: log.info(" - " + host_id + ": " + checker.cpp_versions[host_id])

        # MPI compiler
        if introspection.has_mpi_compiler(): log.info("Local MPI compiler version: " + introspection.mpi_compiler_version())
        if self.nremotes == 1 and self.single_host_id in checker.mpi_versions: log.info("Remote MPI compiler version: " + checker.mpi_versions[self.single_host_id])
        else:
            log.info("Remote MPI compiler versions:")
            for host_id in checker.mpi_versions: log.info(" - " + host_id + ": " + checker.mpi_versions[host_id])

        # Qt version
        log.info("Local Qt version: " + introspection.qmake_version())
        if self.nremotes == 1 and self.single_host_id in checker.qt_versions: log.info("Remote Qt version: " + checker.qt_versions[self.single_host_id])
        else:
            log.info("Remote Qt versions:")
            for host_id in checker.qt_versions: log.info(" - " + host_id + ": " + checker.qt_versions[host_id])

        # Python
        log.info("Local Python version: " + introspection.python_version_long())
        if self.nremotes == 1 and self.single_host_id in checker.python_versions: log.info("Remote Python version: " + checker.python_versions[self.single_host_id])
        else:
            log.info("Remote Python versions:")
            for host_id in checker.python_versions: log.info(" - " + host_id + ": " + checker.python_versions[host_id])

        # SKIRT
        log.info("Local SKIRT version: " + introspection.skirt_version())
        if self.nremotes == 1 and self.single_host_id in checker.skirt_versions: log.info("Remote SKIRT version: " + checker.skirt_versions[self.single_host_id])
        else:
            log.info("Remote SKIRT versions:")
            for host_id in checker.skirt_versions: log.info(" - " + host_id + ": " + checker.skirt_versions[host_id])

        # PTS
        log.info("Local PTS version: " + introspection.pts_version())
        if self.nremotes == 1 and self.single_host_id in checker.pts_versions: log.info("Remote PTS version: " + checker.pts_versions[self.single_host_id])
        else:
            log.info("Remote PTS versions:")
            for host_id in checker.pts_versions: log.info(" - " + host_id + ": " + checker.pts_versions[host_id])

# -----------------------------------------------------------------
