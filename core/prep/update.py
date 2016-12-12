#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.update Contains the SkirtUpdater class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import subprocess
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..remote.remote import Remote
from ..tools.logging import log
from ..tools import introspection
from .installation import get_installation_commands, get_skirt_hpc

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

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SKIRTUpdater, self).__init__(config)

        # The paths to the C++ compiler and MPI compiler
        self.compiler_path = None
        self.mpi_compiler_path = None

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

    # -----------------------------------------------------------------

    def update_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating SKIRT locally ...")

        # Check the compilers (C++ and MPI)
        self.check_compilers_local()

        # Check Qt installation, find qmake
        self.check_qt_local()

        # Pull
        self.pull_local()

        # Build
        self.build_local()

    # -----------------------------------------------------------------

    def check_compilers_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def pull_local(self):

        """
        This function ...
        :return:
        """

        # Get the SKIRT repo directory
        skirt_repo_path = introspection.skirt_repo_dir

        # Debugging
        log.debug("Getting latest version ...")

        # Call the appropriate git command at the SKIRT repository directory
        subprocess.call(["git", "pull", "origin", "master"], cwd=skirt_repo_path)

    # -----------------------------------------------------------------

    def build_local(self):

        """
        This function ...
        :return:
        """

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

        # Check the compilers (C++ and MPI)
        self.check_compilers_remote()

        # Check Qt installation, find qmake
        self.check_qt_remote()

        # Pull
        self.pull_remote()

        # Build
        self.build_remote()

    # -----------------------------------------------------------------

    def check_compilers_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers ...")

        # Get the compiler paths
        self.compiler_path = self.remote.find_and_load_cpp_compiler()
        self.mpi_compiler_path = self.remote.find_and_load_mpi_compiler()

    # -----------------------------------------------------------------

    def check_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking for Qt installation on remote ...")

        # Load Qt module, find the qmake path
        self.qmake_path = self.remote.find_and_load_qmake()

    # -----------------------------------------------------------------

    def pull_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Getting latest version ...")

        # Do HPC UGent in a different way because it seems only SSH is permitted and not HTTPS (but we don't want SSH
        # because of the private/public key thingy, so use a trick
        if self.remote.host.name == "login.hpc.ugent.be":

            # TODO: implement this

            pass

            # First remove the previous SKIRT/git directory
            #self.remote.remove_directory(self.remote.skirt_repo_path)

            # Get the new code
            #get_skirt_hpc(self.remote, url, self.remote.skirt_root_path, self.remote.skirt_repo_path)

        # Else
        else:

            # Change working directory to repository directory
            skirt_git_path = self.remote.skirt_repo_path
            self.remote.change_cwd(skirt_git_path)

            # FOR SSH:
            # Git pull
            #self.remote.ssh.sendline("git pull origin master")
            #self.remote.ssh.expect(":")
            # Expect password
            #if "Enter passphrase for key" in self.remote.ssh.before:
            #    self.remote.execute(self.config.pubkey_password, show_output=True)
            #else: self.remote.prompt()

            # Set the clone command
            command = "git pull origin master"

            # Get the url of the repo from which cloned
            args = ["git", "remote", "show", "origin"]
            output = self.remote.execute(args)
            url = None
            for line in output:
                if "Fetch URL" in line: url = line.split(": ")[1]
            # url = "https://" + host + "/" + user_or_organization + "/" + repo_name + ".git"
            host = url.split("https://")[1].split("/")[0]

            # Find the account file for the repository host (e.g. github.ugent.be)
            username, password = introspection.get_account(host)

            # Set the command lines
            lines = []
            lines.append(command)
            lines.append(("':", username))
            lines.append(("':", password))

            # Clone the repository
            self.remote.execute_lines(*lines, show_output=True)

        # Success
        log.success("SKIRT was successfully updated")

    # -----------------------------------------------------------------

    def build_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Compiling latest version ...")

        # Build SKIRT
        #self.remote.execute("./makeSKIRT.sh", show_output=True, output=False)

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

        # Update PTS
        self.update_pts_local()

        # Update the dependencies
        self.update_dependencies_local()

    # -----------------------------------------------------------------

    def update_pts_local(self):

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

    def update_dependencies_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def update_remote(self):

        """
        This function ...
        :return:
        """

        # Update PTS
        self.update_pts_remote()

        # Update the dependencies
        self.update_dependencies_remote()

    # -----------------------------------------------------------------

    def update_pts_remote(self):

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

        # Set the command
        command = "git pull origin master"

        # Get username and password
        username, password = introspection.get_account("github.ugent.be")

        # Set the command lines
        lines = []
        lines.append(command)
        lines.append(("':", username))
        lines.append(("':", password))

        # Clone the repository
        self.remote.execute_lines(*lines, show_output=True)

        #if "Enter passphrase for key" in self.remote.ssh.before:
        #    self.remote.execute(self.config.pubkey_password, show_output=True)

        #else: self.remote.prompt()

    # -----------------------------------------------------------------

    def update_dependencies_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the PTS dependencies ...")

        # NOTE: this function has completely the same implementation as the function 'get_dependencies_remote' in the PTSInstaller class ! (TODO: accomodate this)

        # Get available conda packages
        output = self.remote.execute("conda search")
        available_packages = []
        for line in output:
            if not line.split(" ")[0]: continue
            available_packages.append(line.split(" ")[0])

        # Get already installed packages
        already_installed = []
        for line in self.remote.execute("conda list"):
            if line.startswith("#"): continue
            already_installed.append(line.split(" ")[0])

        # Use the introspection module on the remote end to get the dependencies and installed python packages
        python = self.remote.start_python_session()
        python.import_package("introspection", from_name="pts.core.tools")
        dependencies = python.get_simple_property("introspection", "get_all_dependencies().keys()")
        packages = python.get_simple_property("introspection", "installed_python_packages()")
        #self.remote.end_python_session()
        # Don't end the python session just yet

        # Get installation commands
        installation_commands, installed, not_installed = get_installation_commands(dependencies, packages, already_installed, available_packages, self.remote)

        # End the python session
        del python

        # Install
        for module in installation_commands:

            # Debugging
            log.debug("Installing '" + module + "' ...")

            command = installation_commands[module]

            if isinstance(command, list): self.remote.execute_lines(*command, show_output=True)
            elif isinstance(command, basestring): self.remote.execute(command, show_output=True)

        # Show installed packages
        log.info("Packages that were installed:")
        for module in installed: log.info(" - " + module)

        # Show not installed packages
        log.info("Packages that could not be installed:")
        for module in not_installed: log.info(" - " + module)

        # Show already present packages
        log.info("Packages that were already present:")
        for module in already_installed: log.info(" - " + module)

        # Update all packages
        update_command = "conda update --all"
        lines = []
        lines.append(update_command)
        lines.append(("Proceed ([y]/n)?", "y", True))
        self.remote.execute_lines(*lines, show_output=True)

# -----------------------------------------------------------------
