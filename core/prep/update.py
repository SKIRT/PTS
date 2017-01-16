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
from .installation import get_installation_commands, find_qmake, build_skirt_on_remote
from ..tools import filesystem as fs
from ..tools import git

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

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Update
        self.update()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Updater, self).setup()

        # Check if a remote instance is passed
        if "remote" in kwargs:

            # Get the remote
            self.remote = kwargs.pop("remote")
            if not self.remote.connected: raise ValueError("Remote has not been setup")

        # Setup the remote execution environment if necessary
        elif self.config.host_id is not None:

            # Create and setup the remote execution environment
            self.remote = Remote()
            self.remote.setup(self.config.host_id)

        # Local
        else: log.info("No remote host is specified, will be updating locally ...")

        # Fix configuration files
        self.remote.fix_configuration_files()

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        if self.remote is None: return self.config.host_id
        else: return self.remote.host_id

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

        # The git version to which SKIRT is updated
        self.git_version = None

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

        # Flag that says whether a rebuild is necessary
        self.rebuild = True

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
        if self.rebuild: self.build_local()

    # -----------------------------------------------------------------

    def check_compilers_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers locally ...")

        # Set C++ compiler path
        if self.remote.has_cpp_compiler: self.compiler_path = self.remote.cpp_compiler_path

        # Set MPI compiler path
        if self.remote.has_mpi_compiler: self.mpi_compiler_path = self.remote.mpi_compiler_path

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        # Find qmake path
        self.qmake_path = find_qmake()

        # Check
        if self.qmake_path is None: raise RuntimeError("qmake was not found")

    # -----------------------------------------------------------------

    def pull_local(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Getting latest version ...")

        # Call the appropriate git command at the SKIRT repository directory
        subprocess.call(["git", "pull", "origin", "master"], cwd=introspection.skirt_repo_dir)

    # -----------------------------------------------------------------

    def build_local(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Compiling latest version ...")

        # Recompile SKIRT
        subprocess.call(["sh", "makeSKIRT.sh"], cwd=introspection.skirt_repo_dir)

    # -----------------------------------------------------------------

    def update_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating SKIRT remotely on host '" + self.host_id + "' ...")

        # Check the compilers (C++ and MPI)
        self.check_compilers_remote()

        # Check Qt installation, find qmake
        self.check_qt_remote()

        # Check the presence of git:
        # no, loading the latest git version interferes with the intel compiler version on HPC UGent
        # self.check_git_remote()

        # Pull
        self.pull_remote()

        # Build
        if self.rebuild: self.build_remote()

    # -----------------------------------------------------------------

    def check_compilers_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers on the remote host ...")

        # Get the compiler paths
        self.compiler_path = self.remote.find_and_load_cpp_compiler()
        self.mpi_compiler_path = self.remote.find_and_load_mpi_compiler()

    # -----------------------------------------------------------------

    def check_git_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of git ...")

        # Find and load git
        path, version = self.remote.find_and_load_git()

        # Debugging
        log.debug("The path of the git installation is '" + path)
        log.debug("The version of git is '" + version + "'")

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

        # Check
        if self.qmake_path is None: raise RuntimeError("qmake was not found")

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

        #if self.remote.host.name == "login.hpc.ugent.be":

            # Get path to the origin.txt file
            #origin_path = fs.join(self.remote.skirt_repo_path, "origin.txt")

            # Installed with PTS
            #if self.remote.is_file(origin_path):

                # Get url
                #url, git_version, git_hash = first_three_items_in_iterator(self.remote.read_lines(origin_path))

                # Get the latest git hash for the URL of the repository
                #latest_git_hash = get_hash_remote_repository(url)

                # Debugging
                #log.debug("Latest git hash for repository: " + latest_git_hash)
                #log.debug("Git hash of version installed: " + git_hash)

                # Compare the hashes
                #if latest_git_hash == git_hash:

                    # Already up to date
                    #log.success("Already up to date")
                    #self.rebuild = False
                    #self.git_version = git_version
                    #return

                #else: # not up to date

                    # Remove the previous SKIRT/git directory
                    #self.remote.remove_directory(self.remote.skirt_repo_path)

                    # Get the new code
                    #self.git_version = get_skirt_hpc(self.remote, url, self.remote.skirt_root_path, self.remote.skirt_repo_path)

            # Not installed with PTS
            #else:

                # Get the url of the "origin"
                #url = get_url_repository(self.remote, self.remote.skirt_repo_path)

                # Get latest hash
                #latest_git_hash = get_hash_remote_repository(url)

                # Get hash of current state
                #git_hash = get_git_hash(self.remote, self.remote.skirt_repo_path)

                # Debugging
                #log.debug("Latest git hash for repository: " + latest_git_hash)
                #log.debug("Git hash of version installed: " + git_hash)

                # Compare git hashes
                #if latest_git_hash == git_hash:
                #    log.success("Already up to date")  # up to date
                #    self.rebuild = False
                #    self.git_version = get_short_git_version(self.remote, self.remote.skirt_repo_path)
                #    return

                # The repository can only be installed through SSH, thus no password is required
                #command = "git pull origin master"

                # Execute the command
                #self.remote.execute(command, cwd=self.remote.skirt_repo_path, show_output=True)

                # Get the git version
                #self.git_version = get_short_git_version(self.remote, self.remote.skirt_repo_path)

        # Else
        #else:

        # Get the url of the "origin"
        url = git.get_url_repository(self.remote, self.remote.skirt_repo_path)

        # Get latest hash
        latest_git_hash = git.get_hash_remote_repository(url)

        # Get hash of current state
        git_hash = git.get_git_hash(self.remote, self.remote.skirt_repo_path)

        # Debugging
        log.debug("Latest git hash for repository: " + latest_git_hash)
        log.debug("Git hash of version installed: " + git_hash)

        # Compare git hashes
        if latest_git_hash == git_hash:

            log.success("Already up to date")  # up to date
            self.rebuild = False
            self.git_version = git.get_short_git_version(self.remote.skirt_repo_path, self.remote)
            return

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

        # Get the url of the "origin"
        url = git.get_url_repository(self.remote, self.remote.skirt_repo_path)

        # Decompose
        host, user_or_organization, repo_name, username, password = git.decompose_https(url)

        # Find the account file for the repository host (e.g. github.ugent.be)
        if username is None and password is None and introspection.has_account(host):

            username, password = introspection.get_account(host)

            # Set the command lines
            lines = []
            lines.append(command)
            lines.append(("':", username))
            lines.append(("':", password))

            # Clone the repository
            self.remote.execute_lines(*lines, show_output=log.is_debug())

        # Pull
        else: self.remote.execute(command, show_output=log.is_debug(), cwd=self.remote.skirt_repo_path)

        # Get git version
        self.git_version = git.get_short_git_version(self.remote.skirt_repo_path, self.remote)

        # Success
        log.success("SKIRT was successfully updated on remote host " + self.remote.host_id)

    # -----------------------------------------------------------------

    def build_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Compiling latest version ...")

        # Execute the build
        build_skirt_on_remote(self.remote, self.remote.skirt_repo_path, self.qmake_path, self.git_version)

        # Success
        log.success("SKIRT was successfully built")

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
        self.pull_local()

        # Update the dependencies
        self.update_dependencies_local()

    # -----------------------------------------------------------------

    def pull_local(self):

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
        self.pull_remote()

        # Update the dependencies
        self.update_dependencies_remote()

    # -----------------------------------------------------------------

    def pull_remote(self):

        """
        This function ...
        :return:
        """

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
        self.remote.execute_lines(*lines, show_output=log.is_debug())

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

            if isinstance(command, list): self.remote.execute_lines(*command, show_output=log.is_debug())
            elif isinstance(command, basestring): self.remote.execute(command, show_output=log.is_debug())

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
        self.remote.execute_lines(*lines, show_output=log.is_debug())

# -----------------------------------------------------------------

def first_item_in_iterator(iterator):

    """
    This function ...
    :param iterator:
    :return:
    """

    for item in iterator:
        output = item
        break
    else: raise ValueError("Iterator contains no elements")

    return output

# -----------------------------------------------------------------

def first_two_items_in_iterator(iterator):

    """
    This function ...
    :param iterator:
    :return:
    """

    output = []
    for item in iterator:
        output.append(item)
        if len(output) == 2: break
    else: raise ValueError("Iterator contains less than 2 elements")

    return output[0], output[1]

# -----------------------------------------------------------------

def first_three_items_in_iterator(iterator):

    """
    This function ...
    :param iterator:
    :return:
    """

    output = []
    for item in iterator:
        output.append(item)
        if len(output) == 3: break
    else: raise ValueError("Iterator contains less than 3 elements")

    return output[0], output[1], output[2]

# -----------------------------------------------------------------
