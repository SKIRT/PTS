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
from .installation import find_qmake, build_skirt_on_remote, build_skirt_local, get_pts_dependencies_remote
from .installation import has_valid_conda_environment_local, has_valid_conda_environment_remote, install_conda_remote
from .installation import create_conda_environment_local, create_conda_environment_remote
from ..tools import git
from ..tools import terminal
from ..tools import filesystem as fs
from ..remote.modules import Modules

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

        # 3. Tes
        self.test()

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

            # Fix configuration files
            self.remote.fix_configuration_files()

        # Local
        else: log.info("No remote host is specified, will be updating locally ...")

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

    def test(self):

        """
        This function ...
        :return:
        """

        # Test locally or remotely
        if self.remote is None: self.test_local()
        else: self.test_remote()

    # -----------------------------------------------------------------

    @abstractmethod
    def test_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def test_remote(self):

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

        # 1. Check the compilers (C++ and MPI)
        self.check_compilers_local()

        # 2. Check Qt installation, find qmake
        self.check_qt_local()

        # 3. Pull
        self.pull_local()

        # 4. Build
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

        # Debugging
        if self.compiler_path is not None: log.debug("The C++ compiler path is '" + self.compiler_path)
        else: log.debug("A C++ compiler is not found")
        if self.mpi_compiler_path is not None: log.debug("The MPI compiler path is '" + self.mpi_compiler_path)
        else: log.debug("An MPI compiler is not found")

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of Qt locally ...")

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

        # Inform the user
        log.info("Building SKIRT ...")

        # Execute the build
        build_skirt_local(introspection.skirt_repo_dir, self.qmake_path, self.git_version)

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

        # Debugging
        log.debug("The C++ compiler path is '" + self.compiler_path)
        if self.mpi_compiler_path is not None: log.debug("The MPI compiler path is '" + self.mpi_compiler_path)
        else: log.debug("An MPI compiler is not found")

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

        # Set the clone command
        command = "git pull origin master"

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
            self.remote.execute_lines(*lines, show_output=log.is_debug(), cwd=self.remote.skirt_repo_path)

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

    def test_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation ...")

        output = terminal.execute("skirt -h")
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were succesfully updated")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation on the remote host ...")

        output = self.remote.execute("skirt -h")
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were successfully updated on the remote host")

# -----------------------------------------------------------------

class PTSUpdater(Updater):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(PTSUpdater, self).__init__(config)

        # The git version to which PTS is updated
        self.git_version = None

        # Paths
        self.conda_installation_path = None
        self.conda_main_executable_path = None
        self.conda_executable_path = None
        self.conda_pip_path = None
        self.conda_activate_path = None
        self.conda_python_path = None
        self.conda_easy_install_path = None

        # The conda environment
        self.conda_environment = None
        self.python_version = None

    # -----------------------------------------------------------------

    @property
    def has_conda(self):

        """
        This function ...
        :return:
        """

        return self.conda_main_executable_path is not None

    # -----------------------------------------------------------------

    def update_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating PTS and its dependencies locally ...")

        # 1. Update PTS
        self.pull_local()

        # 2. Check conda
        self.check_conda_local()

        # 2. Get conda locally
        if not self.has_conda: self.get_conda_local()

        # 3. Create python environment
        if not has_valid_conda_environment_local(self.conda_environment, self.conda_main_executable_path, self.conda_installation_path): self.create_environment_local()
        else: self.set_environment_paths_local()

        # 4. Update conda
        if self.config.dependencies: self.update_conda_local()

        # 5. Update the dependencies
        if self.config.dependencies: self.update_dependencies_local()

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

    def check_conda_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking conda locally ...")

    # -----------------------------------------------------------------

    def get_conda_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_environment_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a fresh python environment ...")

        # Create the environment
        # conda_executable_path, conda_pip_path, conda_activate_path, conda_python_path, conda_easy_install_path
        self.conda_executable_path, self.conda_pip_path, self.conda_activate_path, self.conda_python_path, self.conda_easy_install_path = \
            create_conda_environment_local(self.conda_environment, self.conda_installation_path, introspection.pts_root_dir,
                                           self.python_version, self.conda_main_executable_path)

    # -----------------------------------------------------------------

    def set_environment_paths_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def update_conda_local(self):

        """
        This function ...
        :return:
        """

        pass

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

        # Inform the user
        log.info("Updating PTS on the remote host ...")

        # 1. Pull
        self.pull_remote()

        # 2. Check conda, get paths
        self.check_conda_remote()

        # 2. Get conda
        if not self.has_conda: self.get_conda_remote()

        # 2. Create a remote python environment
        if not has_valid_conda_environment_remote(self.remote, self.conda_environment, self.conda_main_executable_path, self.conda_installation_path): self.create_environment_remote()
        else: self.set_environment_paths_remote()

        # 3. Update conda
        if self.config.dependencies: self.update_conda_remote()

        # 4. Update the dependencies
        if self.config.dependencies: self.update_dependencies_remote()

    # -----------------------------------------------------------------

    def pull_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Getting latest version ...")

        # Get the url of the "origin"
        url = git.get_url_repository(self.remote, self.remote.pts_package_path)

        # Get latest hash
        latest_git_hash = git.get_hash_remote_repository(url)

        # Get hash of current state
        git_hash = git.get_git_hash(self.remote, self.remote.pts_package_path)

        # Debugging
        log.debug("Latest git hash for repository: " + latest_git_hash)
        log.debug("Git hash of version installed: " + git_hash)

        # Compare git hashes
        if latest_git_hash == git_hash:

            log.success("Already up to date")  # up to date
            self.git_version = git.get_short_git_version(self.remote.pts_package_path, self.remote)
            return

        # Set the command
        command = "git pull origin master"

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
            self.remote.execute_lines(*lines, show_output=log.is_debug(), cwd=self.remote.pts_package_path)

        # Pull
        else: self.remote.execute(command, show_output=log.is_debug(), cwd=self.remote.pts_package_path)

        # Get git version
        self.git_version = git.get_short_git_version(self.remote.pts_package_path, self.remote)

        # Success
        log.success("PTS was successfully updated on remote host " + self.remote.host_id)

    # -----------------------------------------------------------------

    def check_conda_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of a conda installation on the remote host ...")

        # Find conda path
        if self.remote.is_executable("conda"): conda_executable_path = self.remote.find_executable("conda")
        else: conda_executable_path = None

        # Find conda installation in the home directory
        if conda_executable_path is None:

            conda_path = fs.join(self.remote.home_directory, "miniconda", "bin", "conda")
            if self.remote.is_file(conda_path): conda_executable_path = conda_path

        # If conda is present
        if conda_executable_path is not None:

            self.conda_installation_path = fs.directory_of(fs.directory_of(conda_executable_path))
            self.conda_main_executable_path = conda_executable_path

            # Debugging
            log.debug("Determining the conda environment name for PTS ...")

            # Get environment name
            pts_alias = self.remote.resolve_alias("pts")
            pts_python_path = pts_alias.split()[0]
            if pts_python_path == "python": raise Exception("Cannot determine the conda environment used for pts")
            else:

                env_path = fs.directory_of(fs.directory_of(pts_python_path))
                environment_name = fs.name(env_path)
                self.conda_environment = environment_name

        # Conda not found
        else:

            # Warning
            log.warning("Conda is not found: installing")

            # Set the conda environment name
            self.conda_environment = "python_pts"
            self.python_version = "2.7"

    # -----------------------------------------------------------------

    def get_conda_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting a Conda python distribution on the remote host ...")

        # Install conda
        self.conda_installation_path, self.conda_main_executable_path = install_conda_remote(self.remote)

    # -----------------------------------------------------------------

    def create_environment_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a fresh python environment ...")

        # Create the environment
        # conda_executable_path, conda_pip_path, conda_activate_path, conda_python_path, conda_easy_install_path
        self.conda_executable_path, self.conda_pip_path, self.conda_activate_path, self.conda_python_path, self.conda_easy_install_path = \
            create_conda_environment_remote(self.remote, self.conda_environment, self.conda_installation_path,
                                        self.remote.pts_root_path, self.python_version, self.conda_main_executable_path)

    # -----------------------------------------------------------------

    def set_environment_paths_remote(self):

        """
        This function ...
        :return:
        """

        # Set paths
        environment_bin_path = fs.join(self.conda_installation_path, "envs", self.conda_environment, "bin")
        if not self.remote.is_directory(environment_bin_path): raise RuntimeError("The environment directory is not present")
        self.conda_executable_path = fs.join(environment_bin_path, "conda")
        self.conda_pip_path = fs.join(environment_bin_path, "pip")
        self.conda_activate_path = fs.join(environment_bin_path, "activate")
        self.conda_python_path = fs.join(environment_bin_path, "python")
        self.conda_easy_install_path = fs.join(environment_bin_path, "easy_install")

        # Check if paths exist
        assert self.remote.is_file(self.conda_executable_path)
        assert self.remote.is_file(self.conda_pip_path)
        assert self.remote.is_file(self.conda_activate_path)
        assert self.remote.is_file(self.conda_python_path)
        assert self.remote.is_file(self.conda_easy_install_path)

    # -----------------------------------------------------------------

    def update_conda_remote(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def update_dependencies_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the PTS dependencies ...")

        # Install PTS dependencies
        # conda_path="conda", pip_path="pip", python_path="python",
        # easy_install_path="easy_install", conda_environment=None
        installed, not_installed, already_installed = get_pts_dependencies_remote(self.remote, self.conda_executable_path,
                                                                                  self.conda_pip_path, self.conda_python_path,
                                                                                  self.conda_easy_install_path, self.conda_environment, conda_activate_path=self.conda_activate_path)

        # Get dependency version restrictions
        versions = introspection.get_constricted_versions()

        # Get names of packages for import names
        real_names = introspection.get_package_names()

        # Get the versions of the packages currently installed
        conda_versions = self.remote.installed_conda_packages(self.conda_executable_path, self.conda_environment)

        # Loop over the already installed packages and update them if permitted
        for module_name in not_installed:

            # Check if the module name may be different from the import name
            if module_name in real_names.values():
                import_name = real_names.keys()[real_names.values().index(module_name)]
            else: import_name = module_name

            # Check if there is a version restriction
            if import_name in versions:

                version = versions[import_name]

                # Debugging
                log.debug("Package '" + module_name + "' has a version restriction: " + version)

                # Debugging
                log.debug("Checking version ...")

                # Get installed version
                installed_version = conda_versions[module_name]
                if version != installed_version: log.error("The version of the package '" + module_name + "' is " + installed_version + " but it should be " + version)

            # No version restriction: update
            else:

                # Debugging
                log.debug("Updating package '" + module_name + "' to the latest version ...")

                # Set command
                command = self.conda_executable_path + " update " + module_name + " -n " + self.conda_environment

                # Debugging
                log.debug("Update command: " + command)

                # Launch the command
                self.remote.ssh.sendline(command)

                # Expect the prompt or question
                while self.remote.ssh.expect([self.remote.ssh.PROMPT, "Proceed ([y]/n)?"], timeout=None) == 1: self.remote.ssh.sendline("y")

        # Success
        log.success("Succesfully installed and updated the dependencies on the remote host")

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation ...")

        # Test
        output = terminal.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were successfully updated ...")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation ...")

        # Test
        output = self.remote.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were successfully updated ...")

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
