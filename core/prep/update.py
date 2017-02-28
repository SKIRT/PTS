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
from ..tools import conda
from .installation import montage_url, imfit_macos_binary_url, imfit_linux_binary_url

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

            # Fix configuration files
            self.remote.fix_configuration_files()

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

        # Modules
        self.modules = None

        # The paths to the C++ compiler and MPI compiler
        self.compiler_path = None
        self.mpi_compiler_path = None

        # The modules
        self.compiler_module = None
        self.mpi_compiler_module = None

        # The git version to which SKIRT is updated
        self.git_version = None

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

        # Qmake modules
        self.qmake_module = None

        # Git path and module
        self.git_path = None
        self.git_module = None

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

        # Check modules
        self.check_modules_remote()

        # Pull
        self.pull_remote()

        # Build
        if self.rebuild: self.build_remote()

    # -----------------------------------------------------------------

    def check_modules_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the modules on the remote host ...")

        # Create the modules
        self.modules = Modules(self.remote)

        # Check compilers
        self.check_compilers_remote()

        # Check Qt
        self.check_qt_remote()

        # Check git
        self.check_git_remote()

    # -----------------------------------------------------------------

    def check_compilers_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers ...")

        # C++ compiler
        cpp_path = self.modules.paths["cpp"]
        cpp_version = self.modules.versions["cpp"]
        cpp_module = self.modules.names["cpp"]

        # Check if module has to be loaded
        if cpp_module is not None:

            # Inform the user
            log.info("The module '" + cpp_module + "' has to be loaded for qmake")

            # Set the module name and the path
            self.compiler_module = cpp_module
            self.compiler_path = cpp_path

        # No module required
        else:

            # Inform the user
            log.info("The C++ compiler is present at '" + cpp_path + "'")

            # Set the path
            self.compiler_path = cpp_path

        # MPI compiler
        mpi_path = self.modules.paths["mpi"]
        mpi_version = self.modules.versions["mpi"]
        mpi_module = self.modules.names["mpi"]

        # Check if MPI compiler is present
        if mpi_path is not None:

            # Check if module has to be loaded
            if mpi_module is not None:

                # Inform the user
                log.info("The module '" + mpi_module + "' has to be loaded for MPI")

                # Set the module name and the path
                self.mpi_compiler_module = mpi_module
                self.mpi_compiler_path = mpi_path

            # No module required
            else:

                # Inform the user
                log.info("The MPI compiler is present at '" + mpi_path + "'")

                # Set the path
                self.mpi_compiler_path = mpi_path

        # No MPI compiler found
        else: log.info("An MPI compiler is not present")

    # -----------------------------------------------------------------

    def check_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking for Qt installation on remote ...")

        # Qt is present
        if self.modules.paths["qmake"] is not None:

            # Get the version
            version = self.modules.versions["qmake"]

            # Check if module has to be loaded
            if self.modules.names["qmake"] is not None:

                # Inform the user
                log.info("The module '" + self.modules.names["qmake"] + "' has to be loaded for qmake ...")

                # Set the module name and the qmake path
                self.qmake_module = self.modules.names["qmake"]
                self.qmake_path = self.modules.paths["qmake"]

            else:

                # Inform the user
                log.info("Qmake is present at '" + self.modules.paths["qmake"] + "'")

                # Set qmake path
                self.qmake_path = self.modules.paths["qmake"]

        # Qt is not present
        else: log.info("Qt is not present, will be installing ...")

    # -----------------------------------------------------------------

    def check_git_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of git on the remote host ...")

        # Get info
        path = self.modules.paths["git"]
        version = self.modules.versions["git"]
        module = self.modules.names["git"]

        if path is None: raise RuntimeError("Git could not be found")

        # Set path and module name
        self.git_path = path
        if module is not None: self.git_module = module

    # -----------------------------------------------------------------

    def pull_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting latest version ...")

        # Load the git module
        if self.git_module is not None: self.remote.load_module(self.git_module, show_output=log.is_debug())

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

            log.success("SKIRT is already up to date on remote host '" + self.remote.host_id + "'")  # up to date
            self.rebuild = False
            self.git_version = git.get_short_git_version(self.remote.skirt_repo_path, self.remote)
            self.remote.unload_all_modules()
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

        # Unload all modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def build_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Compiling latest version ...")

        # Load modules
        if self.qmake_module is not None: self.remote.load_module(self.qmake_module, show_output=log.is_debug())
        else:
            if self.compiler_module is not None: self.remote.load_module(self.compiler_module, show_output=log.is_debug())
            if self.mpi_compiler_module is not None and self.mpi_compiler_module != self.compiler_module: self.remote.load_module(self.mpi_compiler_module, show_output=log.is_debug())

        # Execute the build
        build_skirt_on_remote(self.remote, self.remote.skirt_repo_path, self.qmake_path, self.git_version)

        # Success
        log.success("SKIRT was successfully built")

        # Unload modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation locally ...")

        output = terminal.execute("skirt -h")
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were succesfully updated locally")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation on remote host '" + self.remote.host_id + "'...")

        # Load modules
        if self.qmake_module is not None: self.remote.load_module(self.qmake_module, show_output=log.is_debug())
        else:
            if self.compiler_module is not None: self.remote.load_module(self.compiler_module, show_output=log.is_debug())
            if self.mpi_compiler_module is not None and self.mpi_compiler_module != self.compiler_module: self.remote.load_module(self.mpi_compiler_module, show_output=log.is_debug())

        output = self.remote.execute("skirt -h")
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were successfully updated on remote host '" + self.remote.host_id + "'")

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

        # Path to the Montage installation directory
        self.montage_path = None

        # Path to the Imfit installation directory
        self.imfit_path = None

    # -----------------------------------------------------------------

    @property
    def has_conda(self):

        """
        This function ...
        :return:
        """

        return self.conda_main_executable_path is not None

    # -----------------------------------------------------------------

    @property
    def has_montage(self):

        """
        This function ...
        :return:
        """

        return self.montage_path is not None

    # -----------------------------------------------------------------

    @property
    def has_imfit(self):

        """
        This function ...
        :return:
        """

        return self.imfit_path is not None

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
        if self.config.conda: self.update_conda_local()

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

        # Find conda locally
        conda_executable_path = conda.find_conda()

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

        # Check imfit
        self.check_imfit_remote()

        # Get imfit
        if not self.has_imfit: self.get_imfit_remote()

        # Check Montage
        self.check_montage()

        # Get Montage
        if not self.has_montage: self.get_montage()

        # 2. Check conda, get paths
        self.check_conda_remote()

        # 2. Get conda
        if not self.has_conda: self.get_conda_remote()

        # 2. Create a remote python environment
        if not has_valid_conda_environment_remote(self.remote, self.conda_environment, self.conda_main_executable_path, self.conda_installation_path): self.create_environment_remote()
        else: self.set_environment_paths_remote()

        # 3. Update conda
        if self.config.conda: self.update_conda_remote()

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

            log.success("PTS is already up to date on remote host '" + self.remote.host_id + "'")  # up to date
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
        #log.success("PTS was successfully updated on remote host " + self.remote.host_id)

    # -----------------------------------------------------------------

    def check_imfit_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of Imfit ...")

        # Look for Imfit
        executable_path = self.remote.find_executable("imfit")
        if executable_path is not None: self.imfit_path = fs.directory_of(executable_path)

        # Debugging
        if executable_path is not None: log.debug("Imfit installation found in '" + self.imfit_path + "'")

    # -----------------------------------------------------------------

    def get_imfit_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting Montage on the remote host ...")

        # Download the tar.gz file
        if self.remote.is_macos: url = imfit_macos_binary_url
        elif self.remote.is_linux: url = imfit_linux_binary_url
        else: raise NotImplementedError("Platform not supported")
        filepath = self.remote.download_from_url_to(url, self.remote.home_directory, show_output=log.is_debug())

        # Decompress into "~/Imfit"
        imfit_installation_path = fs.join(self.remote.home_directory, "Imfit")
        self.remote.decompress_directory_to(filepath, imfit_installation_path, show_output=log.is_debug(), remove=True)

        # Set the imfit path
        self.imfit_path = imfit_installation_path

        # Add to path
        self.remote.add_to_path_variable(self.imfit_path, in_shell=True)

    # -----------------------------------------------------------------

    def check_montage(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of Montage ...")

        # Look for Montage
        executable_path = self.remote.find_executable("mArchiveGet")
        if executable_path is not None: self.montage_path = fs.directory_of(executable_path)

        # Debugging
        if executable_path is not None: log.debug("Montage installation found in '" + self.montage_path + "'")

    # -----------------------------------------------------------------

    def get_montage(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting Montage on the remote host ...")

        # Download the tar.gz file
        filepath = self.remote.download_from_url_to(montage_url, self.remote.home_directory, show_output=log.is_debug())

        # Extract the file
        montage_path = self.remote.decompress_directory_in_place(filepath, show_output=log.is_debug(), remove=True)

        # Determine the path to the montage source code directory
        #montage_montage_path = fs.join(montage_path, "Montage")

        #self.remote.execute("./configure", cwd=montage_montage_path, show_output=log.is_debug())

        # Compile
        self.remote.execute("make", cwd=montage_path, show_output=log.is_debug())

        montage_bin_path = fs.join(montage_path, "bin")
        #montage_exec_path = fs.join(montage_bin_path, "montage")
        self.montage_path = montage_bin_path

        # Add to path
        self.remote.add_to_path_variable(montage_bin_path, in_shell=True)

    # -----------------------------------------------------------------

    def check_conda_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of a conda installation on the remote host ...")

        # Find conda
        self.conda_installation_path, self.conda_main_executable_path = self.remote.find_conda()

        # If conda is present
        if self.conda_main_executable_path is not None:

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
        installed, not_installed, already_installed = get_pts_dependencies_remote(self.remote, self.remote.pts_package_path, self.conda_executable_path,
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
        log.info("Testing the PTS installation locally ...")

        # Test
        output = terminal.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were successfully updated locally ...")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation on remote host '" + self.remote.host_id + "'...")

        # Test
        output = self.remote.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were successfully updated on remote host '" + self.remote.host_id + "' ...")

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
