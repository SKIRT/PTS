#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.installation Contains the SKIRTInstaller and PTSInstaller classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys
import traceback
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.log import log
from ..tools import network, archive
from ..tools import terminal
from ..tools import git
from ..tools import conda
from ..tools import time
from ..remote.modules import Modules
from ..units.parsing import parse_unit as u
from ..tools import types

# -----------------------------------------------------------------

# For SSH key:
# eval $(ssh-agent)
# ssh-add

# -----------------------------------------------------------------

skirt_directories = ["git", "run", "doc", "release", "debug"]
pts_directories = ["pts", "run", "doc", "temp", "remotes", "user", "ext"]

# -----------------------------------------------f------------------

# Qt URL
qt_55_url = "https://download.qt.io/archive/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"
qt_56_url = "https://download.qt.io/archive/qt/5.6/5.6.3/single/qt-everywhere-opensource-src-5.6.3.tar.xz"
qt_57_url = "https://download.qt.io/archive/qt/5.7/5.7.1/single/qt-everywhere-opensource-src-5.7.1.tar.gz"
qt_58_url = "https://download.qt.io/archive/qt/5.8/5.8.0/single/qt-everywhere-opensource-src-5.8.0.tar.gz"
qt_59_url = "https://download.qt.io/archive/qt/5.9/5.9.2/single/qt-everywhere-opensource-src-5.9.2.tar.xz"
qt_url = qt_59_url

# -----------------------------------------------------------------

# The approximate size for a full installation of Qt
full_qt_installation_size = 1.5 * u("GB")

# -----------------------------------------------------------------

qt_installation_dirname = "Qt"

# -----------------------------------------------------------------

class Installer(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Installer, self).__init__(*args, **kwargs)

        # The remote execution environment
        self.remote = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # 2. Create the necessary directories
        self.create_directories()

        # 2. Install
        self.install()

        # 3. Test the installation
        self.test()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Installer, self).setup()

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
            from ..remote.remote import Remote
            self.remote = Remote(log_conda=True)
            self.remote.setup(self.config.host_id)

            # Fix configuration files
            self.remote.fix_configuration_files()

        # Local
        else: log.info("No remote host is specified, will be installing locally ...")

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

    def create_directories(self):

        """
        THis function ...
        :return:
        """

        # Install locally or remotely
        if self.remote is None: self.create_directories_local()
        else: self.create_directories_remote()

    # -----------------------------------------------------------------

    def install(self):

        """
        This function ...
        :return:
        """

        # Install locally or remotely
        if self.remote is None: self.install_local()
        else: self.install_remote()

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
    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def create_directories_remote(self):

        """
        This fucntion ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def install_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def install_remote(self):

        """
        This function ...
        :return:
        """

        pass

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

# Determine Qt configure options
qt_configure_options = []
#qt_configure_options.append("-prefix '$HOME/Qt'")
qt_configure_options.append("-prefix '$INSTALLATION_PATH'")
qt_configure_options.append("-opensource")
qt_configure_options.append("-confirm-license")
#qt_configure_options.append("-c++11") # not valid anymore
#qt_configure_options.append("-no-javascript-jit") # not valid anymore
qt_configure_options.append("-no-qml-debug")
qt_configure_options.append("-no-gif")
qt_configure_options.append("-no-libpng")
qt_configure_options.append("-no-libjpeg")
qt_configure_options.append("-no-freetype")
qt_configure_options.append("-no-harfbuzz")
qt_configure_options.append("-no-openssl")
qt_configure_options.append("-no-xinput2")
qt_configure_options.append("-no-xcb-xlib")
qt_configure_options.append("-no-glib")
qt_configure_options.append("-no-gui")
qt_configure_options.append("-no-widgets")
#qt_configure_options.append("-no-nis") # not valid anymore
qt_configure_options.append("-no-cups")
qt_configure_options.append("-no-fontconfig")
qt_configure_options.append("-no-dbus")
qt_configure_options.append("-no-xcb")
qt_configure_options.append("-no-eglfs")
qt_configure_options.append("-no-directfb")
qt_configure_options.append("-no-linuxfb")
qt_configure_options.append("-no-kms")
qt_configure_options.append("-no-opengl")
qt_configure_options.append("-nomake tools")
qt_configure_options.append("-nomake examples")

# -----------------------------------------------------------------

class SKIRTInstaller(Installer):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SKIRTInstaller, self).__init__(*args, **kwargs)

        # The paths to the C++ compiler and MPI compiler
        self.compiler_path = None
        self.mpi_compiler_path = None

        # Modules
        self.compiler_module = None
        self.mpi_compiler_module = None

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

        # The module for qmake
        self.qmake_module = None

        # Path and module for git
        self.git_path = None
        self.git_module = None

        # Path to SKIRT root directory
        self.skirt_root_path = None

        # Path to SKIRT/git
        self.skirt_repo_path = None

        # Path to SKIRT/release
        self.skirt_release_path = None

        # Path to the SKIRT executable
        self.skirt_path = None

        # Git version of SKIRT to be installed
        self.git_version = None

        # The modules
        self.modules = None

    # -----------------------------------------------------------------

    @property
    def qt_url(self):

        """
        This function ...
        :return:
        """

        if self.config.qt_version is None: return qt_url
        elif self.config.qt_version == "5.5": return qt_55_url
        elif self.config.qt_version == "5.6": return qt_56_url
        elif self.config.qt_version == "5.7": return qt_57_url
        elif self.config.qt_version == "5.8": return qt_58_url
        elif self.config.qt_version == "5.9": return qt_59_url
        else: raise ValueError("Invalid version indicator")

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the directory structure ...")

        # Set paths
        self.skirt_root_path = fs.join(fs.home, "SKIRT")
        self.skirt_repo_path = fs.join(self.skirt_root_path, "git")
        self.skirt_release_path = fs.join(self.skirt_root_path, "release")

        # Check if already present
        if fs.is_directory(self.skirt_root_path):
            if self.config.force: fs.remove_directories_in_path(self.skirt_root_path, exact_not_name="run")
            #fs.remove_directory(self.skirt_root_path)
            elif self.config.finish: pass
            else: raise RuntimeError("SKIRT is already installed (or partly present)")

        # Make the root directory
        if not self.config.finish: fs.create_directory(self.skirt_root_path)

        # Create the other directories
        for name in skirt_directories:

            # Determine path
            path = fs.join(self.skirt_root_path, name)

            # Create the directory
            if not fs.is_directory(path): fs.create_directory(path)

    # -----------------------------------------------------------------

    def create_directories_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the directory structure ...")

        # Set root path and pacakge path
        self.skirt_root_path = fs.join(self.remote.home_directory, "SKIRT")
        self.skirt_repo_path = fs.join(self.skirt_root_path, "git")
        self.skirt_release_path = fs.join(self.skirt_root_path, "release")

        # Check if already present
        if self.remote.is_directory(self.skirt_root_path):
            if self.config.force: self.remote.remove_directories_in_path(self.skirt_root_path, exact_not_name="run")
            #self.remote.remove_directory(self.skirt_root_path)
            elif self.config.finish: pass
            else: raise RuntimeError("SKIRT is already installed (or partly present) on the remote host")

        # Make the root directory
        if not self.config.finish: self.remote.create_directory(self.skirt_root_path)

        # Create the other directories
        for name in skirt_directories:

            # Determine path
            path = fs.join(self.skirt_root_path, name)

            # Create the directory
            if not self.remote.is_directory(path): self.remote.create_directory(path)

    # -----------------------------------------------------------------

    @property
    def has_skirt_source_local(self):

        """
        This function ...
        :return:
        """

        has_source = not fs.is_empty(self.skirt_repo_path)
        if not has_source: return False

        # Get git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path)

        # Determine path of SKIRT and FitSKIRT main directories with executable
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        #fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")
        if skirt_main_path in terminal.paths_in_path_variable(): terminal.add_to_path_variable(skirt_main_path, comment=comment, in_shell=True)
        #if fitskirt_main_path in terminal.paths_in_path_variable(): terminal.add_to_path_variable(fitskirt_main_path, comment=comment, in_shell=True)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(skirt_main_path, "skirt")

        # Return
        return True

    # -----------------------------------------------------------------

    def install_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing SKIRT locally ...")

        # Check compilers (C++ and mpi)
        self.check_compilers_local()

        # Check if Qt is installed
        self.check_qt_local()

        # Install Qt
        if not self.has_qt: self.install_qt_local()

        # Get the SKIRT code
        if not self.has_skirt_source_local: self.get_skirt_local()

        # Build SKIRT
        self.build_skirt_local()

    # -----------------------------------------------------------------

    def check_compilers_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers locally ...")

        # Look for C++ compiler
        if introspection.has_cpp_compiler(): self.compiler_path = introspection.cpp_compiler_path()
        else: raise RuntimeError("C++ compilers not present on this system")

        # Look for MPI compiler, not necessary
        if introspection.has_mpi_compiler(): self.mpi_compiler_path = introspection.mpi_compiler_path()

        # Debugging
        log.debug("The C++ compiler path is '" + self.compiler_path)
        if self.mpi_compiler_path is not None: log.debug("The MPI compiler path is '" + self.mpi_compiler_path)
        else: log.debug("An MPI compiler is not found")

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        # Find qmake path
        path = find_qmake()

        if self.config.reinstall_qt:

            # Debugging
            log.debug("Removing files and directories for re-installation of Qt ...")

            source_path = fs.join(fs.home, "qt.tar.gz")
            other_source_path = fs.join(fs.home, "qt.tar.xz")

            decompress_path = fs.join(fs.home, "Qt-install")
            installation_path = fs.join(fs.home, "Qt")

            if fs.is_file(source_path): fs.remove_file(source_path)
            if fs.is_file(other_source_path): fs.remove_file(other_source_path)

            if fs.is_directory(decompress_path): fs.remove_directory(decompress_path)
            if fs.is_directory(installation_path): fs.remove_directory(installation_path)

        elif path is not None: self.qmake_path = path

    # -----------------------------------------------------------------

    def install_qt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing Qt ...")

        # Get the extension of the file
        extension = fs.get_extension(self.qt_url, double=True)

        # Determine the path for the Qt source code
        path = fs.join(fs.home, "qt." + extension)

        # Download tar.gz file
        network.download_file(self.qt_url, path)

        # Decompress tar.gz file
        decompress_path = fs.join(fs.home, "Qt-install")
        fs.create_directory(decompress_path)
        archive.decompress_file(path, decompress_path)

        installation_path = fs.join(fs.home, "Qt")
        configure_options = qt_configure_options[:]
        configure_options[0] = configure_options[0].replace("$INSTALLATION_PATH", installation_path)

        # Determine commands
        configure_command = "./configure " + " ".join(configure_options)
        make_command = "make"
        install_command = "make install"

        # Get the only directory in the Qt-install directory
        qt_everywhere_opensource_path = fs.directories_in_path(decompress_path)[0]

        # Configure
        terminal.execute(configure_command, cwd=qt_everywhere_opensource_path, show_output=log.is_debug)

        # Make
        terminal.execute(make_command, cwd=qt_everywhere_opensource_path, show_output=log.is_debug)

        # Install
        terminal.execute(install_command, cwd=qt_everywhere_opensource_path, show_output=log.is_debug)

        # Remove the tar.gz file and the temporary directory
        fs.remove_file(path)
        fs.remove_directory(decompress_path)

    # -----------------------------------------------------------------

    def get_skirt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the SKIRT source code ...")

        # Get repository link
        if self.config.repository is not None: raise ValueError("Repository cannot be specified for local installation")
        elif self.config.private: url = introspection.private_skirt_https_link
        else: url = introspection.public_skirt_https_link

        # Get host (github.ugent.be)
        host = url.split("//")[1].split("/")[0]

        # Set the clone command
        command = "git clone " + url + " " + self.skirt_repo_path

        # Find the account file for the repository host (e.g. github.ugent.be)
        if introspection.has_account(host):

            username, password = introspection.get_account(host)

            # Set the command lines
            lines = []
            lines.append(command)
            lines.append(("':", username))
            lines.append(("':", password))

            # Clone the repository
            terminal.execute_lines(*lines, cwd=self.skirt_root_path, show_output=log.is_debug)

        else: terminal.execute(command, cwd=self.skirt_root_path, show_output=log.is_debug)

        # Get git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Determine path of SKIRT and FitSKIRT main directories with executable
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        #fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")
        if skirt_main_path in terminal.paths_in_path_variable(): terminal.remove_from_path_variable(skirt_main_path)
        #if fitskirt_main_path in terminal.paths_in_path_variable(): terminal.remove_from_path_variable(fitskirt_main_path)

        # Add to path, also execute in current shell to make SKIRT (and FtiSKIRT) visible
        terminal.add_to_path_variable(skirt_main_path, comment=comment, in_shell=True)
        #terminal.add_to_path_variable(fitskirt_main_path, comment=comment, in_shell=True)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(skirt_main_path, "skirt")

    # -----------------------------------------------------------------

    def build_skirt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building SKIRT ...")

        # Execute the build
        build_skirt_local(self.skirt_repo_path, self.qmake_path, self.git_version)

        # Success
        log.success("SKIRT was successfully built")

    # -----------------------------------------------------------------

    @property
    def has_qt(self):

        """
        This function ...
        :return:
        """

        return self.qmake_path is not None

    # -----------------------------------------------------------------

    @property
    def has_skirt_source_remote(self):

        """
        This function ...
        :return:
        """

        has_source = not self.remote.is_empty(self.skirt_repo_path)
        if not has_source: return False

        # Get the git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path, self.remote)

        # Determine SKIRT and FitSKIRT main paths
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        #fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")

        # Add to environment variables
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        if skirt_main_path not in self.remote.paths_in_path_variable: self.remote.add_to_path_variable(skirt_main_path, comment=comment, in_shell=True)
        #if fitskirt_main_path not in self.remote.paths_in_path_variable: self.remote.add_to_path_variable(fitskirt_main_path, comment=comment, in_shell=True)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(skirt_main_path, "skirt")

        # Return
        return True

    # -----------------------------------------------------------------

    def install_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing SKIRT remotely ...")

        # Check the modules
        self.check_modules_remote()

        # Install Qt if necessary
        if not self.has_qt: self.install_qt_remote()

        # Check the presence of git:
        # no, loading the latest git version interferes with the intel compiler version on HPC UGent
        self.check_git_remote()

        # Get the SKIRT code
        if not self.has_skirt_source_remote: self.get_skirt_remote()

        # Build SKIRT
        self.build_skirt_remote()

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

        if not self.modules.has_cpp: raise RuntimeError("C++ compiler could not be found")

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

        # Check if MPI compiler is present
        if self.modules.has_mpi:

            # MPI compiler
            mpi_path = self.modules.paths["mpi"]
            mpi_version = self.modules.versions["mpi"]
            mpi_module = self.modules.names["mpi"]

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

        # Re-install?
        if self.config.reinstall_qt:

            # Debugging
            log.debug("Removing remote files and directories for re-installation of Qt ...")

            # Get qmake path
            if self.modules.has_qmake: qmake_path = self.modules.paths["qmake"]
            else: qmake_path = None

            temp_path = self.remote.home_directory
            alternative_temp_path = self.remote.scratch_path

            source_path = fs.join(temp_path, "qt.tar.gz")
            other_source_path = fs.join(temp_path, "qt.tar.xz")

            installation_path = fs.join(temp_path, qt_installation_dirname)
            if alternative_temp_path is not None: alternative_installation_path = fs.join(alternative_temp_path, qt_installation_dirname)
            else: alternative_installation_path = None

            decompress_path = fs.join(temp_path, "Qt-install")
            if alternative_temp_path is not None: alternative_decompress_path = fs.join(alternative_temp_path, "Qt-install")
            else: alternative_decompress_path = None

            if self.remote.is_file(source_path): self.remote.remove_file(source_path)
            if self.remote.is_file(other_source_path): self.remote.remove_file(other_source_path)

            if self.remote.is_directory(decompress_path): self.remote.remove_directory(decompress_path)
            if alternative_decompress_path is not None and self.remote.is_directory(alternative_decompress_path): self.remote.remove_directory(alternative_decompress_path)

            if self.remote.is_directory(installation_path): self.remote.remove_directory(installation_path)
            if alternative_installation_path is not None and self.remote.is_directory(alternative_installation_path): self.remote.remove_directory(alternative_installation_path)

        # Qt is present
        elif self.modules.has_qmake:

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

        if not self.modules.has_git: raise RuntimeError("Git could not be found")

        # Get info
        path = self.modules.paths["git"]
        version = self.modules.versions["git"]
        module = self.modules.names["git"]

        # Set path and module name
        self.git_path = path
        if module is not None: self.git_module = module

    # -----------------------------------------------------------------

    def install_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing Qt ...")

        # Load compiler modules
        if self.compiler_module is not None: self.remote.load_module(self.compiler_module, show_output=log.is_debug)
        if self.mpi_compiler_module is not None and self.mpi_compiler_module != self.compiler_module: self.remote.load_module(self.mpi_compiler_module, show_output=log.is_debug)

        # Determine location based on how much space is available
        if self.remote.free_space_home_directory < full_qt_installation_size:

            # Debugging
            log.debug("Not enough space on the home directory, so installing Qt on scratch system ...")

            if self.remote.scratch_path is None: raise RuntimeError("Not enough space on home directory for installing Qt and scratch path not defined")
            else:
                temp_path = self.remote.scratch_path
                installation_path = fs.join(self.remote.scratch_path, qt_installation_dirname)
        else:

            # Debugging
            log.debug("Installing Qt in the home directory ...")

            temp_path = self.remote.home_directory
            installation_path = fs.join(self.remote.home_directory, qt_installation_dirname)

        # Get the extension of the file
        extension = fs.get_extension(self.qt_url, double=True)

        # Determine the path for the Qt source code
        path = fs.join(temp_path, "qt." + extension)

        # Download Qt
        if not self.remote.is_file(path): self.remote.download_from_url_to(self.qt_url, path)

        # Unarchive
        decompress_path = self.remote.create_directory_in(temp_path, "Qt-install", show_output=log.is_debug)

        if not self.remote.is_directory(decompress_path) or self.remote.is_empty(decompress_path):

            self.remote.decompress_file(path, decompress_path)

            if not self.remote.is_directory(installation_path) or self.remote.is_empty(installation_path):

                # Get the only directory in the Qt-install directory
                qt_everywhere_opensource_path = self.remote.directories_in_path(decompress_path)[0]

                # Take a copy of the Qt configure options
                configure_options = qt_configure_options[:]

                # Set the installation path
                configure_options[0] = configure_options[0].replace("$INSTALLATION_PATH", installation_path)

                # Determine commands
                configure_command = "./configure " + " ".join(configure_options)
                make_command = "make"
                install_command = "make install"

                # Execute the commands
                self.remote.execute_lines(configure_command, make_command, install_command, show_output=log.is_debug, cwd=qt_everywhere_opensource_path)

        # Remove decompressed folder and the tar.gz file
        if self.remote.is_file(path): self.remote.remove_file(path)
        if self.remote.is_directory(decompress_path): self.remote.remove_directory(decompress_path)

        # Success
        log.success("Qt was succesfully installed")

        # Set the qmake path
        all_files_in_installation_directory = self.remote.files_in_path(installation_path, recursive=True)
        for path in all_files_in_installation_directory:
            name = fs.name(path)
            if name == "qmake" and self.remote.is_file(path):
                self.qmake_path = path
                break

        # Check if qmake path was set
        if self.qmake_path is None: raise RuntimeError("Qt was installed but qmake could not be located")

        # Unload modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def get_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the SKIRT source code ...")

        # Load git module
        if self.git_module is not None: self.remote.load_module(self.git_module, show_output=log.is_debug)

        # Get repository link
        if self.config.repository is not None: url = introspection.skirt_git_remote_url(self.config.repository)
        elif self.config.private: url = introspection.private_skirt_https_link
        else: url = introspection.public_skirt_https_link

        # Do HPC UGent in a different way because it seems only SSH is permitted and not HTTPS (but we don't want SSH
        # because of the private/public key thingy, so use a trick
        # DON'T NEED THIS AFTER ALL! WE CAN JUST ADD THE USERNAME AND PASSWORD TO THE HTTPS LINK AND IT WORKS!
        #if self.remote.host.name == "login.hpc.ugent.be":
        #    self.git_version = get_skirt_hpc(self.remote, url, self.skirt_root_path, self.skirt_repo_path)
        #else:

        # Decompose
        host, user_or_organization, repo_name, _, _ = git.decompose_repo_url(url)

        # Find the account file for the repository host (e.g. github.ugent.be)
        if introspection.has_account(host): username, password = introspection.get_account(host)
        else: username = password = None

        # Compose HTTPS link
        url = git.compose_https(host, user_or_organization, repo_name, username, password)

        # Set the clone command
        command = "git clone " + url + " " + self.skirt_repo_path

        # Clone
        self.remote.execute(command, show_output=log.is_debug)

        # Get the git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path, self.remote)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Determine SKIRT and FitSKIRT main paths
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        #fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")

        # Add
        comment = "For SKIRT, added by PTS (Python Toolkit for SKIRT)"
        if skirt_main_path in self.remote.paths_in_path_variable: self.remote.remove_from_path_variable(skirt_main_path)
        #if fitskirt_main_path in self.remote.paths_in_path_variable: self.remote.remove_from_path_variable(fitskirt_main_path)
        self.remote.add_to_path_variable(skirt_main_path, comment=comment, in_shell=True)
        #self.remote.add_to_path_variable(fitskirt_main_path, comment=comment, in_shell=True)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(skirt_main_path, "skirt")

        # Success
        log.success("SKIRT was successfully downloaded")

        # Unload all modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def build_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building SKIRT ...")

        # Load modules
        if self.qmake_module is not None: self.remote.load_module(self.qmake_module, show_output=log.is_debug)
        else:
            if self.compiler_module is not None: self.remote.load_module(self.compiler_module, show_output=log.is_debug)
            if self.mpi_compiler_module is not None and self.mpi_compiler_module != self.compiler_module: self.remote.load_module(self.mpi_compiler_module, show_output=log.is_debug)

        # Execute the build
        build_skirt_on_remote(self.remote, self.skirt_repo_path, self.qmake_path, self.git_version)

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

        output = terminal.execute("skirt -h", show_output=log.is_debug)
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were succesfully installed locally")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation on remote host '" + self.remote.host_id + "'...")

        # Load modules
        if self.qmake_module is not None: self.remote.load_module(self.qmake_module, show_output=log.is_debug)
        else:
            if self.compiler_module is not None: self.remote.load_module(self.compiler_module, show_output=log.is_debug)
            if self.mpi_compiler_module is not None and self.mpi_compiler_module != self.compiler_module: self.remote.load_module(self.mpi_compiler_module, show_output=log.is_debug)

        output = self.remote.execute("skirt -h", show_output=log.is_debug)
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("SKIRT and its dependencies were succesfully installed on remote host '" + self.remote.host_id + "'")

# -----------------------------------------------------------------

def build_skirt_local(skirt_repo_path, qmake_path, git_version):

    """
    This function ...
    :param skirt_repo_path:
    :param qmake_path:
    :param git_version:
    :return:
    """

    # Import here instead of at module level to accomodate clean python installs
    from ..tools import parallelization

    # Create command strings
    make_make_command = qmake_path + " BuildSKIRT.pro -o ../release/Makefile CONFIG+=release"
    nthreads = parallelization.ncores() * parallelization.nthreads_per_core()
    make_command = "make -j " + str(nthreads) + " -w -C ../release"

    # Debugging
    log.debug("Make commands:")
    log.debug(" 1) " + make_make_command)
    log.debug(" 2) " + make_command)
    log.debug("in directory " + skirt_repo_path)

    # Execute
    output = terminal.execute(make_make_command, show_output=log.is_debug, cwd=skirt_repo_path)
    for line in output:
        if "Error processing project file" in line: raise RuntimeError(line)

    # Overwrite the git version
    git_version_content = 'const char* git_version = " ' + git_version + ' " ;'
    git_version_path = fs.join(skirt_repo_path, "SKIRTmain", "git_version.h")
    fs.write_line(git_version_path, git_version_content)

    # Make
    terminal.execute(make_command, cwd=skirt_repo_path, show_output=log.is_debug)

    # Success
    log.success("SKIRT was successfully built")

# -----------------------------------------------------------------

def build_skirt_on_remote(remote, skirt_repo_path, qmake_path, git_version):

    """
    This function ...
    :param remote:
    :param skirt_repo_path:
    :param qmake_path:
    :param git_version:
    :return:
    """

    # Create command strings
    make_make_command = qmake_path + " BuildSKIRT.pro -o ../release/Makefile CONFIG+=release"
    nthreads = remote.cores_per_socket
    make_command = "make -j " + str(nthreads) + " -w -C ../release"

    # Debugging
    log.debug("Make commands:")
    log.debug(" 1) " + make_make_command)
    log.debug(" 2) " + make_command)
    log.debug("in directory " + skirt_repo_path)

    # Configure
    output = remote.execute(make_make_command, show_output=log.is_debug, cwd=skirt_repo_path)
    for line in output:
        if "Error processing project file" in line: raise RuntimeError(line)

    # Overwrite the git version
    git_version_content = 'const char* git_version = " ' + git_version + ' " ;'
    git_version_path = fs.join(skirt_repo_path, "SKIRTmain", "git_version.h")
    write_command = 'echo "' + git_version_content + '" > ' + git_version_path
    remote.execute(write_command, show_output=log.is_debug)

    # Make
    output = remote.execute(make_command, show_output=log.is_debug, cwd=skirt_repo_path)

    # Check the output
    for line in output:
        if "No targets specified and no makefile found.  Stop." in line:
            raise RuntimeError("Building SKIRT on the remote " + remote.host_id + " failed")

# -----------------------------------------------------------------

# Anaconda 4.2.0

# LINUX

# Python 3.5 version
# 64bit: https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86_64.sh
# 32bit: https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86.sh

# Python 2.7 version
# 64bit: https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
# 32bit: https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86.sh

anaconda_linux_url = "https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh"

# MACOS

# Python 3.5 version
# https://repo.continuum.io/archive/Anaconda3-4.2.0-MacOSX-x86_64.sh

# Python 2.7 version
# https://repo.continuum.io/archive/Anaconda2-4.2.0-MacOSX-x86_64.sh

anaconda_macos_url = "https://repo.continuum.io/archive/Anaconda2-4.2.0-MacOSX-x86_64.sh"

# -----------------------------------------------------------------

# MINICONDA

# LINUX

# Python 3.5 version
# 64bit: https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# 32bit: https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86.sh

# Python 2.7 version
# 64bit: https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
# 32bit: https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86.sh

miniconda_linux_url = "https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"

# MACOS

# Python 3.5 version
# https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# Python 2.7 version
# https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh

miniconda_macos_url = "https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh"

# -----------------------------------------------------------------

# Silent installation of Miniconda for Linux and OS X is a simple as specifying the -b and -p arguments of the bash installer. The following arguments are supported:

# -b, batch mode
# -p, installation prefix/path
# -f, force installation even if prefix -p already exists
# Batch mode assumes that you agree to the license agreement, and it does not edit the .bashrc or .bash_profile files.

# A complete example:

# wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
# bash ~/miniconda.sh -b -p $HOME/miniconda
# export PATH="$HOME/miniconda/bin:$PATH"

# -----------------------------------------------------------------

# These Miniconda installers contain the conda package manager and Python. Once Miniconda is installed, you can use the conda command to install any other packages and create environments, etc. For example:

# $ conda install numpy
# ...
# $ conda create -n py3k anaconda python=3

# There are two variants of the installer: Miniconda is Python 2 based and Miniconda3 is Python 3 based. Note that the choice of which Miniconda is installed only affects the root environment. Regardless of which version of Miniconda you install, you can still install both Python 2.x and Python 3.x environments.

# The other difference is that the Python 3 version of Miniconda will default to Python 3 when creating new environments and building packages. So for instance, the behavior of

# $ conda create -n myenv python
# will be to install Python 2.7 with the Python 2 Miniconda and to install Python 3.5 with the Python 3 Miniconda. You can override the default by explicitly setting python=2 or python=3. It also determines the default value of CONDA_PY when using conda build.

# We have 32-bit Mac OS X binaries available, please contact us for more details at sales@continuum.io.

# Note: If you already have Miniconda or Anaconda installed, and you just want to upgrade, you should not use the installer. Just use conda update. For instance

# $ conda update conda
# will update conda.


###

# Anaconda has all that plus over 720 open source packages that install with Anaconda or can be installed with the simple conda install command.


###

# In your browser download the Miniconda installer for Linux, then in your terminal window type the following and follow the prompts on the installer screens. If unsure about any setting, simply accept the defaults as they all can be changed later:

# bash Miniconda3-latest-Linux-x86_64.sh
# Now close and re-open your terminal window for the changes to take effect.

# To test your installation, enter the command conda list. If installed correctly, you will see a list of packages that were installed.

# -----------------------------------------------------------------

# Montage
montage_url = "http://montage.ipac.caltech.edu/download/Montage_v5.0.tar.gz"

# Imfit source code
imfit_url = "imfit-1.4.0-source.tar.gz"

# Imfit binaries
imfit_macos_binary_url = "http://www.mpe.mpg.de/~erwin/resources/imfit/binaries/imfit-1.4.0-macintel.tar.gz"
imfit_linux_binary_url = "http://www.mpe.mpg.de/~erwin/resources/imfit/binaries/imfit-1.4.0-linux-64.tar.gz"

# -----------------------------------------------------------------

# Approximate size of a full conda installation (all PTS dependencies) in GB
full_conda_installation_size = 5. * u("GB")

# -----------------------------------------------------------------

class PTSInstaller(Installer):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PTSInstaller, self).__init__(*args, **kwargs)

        # Path to python executable
        self.python_path = None

        # Path to PTS root directory
        self.pts_root_path = None

        # Path to PTS/pts
        self.pts_package_path = None

        # Conda
        self.conda_main_executable_path = None
        self.conda_installation_path = None
        self.conda_executable_path = None
        self.conda_pip_path = None
        self.conda_jupyter_path = None
        self.conda_activate_path = None
        self.conda_python_path = None
        self.conda_easy_install_path = None

        # Path to the PTS executable
        self.pts_path = None

        # The git version
        self.git_version = None

        # The path to the montage installation directory
        self.montage_path = None

        # The path to the Imfit installation directory
        self.imfit_path = None

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        # Set root path and pacakge path
        self.pts_root_path = introspection.pts_root_dir
        self.pts_package_path = introspection.pts_package_dir
        self.pts_path = fs.join(self.pts_package_path, "do", "__main__.py")

        # Force flags are not specified
        #if not (self.config.force or self.config.force_conda):

            # Give warning
            #log.warning("You are asking to install PTS locally, but you are already using PTS at the moment")
            #log.warning("Add the '--force' or '--force_conda' option to force PTS installation or conda python distribution installation")
            #log.warning("Quitting ...")
            #exit()

        # Ask to proceed with installing python environment and dependencies
        #definition = ConfigurationDefinition()
        #definition.add_flag("proceed", "Proceed by installing a clean Python 2 environment for PTS and automatically installing all dependencies?", False)
        #setter = InteractiveConfigurationSetter("Install PTS local", add_logging=False, add_cwd=False)
        #config = setter.run(definition, prompt_optional=True)

        # Proceed or quit
        #if config.proceed: log.info("Proceeding ...")
        #else:
        #    log.info("Quitting ...")
        #    exit()

    # -----------------------------------------------------------------

    def create_directories_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the directory structure ...")

        # Set root path and pacakge path
        self.pts_root_path = fs.join(self.remote.home_directory, "PTS")
        self.pts_package_path = fs.join(self.pts_root_path, "pts")

        # Check if already present
        if self.remote.is_directory(self.pts_root_path):
            if self.config.force: self.remote.remove_directory(self.pts_root_path)
            else: raise RuntimeError("PTS is already installed (or partly present) on the remote host")

        # Make the root directory
        self.remote.create_directory(self.pts_root_path)

        # Create the other directories
        for name in pts_directories:

            # Determine path
            path = fs.join(self.pts_root_path, name)

            # Create the directory
            self.remote.create_directory(path)

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

    def install_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing PTS locally ...")

        # 1. Check the presence of conda
        self.check_conda_local()

        # 2. Get conda locally
        if not self.has_conda: self.get_conda_local()

        # 3. Create python environment
        if not has_valid_conda_environment_local(self.config.python_name, self.conda_main_executable_path, self.conda_installation_path): self.create_environment_local()
        else: self.set_environment_paths_local()

        # 4. Get PTS
        if self.config.force: self.get_pts_local()

        # 5. Get PTS dependencies
        self.get_dependencies_local()

    # -----------------------------------------------------------------

    def check_conda_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of a conda installation ...")

        # Find conda locally
        conda_executable_path = conda.find_conda()

        # If conda is present
        if conda_executable_path is not None:

            self.conda_installation_path = fs.directory_of(fs.directory_of(conda_executable_path))
            if self.config.force_conda:
                fs.remove_directory(self.conda_installation_path)
                python_executable_path = sys.executable
                if self.conda_installation_path in python_executable_path:
                    log.warning("Conda is removed for a clean install but this instance of PTS was launched using conda python. Restart the installation procedure.")
                    log.warning("Quitting ...")
                    exit()
            else: self.conda_main_executable_path = conda_executable_path

    # -----------------------------------------------------------------

    def get_conda_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting a Conda python distribution locally ...")

        # Install conda locally and set paths
        self.conda_installation_path, self.conda_main_executable_path = install_conda_local()

    # -----------------------------------------------------------------

    def create_environment_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a fresh python environment ...")

        # Create the environment
        self.conda_executable_path, self.conda_pip_path, self.conda_activate_path, self.conda_python_path, self.conda_easy_install_path = create_conda_environment_local(self.config.python_name, self.conda_installation_path, self.pts_root_path, self.config.python_version, self.conda_main_executable_path)

    # -----------------------------------------------------------------

    def set_environment_paths_local(self):

        """
        This function ...
        :return:
        """

        # Set paths
        environment_bin_path = fs.join(self.conda_installation_path, "envs", self.config.python_name, "bin")
        if not fs.is_directory(environment_bin_path): raise RuntimeError("Creating the environment failed")
        self.conda_executable_path = fs.join(environment_bin_path, "conda")
        self.conda_pip_path = fs.join(environment_bin_path, "pip")
        self.conda_jupyter_path = fs.join(environment_bin_path, "jupyter")
        self.conda_activate_path = fs.join(environment_bin_path, "activate")
        self.conda_python_path = fs.join(environment_bin_path, "python")
        self.conda_easy_install_path = fs.join(environment_bin_path, "easy_install")

        # Check if paths exist
        assert fs.is_file(self.conda_executable_path)
        assert fs.is_file(self.conda_pip_path)
        assert fs.is_file(self.conda_activate_path)
        assert fs.is_file(self.conda_python_path)
        assert fs.is_file(self.conda_easy_install_path)

        # Jupyter: optional
        if not fs.is_file(self.conda_jupyter_path): self.conda_jupyter_path = None

        # Setup the environment
        # environment_name, pip_name, pts_root_path, python_path
        # environment_name, pip_name, jupyter_name, pts_root_path, python_path, pip_path, jupyter_path
        setup_conda_environment_local(self.config.python_name, self.config.pip_name, self.config.jupyter_name, self.pts_root_path, self.conda_python_path, self.conda_pip_path, self.conda_jupyter_path)

    # -----------------------------------------------------------------

    def get_pts_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading PTS ...")

        # Determine repo url
        if self.config.repository is not None: raise ValueError("Repository cannot be specified for local installation")
        elif self.config.private: url = introspection.private_pts_https_link
        else: url = introspection.public_pts_https_link

        # Get host (github.ugent.be)
        host = url.split("//")[1].split("/")[0]

        # Set the clone command
        command = "git clone " + url + " " + self.pts_package_path

        # Find the account file for the repository host (e.g. github.ugent.be)
        if introspection.has_account(host):

            username, password = introspection.get_account(host)

            # Set the command lines
            lines = []
            lines.append(command)
            lines.append(("':", username))
            lines.append(("':", password))

            # Clone the repository
            terminal.execute_lines(*lines, show_output=log.is_debug)

        else: terminal.execute(command, show_output=log.is_debug)

        # Get the git version
        self.git_version = git.get_short_git_version(self.pts_package_path)

        # Show the git version
        log.info("The version that was installed is '" + self.git_version + "'")

        # Set the path to the main PTS executable
        self.pts_path = fs.join(self.pts_package_path, "do", "__main__.py")

    # -----------------------------------------------------------------

    def get_dependencies_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the PTS dependencies ...")

        # Get available conda packages
        #output = terminal.execute_no_pexpect(self.conda_executable_path + " search", show_output=log.is_debug)
        output = terminal.execute_no_pexpect(self.conda_executable_path + " search")
        available_packages = []
        for line in output:
            if not line.split(" ")[0]: continue
            available_packages.append(line.split(" ")[0])

        # Get already installed packages
        already_installed = []
        #for line in terminal.execute_no_pexpect(self.conda_executable_path + " list --name " + self.config.python_name, show_output=log.is_debug):
        for line in terminal.execute_no_pexpect(self.conda_executable_path + " list --name " + self.config.python_name):
            if line.startswith("#"): continue
            already_installed.append(line.split(" ")[0])

        # Use the introspection module on the remote end to get the dependencies and installed python packages
        # session = self.remote.start_python_session()
        # session.import_package("introspection", from_name="pts.core.tools")

        dependencies = introspection.get_all_dependencies().keys()

        from ..tools import stringify

        # Debugging
        log.debug("")
        log.debug("Dependencies:")
        log.debug("")
        for line in stringify.stringify_list_fancy(dependencies, lines_prefix="   ")[1].split("\n"): log.debug(line)
        log.debug("")

        #packages = introspection.installed_python_packages()
        packages = []

        # Inform the user
        log.info("Activating the '" + self.config.python_name + "' conda environment ...")

        # Activate the correct environment
        previous_environment = conda.activate_environment(self.config.python_name, self.conda_executable_path, self.conda_activate_path)

        # Get installation commands
        installation_commands, installed, not_installed, std_lib, real_names = get_installation_commands(dependencies, packages,
                                                                                    already_installed, available_packages,
                                                                                    conda_path=self.conda_executable_path,
                                                                                    pip_path=self.conda_pip_path,
                                                                                    conda_environment=self.config.python_name,
                                                                                    python_path=self.conda_python_path,
                                                                                    easy_install_path=self.conda_easy_install_path)

        # Show commands
        log.debug("")
        if len(installation_commands) > 0:
            log.debug("Installation commands:")
            log.debug("")
            for module in installation_commands:
                if isinstance(installation_commands[module], list): string = installation_commands[module][0]
                else: string = installation_commands[module]
                log.debug(" - " + module + ": " + string)
            log.debug("")

        # Debugging
        log.debug("In standard library:")
        log.debug("")
        for line in stringify.stringify_list_fancy(std_lib, lines_prefix="   ")[1].split("\n"): log.debug(line)
        log.debug("")

        # Debugging
        log.debug("Already installed:")
        log.debug("")
        for line in stringify.stringify_list_fancy(already_installed, lines_prefix="   ")[1].split("\n"): log.debug(line)
        log.debug("")

        # Debugging
        log.debug("Cannot be installed:")
        log.debug("")
        for line in stringify.stringify_list_fancy(not_installed, lines_prefix="   ")[1].split("\n"): log.debug(line)
        log.debug("")

        # Internal check
        for dependency in dependencies:
            #print(real_names)
            if dependency in real_names: dependency = real_names[dependency]
            if dependency in installed: continue
            if dependency in not_installed: continue
            if dependency in already_installed: continue
            if dependency in std_lib: continue
            log.error("Dependency '" + dependency + "' not in 'installed', 'not_installed', 'already_installed', or 'std_lib'")

        # Install
        for module in installation_commands:

            # Inform the user
            log.info("Installing '" + module + "' ...")

            # Check which command can be used
            command = installation_commands[module]

            # Debugging
            if isinstance(command, list): log.debug("Installation command: '" + command[0] + "'")
            elif types.is_string_type(command): log.debug("Installation_command: '" + command + "'")
            else: raise ValueError("Invalid installation command: " + str(command))

            # Launch the command
            import subprocess
            try:
                if isinstance(command, list):
                    # Skip module if already installed together with another package in the meantime
                    if command[0].startswith(self.conda_executable_path):
                        if conda.is_present_package(module, self.config.python_name, self.conda_executable_path): continue
                    output = terminal.execute_lines_expect_clone(*command, show_output=log.is_debug)
                elif types.is_string_type(command):
                    if "setup.py" in command:
                        # special: for python setup.py, we must be in the directory or it won't work
                        dir_path = fs.directory_of(command.split()[1])
                        setup_path = fs.join(dir_path, "setup.py")
                        command.replace(setup_path, "setup.py")
                        output = terminal.execute(command, show_output=log.is_debug, cwd=dir_path)
                    else: output = terminal.execute_no_pexpect(command, show_output=log.is_debug)
                else: raise ValueError("Invalid installation command: " + str(command))

                # Check the output
                for line in output:
                    if "Exception:" in line:
                        failed = True
                        break
                else: failed = False

                if failed:
                    installed.remove(module)
                    log.warning("Something went wrong installing '" + module + "'")
                    for line in output:
                        log.warning(line)
                    not_installed.append(module)
                else: log.success("Installation of '" + module + "' was succesful")

            # Failed
            except subprocess.CalledProcessError:

                installed.remove(module)
                log.warning("Something went wrong installing '" + module + "'")
                traceback.print_exc()
                not_installed.append(module)

        from ..tools import stringify

        # Show installed packages
        if len(installed) > 0:
            log.info("Packages that were installed:")
            print("")
            print(stringify.stringify_list_fancy(installed, width=100, delimiter=", ", lines_prefix="    ")[1])
            print("")

        # Show not installed packages
        if len(not_installed) > 0:
            log.info("Packages that could not be installed:")
            print("")
            print(stringify.stringify_list_fancy(not_installed, width=100, delimiter=", ", lines_prefix="    ")[1])
            print("")

        # Show already present packages
        if len(already_installed) > 0:
            log.info("Packages that were already present:")
            print("")
            print(stringify.stringify_list_fancy(already_installed, width=100, delimiter=", ", lines_prefix="    ")[1])
            print("")

        # Inform he user
        log.info("Deactivating the conda environment ...")

        # Activate the previous environment
        conda.activate_environment(previous_environment, self.conda_executable_path, self.conda_activate_path)

        # Deactivate
        conda.deactivate(self.conda_activate_path.replace("activate", "deactivate"))

    # -----------------------------------------------------------------

    def install_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing PTS remotely ...")

        # Check the presence of imfit
        self.check_imfit_remote()

        # Install imfit
        if not self.has_imfit: self.get_imfit_remote()

        # Check the presence of Montage
        self.check_montage()

        # Install montage
        if not self.has_montage: self.get_montage()

        # Check presence of conda remotely
        self.check_conda_remote()

        # 1. Get the conda python distribution
        if not self.has_conda: self.get_conda_remote()

        # 2. Create a remote python environment
        if not has_valid_conda_environment_remote(self.remote, self.config.python_name, self.conda_main_executable_path, self.conda_installation_path): self.create_environment_remote()
        else: self.set_environment_paths_remote()

        # 3. Get PTS
        self.get_pts_remote()

        # 4. Get PTS dependencies
        self.get_dependencies_remote()

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
        filepath = self.remote.download_from_url_to(url, self.remote.home_directory, show_output=log.is_debug)

        # Decompress into "~/Imfit"
        imfit_installation_path = fs.join(self.remote.home_directory, "Imfit")
        if self.remote.is_directory(imfit_installation_path): raise RuntimeError("There is already a directory '" + imfit_installation_path + "'")
        else: self.remote.create_directory(imfit_installation_path)
        self.remote.decompress_directory_to(filepath, imfit_installation_path, show_output=log.is_debug, remove=True)

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
        filepath = self.remote.download_from_url_to(montage_url, self.remote.home_directory, show_output=log.is_debug)

        # Extract the file
        montage_path = self.remote.decompress_directory_in_place(filepath, show_output=log.is_debug, remove=True)

        # Compile
        self.remote.execute("make", cwd=montage_path, show_output=log.is_debug)

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
        self.conda_installation_path, conda_main_executable_path = self.remote.find_conda()

        # If force installation
        if self.conda_installation_path is not None and self.config.force_conda: self.remote.remove_directory(self.conda_installation_path)
        else: self.conda_main_executable_path = conda_main_executable_path

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
        # # conda_executable_path, conda_pip_path, conda_activate_path, conda_python_path, conda_easy_install_path
        self.conda_executable_path, self.conda_pip_path, self.conda_activate_path, self.conda_python_path, \
        self.conda_easy_install_path, self.conda_jupyter_path = create_conda_environment_remote(self.remote, self.config.python_name,
                                                                       self.conda_installation_path, self.pts_root_path,
                                                                       self.config.python_version,
                                                                       self.conda_main_executable_path)

    # -----------------------------------------------------------------

    def set_environment_paths_remote(self):

        """
        This function ...
        :return:
        """

        # Set paths
        environment_bin_path = fs.join(self.conda_installation_path, "envs", self.config.python_name, "bin")
        if not self.remote.is_directory(environment_bin_path): raise RuntimeError("Creating the environment failed")
        self.conda_executable_path = fs.join(environment_bin_path, "conda")
        self.conda_pip_path = fs.join(environment_bin_path, "pip")
        self.conda_jupyter_path = fs.join(environment_bin_path, "jupyter")
        self.conda_activate_path = fs.join(environment_bin_path, "activate")
        self.conda_python_path = fs.join(environment_bin_path, "python")
        self.conda_easy_install_path = fs.join(environment_bin_path, "easy_install")

        # Check if paths exist
        assert self.remote.is_file(self.conda_executable_path)
        assert self.remote.is_file(self.conda_pip_path)
        assert self.remote.is_file(self.conda_activate_path)
        assert self.remote.is_file(self.conda_python_path)
        #assert self.remote.is_file(self.conda_easy_install_path)
        if not self.remote.is_file(self.conda_easy_install_path): log.warning("easy_install executable is not present")

        # Jupyter is not required
        if not self.remote.is_file(self.conda_jupyter_path): self.conda_jupyter_path = None

        # Setup
        # remote, environment_name, pip_name, pts_root_path, python_path, pip_path
        setup_conda_environment_remote(self.remote, self.config.python_name, self.config.pip_name, self.config.jupyter_name, self.pts_root_path, self.conda_python_path, self.conda_pip_path, self.conda_jupyter_path)

    # -----------------------------------------------------------------

    def get_pts_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading PTS on the remote host ...")

        if self.config.repository is not None: url = introspection.pts_git_remote_url(self.config.repository)
        elif self.config.private: url = introspection.private_pts_https_link
        else: url = introspection.public_pts_https_link

        # Decompose
        host, user_or_organization, repo_name, _, _ = git.decompose_repo_url(url)

        # Find the account file for the repository host (e.g. github.ugent.be)
        if introspection.has_account(host): username, password = introspection.get_account(host)
        else: username = password = None

        # Compose url
        url = git.compose_https(host, user_or_organization, repo_name, username, password)

        # Set the clone command
        command = "git clone " + url + " " + self.pts_package_path

        # Clone
        self.remote.execute(command, show_output=log.is_debug)

        # Get the git version
        self.git_version = git.get_short_git_version(self.pts_package_path, self.remote)

        # Show the git version
        log.info("The version that was installed is '" + self.git_version + "'")

        # Set the path to the main PTS executable
        self.pts_path = fs.join(self.pts_package_path, "do", "__main__.py")

    # -----------------------------------------------------------------

    def get_dependencies_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the PTS dependencies on the remote host ...")

        # Install PTS dependencies
        installed, not_installed, already_installed = get_pts_dependencies_remote(self.remote, self.pts_package_path, conda_path=self.conda_executable_path, pip_path=self.conda_pip_path,
                                    python_path=self.conda_python_path, easy_install_path=self.conda_easy_install_path,
                                    conda_environment=self.config.python_name, conda_activate_path=self.conda_activate_path)

        # Success
        log.success("Succesfully installed the dependencies on the remote host")

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation locally ...")

        output = terminal.execute_no_pexpect(self.pts_path + " -h", show_output=log.is_debug)
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were succesfully installed locally")

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation on remote host '" + self.remote.host_id + "'...")

        output = self.remote.execute(self.pts_path + " -h", show_output=log.is_debug)
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

        # Success
        log.success("PTS and its dependencies were succesfully installed on remote host '" + self.remote.host_id + "'")

# -----------------------------------------------------------------

def find_real_name(module_name, available_packages, real_names):

    """
    This function ...
    :param module_name:
    :param available_packages:
    :param real_names:
    :return:
    """

    if module_name in available_packages: return module_name, "conda", None
    if module_name in real_names: return real_names[module_name], "conda", None

    # Search on pypi
    from pts.core.tools.pypi import search
    results = list(search(module_name))
    in_pip = False
    version = None
    for result in results:
        if result["name"] == module_name:
            in_pip = True
            version = result["version"]
            break
    if in_pip: return module_name, "pip", version

    if "_" in module_name:

        new_name = module_name.replace("_", "-")
        results = list(search(new_name))
        #print(results)
        in_pip = False
        version = None
        for result in results:
            if result["name"] == new_name:
                in_pip = True
                version = result["version"]
                break
        if in_pip: return new_name, "pip", version

    # Not found
    return None, None, None

# -----------------------------------------------------------------

def get_installation_commands(dependencies, packages, already_installed, available_packages, conda_path="conda",
                              pip_path="pip", python_path="python", easy_install_path="easy_install",
                              conda_environment=None, remote=None, check_by_importing=True, repositories=None):

    """
    This function ...
    :param dependencies:
    :param packages:
    :param already_installed:
    :param available_packages:
    :param conda_path:
    :param pip_path:
    :param python_path:
    :param easy_install_path:
    :param conda_environment:
    :param remote:
    :param check_by_importing:
    :param repositories:
    :return:
    """

    # Get dependency version restrictions
    versions = introspection.get_constricted_versions()

    # Get names of packages for import names
    real_names = introspection.get_package_names()

    # Get repositories for import names
    if repositories is None: repositories = introspection.get_package_repositories()

    installed = []
    not_installed = []
    std_lib = []

    commands = dict()

    # Loop over the dependencies
    for module in dependencies:

        # Debugging
        log.debug("Checking dependency '" + module + "' ...")

        module_name = module

        if module_name in packages: continue

        # Skip packages from the standard library
        if introspection.is_std_lib(module_name):
            std_lib.append(module_name)
            continue

        # Check if already installed
        if module_name in already_installed: continue
        if module_name in real_names and real_names[module_name] in already_installed: continue

        # Double check
        if module_name.replace("_", "-") in already_installed:
            already_installed.append(module_name)
            continue

        # Check if a repository link is defined for this package
        if module_name in repositories:

            # Check explicitly by trying to import the module
            if check_by_importing and can_import_module(module_name, remote, python_path):
                already_installed.append(module_name)
                continue

            via = repositories[module_name]
            install_module_name = module_name
            version = None

        else:

            # Find name, check if available
            install_module_name, via, version = find_real_name(module_name, available_packages, real_names)

            # Triple check
            if install_module_name in already_installed:
                already_installed.append(module_name)
                continue

            # Check explicitly by trying to import the module
            if check_by_importing and can_import_module(module_name, remote, python_path):
                already_installed.append(module_name)
                continue

            # Module is not found
            if install_module_name is None:

                log.warning("Package '" + module + "' can not be installed")
                not_installed.append(module)
                continue

            # Set real name
            else: real_names[module_name] = install_module_name

        # Debugging
        log.debug("Checking whether a specific version of the package is required ...")

        # Checking the version restriction
        if module_name in versions:
            version = versions[module_name]
            log.debug("Version '" + version + "' is required")

        # Debugging
        log.debug("Determining installation command ...")

        # Installable via conda
        if via == "conda":

            if conda_environment is not None: command = conda_path + " install --name " + conda_environment + " " + install_module_name
            else: command = conda_path + " install " + install_module_name
            if version is not None: command += "=" + version

            lines = []
            lines.append(command)
            lines.append(("Proceed ([y]/n)?", "y", True))

            # Set the commands
            commands[install_module_name] = lines

            # Add to installed
            installed.append(install_module_name)

        # Install with pip
        elif via == "pip":

            command = pip_path + " install " + install_module_name
            if version is not None: command += "==" + version
            commands[install_module_name] = command

            # Add to installed
            installed.append(install_module_name)

        # GitHub url
        elif "github.com" in via:

            #log.warning("Obtaining packages from a remote repository on GitHub is not supported yet")
            #not_installed.append(module_name)
            #continue

            command = pip_path + " install git+" + via
            commands[install_module_name] = command

            # Add to installed
            installed.append(install_module_name)

        # Link with source code
        elif via.startswith("http"):

            # Determine download path
            from ..tools import time
            if remote is not None: path = remote.create_directory_in(remote.pts_temp_path, time.unique_name(module_name))
            else: path = introspection.create_temp_dir(time.unique_name(module_name))

            # Download
            if remote is not None: filepath = remote.download_from_url_to(via, path)
            else: filepath = network.download_file_no_requests(via, path)

            # Debugging
            log.debug("Decompressing the file '" + filepath + "' ...")

            # Decompress
            if remote is not None: decompress_directory_path = remote.decompress_directory_in_place(filepath, remove=True)
            else: decompress_directory_path = archive.decompress_directory_in_place(filepath, remove=True)

            # Set the installation command
            #setup_path = fs.join(decompress_directory_path, "setup.py")
            #command = python_path + " " + setup_path + " install"

            # Use pip install ./directory
            command = pip_path + " install " + decompress_directory_path
            commands[install_module_name] = command

            # Add to installed
            installed.append(install_module_name)

        # Not recognized
        else: not_installed.append(install_module_name)

    # Return ...
    return commands, installed, not_installed, std_lib, real_names

# -----------------------------------------------------------------

def find_qmake():

    """
    This function ...
    :return:
    """

    # Translated from ./makeSKIRT.sh

    # Get a list of qmake paths installed on this system
    qmake_paths = []

    for qt_dir in fs.directories_in_path(fs.home, startswith="Qt"):
        qmake_paths = fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

    for qt_dir in fs.directories_in_path("/usr/local", startswith="Qt"):
        qmake_paths += fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

    qmake_path = introspection.qmake_path()
    qmake_paths += [qmake_path] if qmake_path is not None else []

    # Get the most recent installation

    latest_qmake_path = None

    latest_version = None

    # Loop over the qmake paths: SAME IMPLEMENTATION AS IN Remote class -> _check_qt_remote_no_lmod
    for qmake_path in qmake_paths:

        # Debugging
        log.debug("Found a qmake executable at '" + qmake_path + "'")

        # Get the version
        output = terminal.execute(qmake_path + " -v", show_output=log.is_debug)
        for line in output:
            if "Qt version" in line:
                qt_version = line.split("Qt version ")[1].split(" in")[0]
                break
        else: raise RuntimeError("Qt version could not be determined")

        if qt_version < "5.2.0": continue  # oldest supported version
        if "conda" in qmake_path: continue
        if "canopy" in qmake_path: continue
        if "epd" in qmake_path: continue
        if "enthought" in qmake_path: continue

        if latest_version is None or qt_version > latest_version:
            latest_version = qt_version
            latest_qmake_path = qmake_path

    # Return the path
    return latest_qmake_path

# -----------------------------------------------------------------

def get_pts_dependencies_remote(remote, pts_package_path, conda_path="conda", pip_path="pip", python_path="python",
                                easy_install_path="easy_install", conda_environment=None, conda_activate_path="activate"):

    """
    This fucntion ...
    :param remote:
    :param pts_package_path:
    :param conda_path:
    :param pip_path:
    :param python_path:
    :param easy_install_path:
    :param conda_environment:
    :param conda_activate_path:
    :return:
    """

    # Get available conda packages
    available_packages = remote.available_conda_packages(conda_path)

    # Get already installed packages
    already_installed = []

    #print("Conda environment: " + str(conda_environment))

    # List
    if conda_environment is not None: list_command = conda_path + " list --name " + conda_environment
    else: list_command = conda_path + " list"
    #for line in remote.execute(list_command, show_output=log.is_debug):
    for line in remote.execute(list_command):
        if line.startswith("#"): continue
        already_installed.append(line.split(" ")[0])

    # Debugging
    #log.debug("Already installed packages:")
    #for package in already_installed: log.debug(" - " + package)

    ## Install essential packages: numpy and astropy
    for essential in ["pip", "numpy", "astropy"]:

        # Skip installation if it is already present
        if essential in already_installed: continue

        # Install
        if conda_environment is not None: command = conda_path + " install " + essential + " --name " + conda_environment
        else: command = conda_path + " install " + essential
        lines = []
        lines.append(command)
        lines.append(("Proceed ([y]/n)?", "y", True))
        remote.execute_lines(*lines, show_output=log.is_debug)

    # Use the introspection module on the remote end to get the dependencies and installed python packages
    screen_output_path = remote.create_directory_in(remote.pts_temp_path, time.unique_name("installation"))
    #session = remote.start_python_session(output_path=screen_output_path)
    #session.import_package("introspection", from_name="pts.core.tools")
    #dependencies = session.get_simple_property("introspection", "get_all_dependencies().keys()")

    # Get list of dependencies on the remote host
    dependencies_script_path = fs.join(pts_package_path, "dependencies.py")
    dependencies = remote.execute("python " + dependencies_script_path)[1:] # FIRST LINE = "10/07/2017 17:25:59.973   Welcome to PTS"

    from ..tools import stringify

    # Debugging
    log.debug("")
    log.debug("Dependencies:")
    log.debug("")
    for line in stringify.stringify_list_fancy(dependencies, lines_prefix="   ")[1].split("\n"): log.debug(line)
    log.debug("")

    #packages = session.get_simple_property("introspection", "installed_python_packages()")
    packages = []

    # self.remote.end_python_session()
    # Don't end the python session just yet

    # Inform the user
    log.info("Activating the '" + conda_environment + "' conda environment ...")

    # Change the conda environment
    previous_environment = remote.activate_conda_environment(conda_environment, conda_path, conda_activate_path, show_output=log.is_debug)

    # Get installation commands
    # dependencies, packages, already_installed, available_packages, conda_path="conda",
    # pip_path="pip", python_path="python", easy_install_path="easy_install",
    # conda_environment=None
    installation_commands, installed, not_installed, std_lib, real_names = get_installation_commands(dependencies, packages,
                                                                                already_installed, available_packages,
                                                                                conda_path=conda_path,
                                                                                pip_path=pip_path,
                                                                                conda_environment=conda_environment,
                                                                                python_path=python_path,
                                                                                easy_install_path=easy_install_path,
                                                                                remote=remote)

    # Stop the python session
    #del session

    log.debug("")
    if len(installation_commands) > 0:
        log.debug("Installation commands:")
        log.debug("")
        for module in installation_commands:
            if isinstance(installation_commands[module], list): string = installation_commands[module][0]
            else: string = installation_commands[module]
            log.debug(" - " + module + ": " + string)
        log.debug("")

    # Debugging
    log.debug("In standard library:")
    log.debug("")
    for line in stringify.stringify_list_fancy(std_lib, lines_prefix="   ")[1].split("\n"): log.debug(line)
    log.debug("")

    # Debugging
    log.debug("Already installed:")
    log.debug("")
    for line in stringify.stringify_list_fancy(already_installed, lines_prefix="   ")[1].split("\n"): log.debug(line)
    log.debug("")

    # Debugging
    log.debug("Cannot be installed:")
    log.debug("")
    for line in stringify.stringify_list_fancy(not_installed, lines_prefix="   ")[1].split("\n"): log.debug(line)
    log.debug("")

    # Internal check
    for dependency in dependencies:
        #print(real_names)
        if dependency in real_names: dependency = real_names[dependency]
        if dependency in installed: continue
        if dependency in not_installed: continue
        if dependency in already_installed: continue
        if dependency in std_lib: continue
        log.error("Dependency '" + dependency + "' not in 'installed', 'not_installed', 'already_installed', or 'std_lib'")

    # Install
    for module in installation_commands:

        # Inform the user
        log.info("Installing '" + module + "' ...")

        command = installation_commands[module]

        # Debugging
        if isinstance(command, list): log.debug("Installation command: '" + command[0] + "'")
        elif types.is_string_type(command): log.debug("Installation_command: '" + command + "'")
        else: raise ValueError("Invalid installation command: " + str(command))

        # Install remotely
        # remote, command, conda_path, module, conda_environment
        result = install_module_remote(remote, command, conda_path, module, conda_environment)

        # Fail, stack trace back
        if isinstance(result, list):

            installed.remove(module)
            log.warning("Something went wrong installing '" + module + "'")
            for line in result: print(line)
            not_installed.append(module)

        # Success
        elif result: log.success("Installation of '" + module + "' was succesful")

        # Already installed together with another module
        else: log.success("Module '" + module + "' was already installed")

    # Show installed packages
    if len(installed) > 0:
        log.info("Packages that were installed:")
        print("")
        print(stringify.stringify_list_fancy(installed, width=100, delimiter=", ", lines_prefix="    ")[1])
        print("")

    # Show not installed packages
    if len(not_installed) > 0:
        log.info("Packages that could not be installed:")
        print("")
        print(stringify.stringify_list_fancy(not_installed, width=100, delimiter=", ", lines_prefix="    ")[1])
        print("")

    # Show already present packages
    if len(already_installed) > 0:
        log.info("Packages that were already present:")
        print("")
        print(stringify.stringify_list_fancy(already_installed, width=100, delimiter=", ", lines_prefix="    ")[1])
        print("")

    # Inform he user
    log.info("Deactivating the conda environment ...")

    # Change the environment back
    remote.activate_conda_environment(previous_environment, conda_path, conda_activate_path)

    # Deactivate
    remote.deactivate_conda(conda_activate_path.replace("activate", "deactivate"))

    # Return
    return installed, not_installed, already_installed

# -----------------------------------------------------------------

def install_conda_local():

    """
    This function ...
    :return:
    """

    # Determine installer path
    installer_path = fs.join(fs.home, "conda.sh")

    # Debugging
    log.debug("Downloading the installer ...")

    # Download the installer
    if introspection.is_macos(): network.download_file_no_requests(miniconda_macos_url, installer_path, overwrite=True)
    elif introspection.is_linux(): network.download_file_no_requests(miniconda_linux_url, installer_path, overwrite=True)
    else: raise OSError("Your operating system is not supported")

    # Determine the conda installation directory
    conda_installation_path = fs.join(fs.home, "miniconda")

    # Check if not present
    if fs.is_directory(conda_installation_path): raise RuntimeError("Conda is already installed in '" + conda_installation_path + "'")

    # Determine the options
    options = "-b -p " + conda_installation_path

    # Debugging
    log.debug("Running the installer ...")

    # Run the installation script
    terminal.run_script(installer_path, options, show_output=log.is_debug, no_pexpect=True)

    # Debugging
    log.debug("Removing the installer ...")

    # Remove the installer
    fs.remove_file(installer_path)

    # Set the installation path
    conda_bin_path = fs.join(conda_installation_path, "bin")

    # Set main executable
    conda_main_executable_path = fs.join(conda_installation_path, "bin", "conda")

    # Debugging
    log.debug("Adding the conda executables to the PATH ...")

    # Clear previous things
    comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
    if conda_bin_path in terminal.paths_in_path_variable(): terminal.remove_from_path_variable(conda_bin_path)
    terminal.remove_from_path_variable_containing("miniconda/bin")

    # Run the export command also in the current shell, so that the conda commands can be found
    terminal.add_to_path_variable(conda_bin_path, comment=comment, in_shell=True)

    # Return paths
    return conda_installation_path, conda_main_executable_path

# -----------------------------------------------------------------

def install_conda_remote(remote):

    """
    This function ...
    :param remote:
    :return:
    """

    # Determine the conda installation location based on how much space is available
    if remote.free_space_home_directory < full_conda_installation_size:

        # Debugging
        log.debug("Not enough space on the home directory, so installing Conda on scratch system ...")

        if remote.scratch_path is None:
            raise RuntimeError("Not enough space on home directory for installing Qt and scratch path not defined")
        else:
            installer_path = fs.join(remote.scratch_path, "conda.sh")
            conda_installation_path = fs.join(remote.scratch_path, "miniconda")
    else:

        # Debugging
        log.debug("Installing Conda in the home directory ...")

        installer_path = fs.join(remote.home_directory, "conda.sh")
        conda_installation_path = fs.join(remote.home_directory, "miniconda")

    # Download the installer
    if remote.is_macos: remote.download_from_url_to(miniconda_macos_url, installer_path, overwrite=True)
    elif remote.is_linux: remote.download_from_url_to(miniconda_linux_url, installer_path, overwrite=True)
    else: raise OSError("The operating system on the remote host is not supported")

    # Cleanup leftovers
    if remote.is_directory(conda_installation_path): remote.remove_directory(conda_installation_path)

    # Set options for installer
    options = "-b -p " + conda_installation_path

    # Debugging
    log.debug("Running the '" + installer_path + "' script with options: [" + options + "]")

    # Run the installation script
    remote.run_script(installer_path, options, show_output=log.is_debug)

    # Debugging
    log.debug("Removing the installer script ...")

    # Remove the installer
    remote.remove_file(installer_path)

    # Set the installation path
    conda_installation_path = conda_installation_path
    conda_bin_path = fs.join(conda_installation_path, "bin")

    # Set main executable
    conda_main_executable_path = fs.join(conda_installation_path, "bin", "conda")

    # Debugging
    log.debug("Adding the conda executables to the PATH ...")

    # Remove previous things
    comment = "For Conda, added by PTS (Python Toolkit for SKIRT)"
    if conda_bin_path in remote.paths_in_path_variable: remote.remove_from_path_variable(conda_bin_path)
    remote.remove_from_path_variable_containing("miniconda/bin")

    # Run the export command also in the current shell, so that the conda commands can be found
    remote.add_to_path_variable(conda_bin_path, comment=comment, in_shell=True)

    # Return paths
    return conda_installation_path, conda_main_executable_path

# -----------------------------------------------------------------

def has_valid_conda_environment_local(environment_name, conda_path="conda", conda_installation_path=None):

    """
    This function ...
    :param environment_name:
    :param conda_path:
    :param conda_installation_path:
    :return:
    """

    # Check
    if conda.is_environment(environment_name, conda_path):

        # Check if the bin dir is present
        env_path = fs.join(conda_installation_path, "envs", environment_name)
        bin_path = fs.join(env_path, "bin")
        python_path = fs.join(bin_path, "python")

        if not fs.is_directory(bin_path):

            # Remove the environment
            fs.remove_directory(env_path)
            return False

        elif fs.is_empty(bin_path):

            # Remove the environment
            fs.remove_directory(env_path)
            return False

        elif not fs.is_file(python_path):

            # Remove the environment
            fs.remove_directory(env_path)
            return False

        else: return True

    else: return False

# -----------------------------------------------------------------

def has_valid_conda_environment_remote(remote, environment_name, conda_path="conda", conda_installation_path=None):

    """
    This function ...
    :param remote:
    :param environment_name:
    :param conda_path:
    :param conda_installation_path:
    :return:
    """

    # Is conda environment
    if remote.is_conda_environment(environment_name, conda_path):

        # Check if the bin dir is present
        env_path = fs.join(conda_installation_path, "envs", environment_name)
        bin_path = fs.join(env_path, "bin")
        python_path = fs.join(bin_path, "python")
        if not remote.is_directory(bin_path):

            # Remove the environment
            remote.remove_directory(env_path)
            return False

        elif remote.is_empty(bin_path):

            # Remove the environment
            remote.remove_directory(env_path)
            return False

        elif not remote.is_file(python_path):

            # Remove the environment
            remote.remove_directory(env_path)
            return False

        else: return True

    # Not a conda environment
    else: return False

# -----------------------------------------------------------------

def create_conda_environment_local(environment_name, conda_installation_path, pts_root_path, python_version="2.7", conda_path="conda"):

    """
    This function ...
    :param environment_name:
    :param conda_installation_path:
    :param pts_root_path:
    :param python_version:
    :param conda_path:
    :return:
    """

    # Generate the command
    command = conda_path + " create -n " + environment_name + " python=" + python_version

    # Expect question
    expect = "Proceed ([y]/n)?"
    lines = []
    lines.append(command)
    lines.append((expect, "y", True))

    # Execute the lines
    terminal.execute_lines_expect_clone(*lines, show_output=log.is_debug)

    # Check whether the environment has been made and the executables are present
    environment_bin_path = fs.join(conda_installation_path, "envs", environment_name, "bin")
    if not fs.is_directory(environment_bin_path): raise RuntimeError("Creating the environment failed")
    conda_executable_path = fs.join(environment_bin_path, "conda")
    conda_pip_path = fs.join(environment_bin_path, "pip")
    conda_activate_path = fs.join(environment_bin_path, "activate")
    conda_python_path = fs.join(environment_bin_path, "python")
    conda_easy_install_path = fs.join(environment_bin_path, "easy_install")

    # Clear previous things
    #comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
    #terminal.remove_aliases(environment_name, "pts", "ipts")
    #terminal.remove_aliases_and_variables_with_comment(comment)
    #if pts_root_path in terminal.paths_in_python_path_variable(): terminal.remove_from_python_path_variable(pts_root_path)

    # Add an alias for the PTS python version
    #terminal.define_alias(environment_name, conda_python_path, comment=comment, in_shell=True)

    # Add PTS to shell configuration file
    #terminal.add_to_python_path_variable(pts_root_path, comment=comment, in_shell=True)
    #terminal.define_alias("pts", conda_python_path + " -m pts.do", comment=comment, in_shell=True)
    #terminal.define_alias("ipts", conda_python_path + " -im pts.do", comment=comment, in_shell=True)

    # Return the paths
    return conda_executable_path, conda_pip_path, conda_activate_path, conda_python_path, conda_easy_install_path

# -----------------------------------------------------------------

def setup_conda_environment_local(environment_name, pip_name, jupyter_name, pts_root_path, python_path, pip_path, jupyter_path):

    """
    This function ...
    :param environment_name:
    :param pip_name:
    :param jupyter_name:
    :param pts_root_path:
    :param python_path:
    :param pip_path:
    :param jupyter_path:
    :return:
    """

    # Clear previous things
    comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
    terminal.remove_aliases(environment_name, pip_name, "pts", "ipts")
    terminal.remove_aliases_and_variables_with_comment(comment)
    if pts_root_path in terminal.paths_in_python_path_variable(): terminal.remove_from_python_path_variable(pts_root_path)

    # Add an alias for the PTS python version
    terminal.define_alias(environment_name, python_path, comment=comment, in_shell=True)

    # Add an alias for pip
    terminal.define_alias(pip_name, pip_path, comment=comment, in_shell=True)

    # Add an alias for jupyter
    terminal.define_alias(jupyter_name, jupyter_path, comment=comment, in_shell=True)

    # Add PTS to shell configuration file
    terminal.add_to_python_path_variable(pts_root_path, comment=comment, in_shell=True)
    terminal.define_alias("pts", python_path + " -m pts.do", comment=comment, in_shell=True)
    terminal.define_alias("ipts", python_path + " -im pts.do", comment=comment, in_shell=True)

# -----------------------------------------------------------------

def create_conda_environment_remote(remote, environment_name, conda_installation_path, pts_root_path, python_version="2.7", conda_path="conda"):

    """
    This function ...
    :param remote:
    :param environment_name:
    :param conda_installation_path:
    :param pts_root_path:
    :param python_version:
    :param conda_path:
    :return:
    """

    # Generate the command
    command = conda_path + " create -n " + environment_name + " python=" + python_version

    # Expect question
    expect = "Proceed ([y]/n)?"
    lines = []
    lines.append(command)
    lines.append((expect, "y", True))

    # Execute the commands
    remote.execute_lines(*lines, show_output=log.is_debug)

    # Check whether the environment has been made and the executables are present
    environment_bin_path = fs.join(conda_installation_path, "envs", environment_name, "bin")
    if not remote.is_directory(environment_bin_path): raise RuntimeError("Creating the environment failed")
    conda_executable_path = fs.join(environment_bin_path, "conda")
    conda_pip_path = fs.join(environment_bin_path, "pip")
    conda_activate_path = fs.join(environment_bin_path, "activate")
    conda_python_path = fs.join(environment_bin_path, "python")
    conda_easy_install_path = fs.join(environment_bin_path, "easy_install")
    conda_jupyter_path = fs.join(environment_bin_path, "jupyter")

    # Clear previous things
    #comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
    #remote.remove_aliases(environment_name, "pts", "ipts")
    #remote.remove_aliases_and_variables_with_comment(comment)
    #if pts_root_path in remote.paths_in_python_path_variable: remote.remove_from_python_path_variable(pts_root_path)

    # Add an alias for the PTS python version
    #remote.define_alias(environment_name, conda_python_path, comment=comment, in_shell=True)

    # Add PTS to shell configuration file
    #remote.add_to_python_path_variable(pts_root_path, comment=comment, in_shell=True)
    #remote.define_alias("pts", conda_python_path + " -m pts.do", comment=comment, in_shell=True)
    #remote.define_alias("ipts", conda_python_path + " -im pts.do", comment=comment, in_shell=True)

    # Return the paths
    return conda_executable_path, conda_pip_path, conda_activate_path, conda_python_path, conda_easy_install_path, conda_jupyter_path

# -----------------------------------------------------------------

def setup_conda_environment_remote(remote, environment_name, pip_name, jupyter_name, pts_root_path, python_path, pip_path, jupyter_path=None):

    """
    This function ...
    :param remote:
    :param environment_name:
    :param pip_name:
    :param jupyter_name:
    :param pts_root_path:
    :param python_path:
    :param pip_path:
    :param jupyter_path:
    :return:
    """

    # Clear previous things
    comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
    remote.remove_aliases(environment_name, pip_name, "pts", "ipts")
    remote.remove_aliases_and_variables_with_comment(comment)
    if pts_root_path in remote.paths_in_python_path_variable: remote.remove_from_python_path_variable(pts_root_path)

    # Add an alias for the PTS python version
    remote.define_alias(environment_name, python_path, comment=comment, in_shell=True)

    # Add an alias for pip
    remote.define_alias(pip_name, pip_path, comment=comment, in_shell=True)

    # Add an alias for jupyter
    if jupyter_path is not None: terminal.define_alias(jupyter_name, jupyter_path, comment=comment, in_shell=True)

    # Add PTS to shell configuration file
    remote.add_to_python_path_variable(pts_root_path, comment=comment, in_shell=True)
    remote.define_alias("pts", python_path + " -m pts.do", comment=comment, in_shell=True)
    remote.define_alias("ipts", python_path + " -im pts.do", comment=comment, in_shell=True)

# -----------------------------------------------------------------

def can_import_module(module_name, remote=None, python_path="python"):

    """
    This function ...
    :param module_name:
    :param remote:
    :param python_path:
    :return:
    """

    if remote is not None: output = remote.execute(python_path + " -c 'import " + module_name + "'")
    else: output = terminal.execute(python_path + " -c 'import " + module_name + "'")

    # Check the output
    if len(output) == 0: return True
    else:
        # Look for Import error
        for line in output:
            if "ImportError:" in line: return False
        else: return True

# -----------------------------------------------------------------

def install_module_remote(remote, command, conda_path, module, conda_environment, check_present=True):

    """
    This function ...
    :param remote:
    :param command:
    :param conda_path:
    :param module:
    :param conda_environment:
    :param check_present:
    :return:
    """

    # Launch the installation command
    try:

        # Multiple commands
        if isinstance(command, list):

            # Skip module if already installed together with another package in the meantime
            if command[0].startswith(conda_path) and check_present:
                if remote.is_present_conda_package(module, conda_environment, conda_path): return False
            output = remote.execute_lines(*command, show_output=log.is_debug)

        # Simple command
        elif isinstance(command, basestring):

            if "setup.py" in command:

                # special: for python setup.py, we must be in the directory or it won't work
                dir_path = fs.directory_of(command.split()[1])
                setup_path = fs.join(dir_path, "setup.py")
                command.replace(setup_path, "setup.py")
                output = remote.execute(command, show_output=log.is_debug, cwd=dir_path)

            # Execute the command
            else: output = remote.execute(command, show_output=log.is_debug)

        else: raise ValueError("Invalid command: " + str(command))

        # Check the output
        for line in output:
            if "Exception:" in line:
                failed = True
                break
        else: failed = False

        if failed: return output
        else: return True

    # An error occured
    except Exception, err:

        tb = traceback.extract_tb(sys.exc_traceback)
        lines = traceback.format_list(tb)
        return lines

# -----------------------------------------------------------------
