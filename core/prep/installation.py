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
import requests
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..remote.remote import Remote
from ..tools import introspection
from ..tools import filesystem as fs
from ..tools.logging import log
from ..tools import google
from ..tools import network, archive
from ..tools import terminal
from ..tools import git
from ..tools import parallelization
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from ..remote.modules import Modules

# -----------------------------------------------------------------

# For SSH key:
# eval $(ssh-agent)
# ssh-add

# -----------------------------------------------------------------

skirt_directories = ["git", "run", "doc", "release", "debug"]
pts_directories = ["pts", "run", "doc", "temp", "remotes", "user", "ext"]

# -----------------------------------------------f------------------

# Qt URL
# link = "http://download.qt.io/official_releases/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"
#url = "http://download.qt.io/official_releases/qt/5.7/5.7.0/single/qt-everywhere-opensource-src-5.7.0.tar.gz"
qt_url = "http://download.qt.io/official_releases/qt/5.7/5.7.1/single/qt-everywhere-opensource-src-5.7.1.tar.gz" # latest

# -----------------------------------------------------------------

class Installer(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(Installer, self).__init__(config)

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

        # Setup the remote execution environment if necessary
        elif self.config.host_id is not None:

            # Create and setup the remote execution environment
            self.remote = Remote()
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
qt_configure_options.append("-prefix '$HOME/Qt'")
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

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SKIRTInstaller, self).__init__(config)

        # The paths to the C++ compiler and MPI compiler
        self.compiler_path = None
        self.mpi_compiler_path = None

        # The path to the qmake executable corresponding to the most recent Qt installation
        self.qmake_path = None

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

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the directory structure ...")

        # Set paths
        self.skirt_root_path = fs.join(fs.home(), "SKIRT")
        self.skirt_repo_path = fs.join(self.skirt_root_path, "git")
        self.skirt_release_path = fs.join(self.skirt_root_path, "release")

        # Check if already present
        if fs.is_directory(self.skirt_root_path):
            if self.config.force: fs.remove_directory(self.skirt_root_path)
            else: raise RuntimeError("SKIRT is already installed (or partly present)")

        # Make the root directory
        fs.create_directory(self.skirt_root_path)

        # Create the other directories
        for name in skirt_directories:

            # Determine path
            path = fs.join(self.skirt_root_path, name)

            # Create the directory
            fs.create_directory(path)

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
            if self.config.force: self.remote.remove_directory(self.skirt_root_path)
            else: raise RuntimeError("SKIRT is already installed (or partly present) on the remote host")

        # Make the root directory
        self.remote.create_directory(self.skirt_root_path)

        # Create the other directories
        for name in skirt_directories:

            # Determine path
            path = fs.join(self.skirt_root_path, name)

            # Create the directory
            self.remote.create_directory(path)

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
        self.get_skirt_local()

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

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        # Find qmake path
        self.qmake_path = find_qmake()

    # -----------------------------------------------------------------

    def install_qt_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing Qt ...")

        # Determine the path for the Qt source code
        path = fs.join(fs.home(), "qt.tar.gz")

        # Download tar.gz file
        network.download_file(qt_url, path)

        # Decompress tar.gz file
        decompress_path = fs.join(fs.home(), "Qt-install")
        archive.decompress_file(path, decompress_path)

        # Determine commands
        configure_command = "./configure " + " ".join(qt_configure_options)
        make_command = "make"
        install_command = "make install"

        # Configure
        terminal.execute(configure_command, cwd=decompress_path)

        # Make
        terminal.execute(make_command, cwd=decompress_path)

        # Install
        terminal.execute(install_command, cwd=decompress_path)

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

            # Clone
            #child = pexpect.spawn(command, timeout=30)
            #child.expect([':'])
            #child.sendline(username)
            #child.expect([':'])
            #child.sendline(password)

            # Set the command lines
            lines = []
            lines.append(command)
            lines.append(("':", username))
            lines.append(("':", password))

            # Clone the repository
            terminal.execute_lines(*lines, cwd=self.skirt_root_path)

        else: terminal.execute(command, cwd=self.skirt_root_path)

        # Get git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Determine path of SKIRT and FitSKIRT main directories with executable
        comment = "For SKIRT and FitSKIRT, added by PTS (Python Toolkit for SKIRT)"
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")

        # Add to path, also execute in current shell to make SKIRT (and FtiSKIRT) visible
        add_to_environment_variable("PATH", skirt_main_path, comment=comment, in_shell=True)
        add_to_environment_variable("PATH", fitskirt_main_path, comment=comment, in_shell=True)

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

    def install_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing SKIRT remotely ...")

        # Check the compilers (C++ and MPI)
        self.check_compilers_remote()

        # Check Qt installation
        self.check_qt_remote()

        # Install Qt if necessary
        if not self.has_qt: self.install_qt_remote()

        # Check the presence of git:
        # no, loading the latest git version interferes with the intel compiler version on HPC UGent
        self.check_git_remote()

        # Get the SKIRT code
        self.get_skirt_remote()

        # Build SKIRT
        self.build_skirt_remote()

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

        # Debugging
        log.debug("The C++ compiler path is '" + self.compiler_path)
        log.debug("The MPI compiler path is '" + self.mpi_compiler_path)

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

        # Debugging
        log.debug("The qmake path is '" + self.qmake_path + "'")

    # -----------------------------------------------------------------

    def install_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing Qt ...")

        # Detemrine location
        if self.remote.host.scratch_path is not None: temp_path = self.remote.absolute_path(self.remote.host.scratch_path)
        else: temp_path = self.remote.home_directory

        # Determine the path for the Qt source code
        path = fs.join(temp_path, "qt.tar.gz")

        # Download Qt
        self.remote.download_from_url_to(qt_url, path)

        # Unarchive
        decompress_path = self.remote.create_directory_in(temp_path, "Qt-install")
        self.remote.decompress_file(path, decompress_path)

        decompress_path = fs.join(temp_path, "Qt-install")

        # Get the only directory in the Qt-install directory
        qt_everywhere_opensource_path = self.remote.directories_in_path(decompress_path)[0]

        # Determine commands
        configure_command = "./configure " + " ".join(qt_configure_options)
        make_command = "make"
        install_command = "make install"

        # Execute the commands
        self.remote.execute_lines(configure_command, make_command, install_command, show_output=log.is_debug(), cwd=qt_everywhere_opensource_path)

        # Remove decompressed folder and the tar.gz file
        self.remote.remove_file(path)
        self.remote.remove_directory(decompress_path)

        # Success
        log.success("Qt was succesfully installed")

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

        # Unload all modules to avoid conflicts with the other modules
        #self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def get_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the SKIRT source code ...")

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
        self.remote.execute(command, show_output=log.is_debug())

        # Get the git version
        self.git_version = git.get_short_git_version(self.skirt_repo_path, self.remote)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Determine SKIRT and FitSKIRT main paths
        skirt_main_path = fs.join(self.skirt_release_path, "SKIRTmain")
        fitskirt_main_path = fs.join(self.skirt_release_path, "FitSKIRTmain")

        # Add
        comment = "Added by the Python Toolkit for SKIRT (PTS)"
        self.remote.add_to_environment_variable("PATH", skirt_main_path, comment=comment, in_shell=True)
        self.remote.add_to_environment_variable("PATH", fitskirt_main_path, comment=comment, in_shell=True)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(skirt_main_path, "skirt")

        # Success
        log.success("SKIRT was successfully downloaded")

    # -----------------------------------------------------------------

    def build_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building SKIRT ...")

        # Execute the build
        build_skirt_on_remote(self.remote, self.skirt_repo_path, self.qmake_path, self.git_version)

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
                log.succes("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the SKIRT installation ...")

        output = self.remote.execute("skirt -h")
        for line in output:
            if "Welcome to SKIRT" in line:
                log.success("SKIRT is working")
                break
        else:
            log.error("Something is wrong with the SKIRT installation:")
            for line in output: log.error("   " + line)

# -----------------------------------------------------------------

def build_skirt_local(skirt_repo_path, qmake_path, git_version):

    """
    This function ...
    :param skirt_repo_path:
    :param qmake_path:
    :param git_version:
    :return:
    """

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
    output = terminal.execute(make_make_command, show_output=log.is_debug(), cwd=skirt_repo_path)

    # Overwrite the git version
    git_version_content = 'const char* git_version = " ' + git_version + ' " ;'
    git_version_path = fs.join(skirt_repo_path, "SKIRTmain", "git_version.h")
    fs.write_line(git_version_path, git_version_content)

    # Make
    terminal.execute(make_command, cwd=skirt_repo_path)

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
    output = remote.execute(make_make_command, show_output=log.is_debug(), cwd=skirt_repo_path)

    # Overwrite the git version
    git_version_content = 'const char* git_version = " ' + git_version + ' " ;'
    git_version_path = fs.join(skirt_repo_path, "SKIRTmain", "git_version.h")
    write_command = 'echo "' + git_version_content + '" > ' + git_version_path
    remote.execute(write_command)

    # Make
    output = remote.execute(make_command, show_output=log.is_debug(), cwd=skirt_repo_path)

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

class PTSInstaller(Installer):

    """
    This function ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(PTSInstaller, self).__init__(config)

        # Path to python executable
        self.python_path = None

        # Path to PTS root directory
        self.pts_root_path = None

        # Path to PTS/pts
        self.pts_package_path = None

        # Conda
        self.conda_executable_path = None
        self.conda_pip_path = None
        self.conda_python_path = None

        # Path to the PTS executable
        self.pts_path = None

        # The git version
        self.git_version = None

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        # Give warning
        log.warning("You are asking to install PTS locally, but you are already using PTS at the moment")

        # Ask to proceed with installing python environment and dependencies
        definition = ConfigurationDefinition()
        definition.add_flag("proceed", "Proceed by installing a clean Python 2 environment for PTS and automatically installing all dependencies?", False)
        setter = InteractiveConfigurationSetter("Install PTS local", add_logging=False, add_cwd=False)
        config = setter.run(definition, prompt_optional=True)

        # Proceed or quit
        if config.proceed: log.info("Proceeding ...")
        else:
            log.info("Quitting ...")
            exit()

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

    def install_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing PTS locally ...")

        # 1. Get a python distribution
        self.get_python_distribution_local()

        # 2. Get PTS
        self.get_pts_local()

        # 3. Get PTS dependencies
        self.get_dependencies_local()

    # -----------------------------------------------------------------

    def get_python_distribution_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting a Anaconda Python 2 distribution locally ...")

        # Determine path to profile file
        profile_path = introspection.shell_configuration_path()

        # Determine installer path
        installer_path = fs.join(fs.home(), "conda.sh")

        # Download the installer
        if self.remote.is_macos: network.download_file(miniconda_macos_url, installer_path, overwrite=True)
        elif self.remote.is_linux: network.download_file(miniconda_linux_url, installer_path, overwrite=True)
        else: raise NotImplementedError("OS must be Linux or MacOS")

        # Determine the conda installation directory
        conda_installation_path = fs.join(self.remote.home_directory, "miniconda")

        # Check if not present
        if fs.is_directory(conda_installation_path): raise RuntimeError("Miniconda is already installed in '" + conda_installation_path + "'")

        # Run the installation script
        command = "bash " + conda_installer_path + " -b -p " + conda_installation_path
        self.remote.execute(command, show_output=log.is_debug())

        # CREATE A PYTHON ENVIRONMENT FOR PTS

        conda_bin_path = fs.join(conda_installation_path, "bin")
        self.conda_executable_path = fs.join(conda_bin_path, "conda")
        self.conda_pip_path = fs.join(conda_bin_path, "pip")
        self.conda_python_path = fs.join(conda_bin_path, "python")

        # Debugging
        log.debug("Adding the conda executables to the PATH ...")

        # Add conda bin path to bashrc / profile
        comment = "For Miniconda, added by PTS (Python Toolkit for SKIRT)"

        # Run the export command also in the current shell, so that the conda commands can be found
        self.remote.add_to_environment_variable("PATH", conda_bin_path, comment=comment, shell=True)

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
            terminal.execute_lines(*lines)

        else: terminal.execute(command)

        # Get the git version
        self.git_version = git.get_short_git_version(self.pts_package_path)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Add PTS to shell configuration file
        comment = "For PTS, added by PTS (Python Toolkit for SKIRT)"
        add_to_environment_variable("PYTHONPATH", self.pts_root_path, comment=comment, in_shell=True)
        define_alias("pts", "python -m pts.do", comment=comment, in_shell=True)
        define_alias("ipts", "python -im pts.do", comment=comment, in_shell=True)

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
        output = terminal.execute("conda search")
        available_packages = []
        for line in output:
            if not line.split(" ")[0]: continue
            available_packages.append(line.split(" ")[0])

        # Get already installed packages
        already_installed = []
        for line in terminal.execute("conda list"):
            if line.startswith("#"): continue
            already_installed.append(line.split(" ")[0])

        # Use the introspection module on the remote end to get the dependencies and installed python packages
        # session = self.remote.start_python_session()
        # session.import_package("introspection", from_name="pts.core.tools")

        dependencies = introspection.get_all_dependencies().keys()
        packages = introspection.installed_python_packages()

        # Get installation commands
        installation_commands, installed, not_installed = get_installation_commands(dependencies, packages, already_installed, available_packages)

        # Install
        for module in installation_commands:

            # Debugging
            log.debug("Installing '" + module + "' ...")

            command = installation_commands[module]

            if isinstance(command, list): terminal.execute_lines(*command, show_output=log.is_debug())
            elif isinstance(command, basestring): terminal.execute(command, show_output=log.is_debug())

        # Show installed packages
        log.info("Packages that were installed:")
        for module in installed: log.info(" - " + module)

        # Show not installed packages
        log.info("Packages that could not be installed:")
        for module in not_installed: log.info(" - " + module)

        # Show already present packages
        log.info("Packages that were already present:")
        for module in already_installed: log.info(" - " + module)

    # -----------------------------------------------------------------

    def install_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing PTS remotely ...")

        # Get a python distribution
        self.get_python_distribution_remote()

        # Get PTS
        self.get_pts_remote()

        # Get PTS dependencies
        self.get_dependencies_remote()

    # -----------------------------------------------------------------

    def get_python_distribution_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting a Python distribution on the remote host ...")

        #if self.remote.in_python_virtual_environment(): self.python_path = self.remote.execute("which python")[0]
        #else:

        # MacOS
        if self.remote.is_macos:

            # Add conda path to .profile
            profile_path = fs.join(self.remote.home_directory, ".profile")

            # ...

        # Linux
        elif self.remote.is_linux == "Linux":

            conda_installer_path = fs.join(self.remote.home_directory, "conda.sh")

            # Download anaconda
            #self.remote.download(miniconda_linux_url, conda_installer_path)

            if not self.remote.is_file(conda_installer_path):

                # Download the installer
                self.remote.download_from_url_to(miniconda_linux_url, conda_installer_path)

            #conda_installer_path = fs.join(self.remote.home_directory, fs.name(miniconda_linux_url))

            # Run the installer
            #self.remote.execute("sh " + conda_installer_path, show_output=True)

            conda_installation_path = fs.join(self.remote.home_directory, "miniconda")

            if not self.remote.is_directory(conda_installation_path):

                command = "bash " + conda_installer_path + " -b -p " + conda_installation_path
                self.remote.execute(command, show_output=log.is_debug())

            conda_bin_path = fs.join(conda_installation_path, "bin")
            self.conda_executable_path = fs.join(conda_bin_path, "conda")
            self.conda_pip_path = fs.join(conda_bin_path, "pip")
            self.conda_python_path = fs.join(conda_bin_path, "python")

            # Debugging
            log.debug("Adding the conda executables to the PATH ...")

            # Add conda bin path to bashrc / profile
            comment = "For Miniconda, added by PTS (Python Toolkit for SKIRT)"

            # Run the export command also in the current shell, so that the conda commands can be found
            self.remote.add_to_environment_variable("PATH", conda_bin_path, comment=comment, shell=True)

        # Other
        else: raise NotImplementedError("OS must be MacOS or Linux")

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

        # Do HPC UGent in a different way because it seems only SSH is permitted and not HTTPS (but we don't want SSH
        # because of the private/public key thingy, so use a trick
        # DON'T NEED THIS AFTER ALL! WE CAN JUST ADD THE USERNAME AND PASSWORD TO THE HTTPS LINK AND IT WORKS!
        #if self.remote.host.name == "login.hpc.ugent.be": self.git_version = get_pts_hpc(self.remote, url, self.pts_root_path, self.pts_package_path)
        #else:

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
        self.remote.execute(command, show_output=log.is_debug())

        # Get the git version
        self.git_version = git.get_short_git_version(self.pts_package_path, self.remote)

        # Show the git version
        log.info("The git version to be installed is '" + self.git_version + "'")

        # Run commands in current shell, so that the pts command can be found
        comment = "Added by the Python Toolkit for SKIRT (PTS)"
        self.remote.add_to_environment_variable("PYTHONPATH", self.pts_root_path, comment=comment, in_shell=True)
        self.remote.define_alias("pts", "python -m pts.do", comment=comment, in_shell=True)
        self.remote.define_alias("ipts", "python -im pts.do", comment=comment, in_shell=True)

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
        get_pts_dependencies_remote(self.remote)

        # Success
        log.success("Succesfully installed the dependencies on the remote host")

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation ...")

        output = terminal.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing the PTS installation ...")

        output = self.remote.execute("pts -h")
        for line in output:
            if "usage: pts" in line: break
        else:
            log.error("Something is wrong with the PTS installation:")
            for line in output: log.error("   " + line)

# -----------------------------------------------------------------

def find_real_name(module_name, available_packages, session=None):

    """
    This function ...
    :param module_name:
    :param available_packages:
    :param session:
    :return:
    """

    if module_name in available_packages: return module_name, None

    # Look for real module name
    try:
        module_url = google.lucky(module_name)
    except Exception:
        if session is not None:
            # use google on the remote end, because there are strange errors when using it on the client end when it has
            # a VPN connection open
            module_url = session.get_simple_property("google", "lucky('" + module_name + "')")
        else: return None, None

    # Search for github.com/ name
    session = requests.session()
    r = session.get(module_url)
    page_as_string = r.content

    if "github.com/" in page_as_string:

        module_name = page_as_string.split("github.com/")[1].split("/")[0]

        if module_name in available_packages: return module_name, None
        else: return module_name, "github.com"

    if "pip install" in page_as_string:

        module_name = page_as_string.split("pip install ")[1].split(" ")[0]

        if module_name in available_packages: return module_name, None
        else: return module_name, "pip"

    # Not found
    return None, None

# -----------------------------------------------------------------

def get_installation_commands(dependencies, packages, already_installed, available_packages, session=None):

    """
    This function ...
    :param dependencies:
    :param packages:
    :param already_installed:
    :param available_packages:
    :param session:
    :return:
    """

    installed = []
    not_installed = []

    commands = dict()

    # Loop over the dependencies
    for module in dependencies:

        module_name = module

        if module_name in packages: continue

        # Skip packages from the standard library
        if introspection.is_std_lib(module_name): continue

        # Check if already installed
        if module_name in already_installed: continue

        # Find name, check if available
        module_name, via = find_real_name(module_name, available_packages, session)

        if module_name is None:
            log.warning("Package '" + module + "' can not be installed")
            not_installed.append(module)
            continue

        #log.debug("Installing '" + module + "' ...")

        if via is None:

            command = "conda install " + module_name

            # self.remote.execute(command, show_output=True)

            lines = []
            lines.append(command)
            lines.append(("Proceed ([y]/n)?", "y"))

            #self.remote.execute_lines(*lines, show_output=True)

            commands[module] = lines

        elif via.startswith("pip"):

            command = via

            #self.remote.execute(command, show_output=True)

            commands[module] = command

        else: # not implemented yet

            not_installed.append(module)

        # Add to installed
        installed.append(module_name)

        # Return ...
        return commands, installed, not_installed

# -----------------------------------------------------------------

def find_qmake():

    """
    This function ...
    :return:
    """

    # Translated from ./makeSKIRT.sh

    # Get a list of qmake paths installed on this system
    qmake_paths = []

    for qt_dir in fs.directories_in_path(fs.home(), startswith="Qt"):
        qmake_paths = fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

    for qt_dir in fs.directories_in_path("/usr/local", startswith="Qt"):
        qmake_paths += fs.files_in_path(qt_dir, recursive=True, exact_name="qmake", extension="")

    qmake_path = introspection.qmake_path()
    qmake_paths += [qmake_path] if qmake_path is not None else []

    # Get the most recent installation

    latest_qmake_path = None

    # Loop over the qmake paths: SAME IMPLEMENTATION AS IN Remote class -> _check_qt_remote_no_lmod
    for qmake_path in qmake_paths:

        # Get the version
        output = terminal.execute(qmake_path + " -v")

        qt_version = output[1].split("Qt version ")[1].split(" in")[0]

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

def add_to_environment_variable(variable_name, value, comment=None, in_shell=False):

    """
    This function ...
    :param variable_name:
    :param value:
    :param comment:
    :param in_shell:
    :return:
    """

    # Determine command
    export_command = "export " + variable_name + "=" + value + ":$PATH"

    # Define lines
    lines = []
    lines.append("")
    if comment is not None: lines.append("# " + comment)
    lines.append(export_command)
    lines.append("")

    # Add lines
    fs.append_lines(introspection.shell_configuration_path(), lines)

    # Run export path in the current shell to make variable visible
    if in_shell: terminal.execute(export_command)

# -----------------------------------------------------------------

def define_alias(name, alias_to, comment=None, in_shell=False):

    """
    This function ...
    :param name:
    :param alias_to:
    :param comment:
    :param in_shell:
    :return:
    """

    # Generate the command
    alias_command = 'alias ' + name + '="' + alias_to + '"'

    # Define lines
    lines = []
    lines.append("")
    if comment is not None: lines.append("# " + comment)
    lines.append(alias_command)
    lines.append("")

    # Add lines
    fs.append_lines(introspection.shell_configuration_path(), lines)

    # Execute in shell
    if in_shell: terminal.execute(alias_command)

# -----------------------------------------------------------------

def get_pts_dependencies_remote(remote):

    """
    This fucntion ...
    :param remote:
    :return:
    """

    # Get available conda packages
    output = remote.execute("conda search")
    available_packages = []
    for line in output:
        if not line.split(" ")[0]: continue
        available_packages.append(line.split(" ")[0])

    # Get already installed packages
    already_installed = []
    for line in remote.execute("conda list"):
        if line.startswith("#"): continue
        already_installed.append(line.split(" ")[0])

    ## Install essential packages: numpy and astropy
    for essential in ["numpy", "astropy"]:

        # Skip installation if it is already present
        if essential in already_installed: continue

        # NUMPY
        command = "conda install " + essential
        lines = []
        lines.append(command)
        lines.append(("Proceed ([y]/n)?", "y"))
        remote.execute_lines(*lines, show_output=log.is_debug())

    # Use the introspection module on the remote end to get the dependencies and installed python packages
    session = remote.start_python_session()
    session.import_package("introspection", from_name="pts.core.tools")
    dependencies = session.get_simple_property("introspection", "get_all_dependencies().keys()")
    packages = session.get_simple_property("introspection", "installed_python_packages()")
    # self.remote.end_python_session()
    # Don't end the python session just yet

    # Get installation commands
    session.import_package("google", from_name="pts.core.tools")
    installation_commands, installed, not_installed = get_installation_commands(dependencies, packages,
                                                                                already_installed, available_packages,
                                                                                session)

    # Stop the python session
    del session

    # Install
    for module in installation_commands:

        # Debugging
        log.debug("Installing '" + module + "' ...")

        command = installation_commands[module]

        if isinstance(command, list): remote.execute_lines(*command, show_output=log.is_debug())
        elif isinstance(command, basestring): remote.execute(command, show_output=log.is_debug())

    # Show installed packages
    log.info("Packages that were installed:")
    for module in installed: log.info(" - " + module)

    # Show not installed packages
    log.info("Packages that could not be installed:")
    for module in not_installed: log.info(" - " + module)

    # Show already present packages
    log.info("Packages that were already present:")
    for module in already_installed: log.info(" - " + module)

# -----------------------------------------------------------------
