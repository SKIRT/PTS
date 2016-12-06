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
from ..simulation.execute import SkirtExec
from ..remote.remote import Remote
from ..tools import introspection
from ..tools import filesystem as fs
from ..tools.logging import log
from ..tools import google

# -----------------------------------------------------------------

# For SSH key:
# eval $(ssh-agent)
# ssh-add

# -----------------------------------------------------------------

skirt_directories = ["git", "run", "doc", "release", "debug"]
pts_directories = ["pts", "run", "doc", "temp", "remotes", "user", "ext"]

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

    def run(self):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2. Create the necessary directories
        self.create_directories()

        # 2. Install
        self.install()

        # 3. Test the installation
        self.test()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Installer, self).setup()

        # Setup the remote execution environment if necessary
        if self.config.remote is not None:

            # Create and setup the remote execution environment
            self.remote = Remote()
            self.remote.setup(self.config.remote)

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

# Private repository links
private_skirt_ssh_link = "git@github.ugent.be:SKIRT/SKIRT.git"
private_skirt_https_link = "https://github.ugent.be/SKIRT/SKIRT.git"

# Public repository links
public_skirt_https_link = "https://github.com/SKIRT/SKIRT.git"

# -----------------------------------------------------------------

# Determine Qt configure options
qt_configure_options = []
qt_configure_options.append("-prefix '$HOME/Qt/Desktop/5.2.1'")
qt_configure_options.append("-opensource")
qt_configure_options.append("-confirm-license")
qt_configure_options.append("-c++11")
qt_configure_options.append("-no-javascript-jit")
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
qt_configure_options.append("-no-nis")
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

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_directories_remote(self):

        """
        This function ...
        :return:
        """

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

        # Check if Qt is installed
        self.check_qt_local()

        # Install Qt
        self.install_qt_local()

        # Get the SKIRT code
        self.get_skirt_local()

        # Build SKIRT
        self.build_skirt_local()

    # -----------------------------------------------------------------

    def check_qt_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def install_qt_local(self):

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

        # Check whether the Qt version is supported
        # if [[ $VERSION > '5.2.0' ]]

    # -----------------------------------------------------------------

    def get_skirt_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def build_skirt_local(self):

        """
        This function ...
        :return:
        """

        pass

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

        # Check Qt installation
        self.check_qt_remote()

        # Install Qt if necessary
        if not self.has_qt: self.install_qt_remote()

        # Get the SKIRT code
        self.get_skirt_remote()

        # Build SKIRT
        self.build_skirt_remote()

    # -----------------------------------------------------------------

    def check_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking for Qt installation on remote ...")

        # Installation modules are defined
        if self.remote.host.installation_modules is not None:

            # Unload modules that are already loaded
            self.remote.unload_all_modules()

            # Load modules necessary for installation (this will hopefully bring out the qmake executable)
            self.remote.load_installation_modules()

        # Keep the qmake paths in a list to decide later which one we can use
        qmake_paths = []

        # Search for qmake in the home directory
        command = "find " + self.remote.home_directory + "/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.remote.execute(command)

        # Search for qmake in the /usr/local directory
        command = "find /usr/local/Qt* -name qmake -type f 2>/dev/null"
        qmake_paths += self.remote.execute(command)

        ## TOO SLOW?
        # Search for Qt directories in the home directory
        #for directory_path in self.remote.directories_in_path(self.remote.home_directory, startswith="Qt"):
            # Search for 'qmake' executables
            #for path in self.remote.files_in_path(directory_path, recursive=True):
                # Add the path
                #qmake_paths.append(path)

        ## TOO SLOW?
        # Search for Qt directories in /usr/local
        #for directory_path in self.remote.directories_in_path("/usr/local", startswith="Qt"):
            # Search for 'qmake' executables
            #for path in self.remote.files_in_path(directory_path, recursive=True):
                # Add the path
                #qmake_paths.append(path)

        # Check if qmake can be found by running 'which'
        qmake_path = self.remote.find_executable("qmake")
        if qmake_path is not None: qmake_paths.append(qmake_path)

        latest_version = None

        # Loop over the qmake paths
        for qmake_path in qmake_paths:

            # Get the version
            output = self.remote.execute(qmake_path + " -v")

            qt_version = output[1].split("Qt version ")[1].split(" in")[0]

            if qt_version < "5.2.0": continue # oldest supported version
            if "conda" in qmake_path: continue
            if "canopy" in qmake_path: continue
            if "epd" in qmake_path: continue
            if "enthought" in qmake_path: continue

            if latest_version is None or qt_version > latest_version:

                latest_version = qt_version
                self.qmake_path = qmake_path

    # -----------------------------------------------------------------

    def install_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing Qt ...")

        #link = "http://download.qt.io/official_releases/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"

        # Qt URL
        url = "http://download.qt.io/official_releases/qt/5.7/5.7.0/single/qt-everywhere-opensource-src-5.7.0.tar.gz"

        # Determine the path for the Qt source code
        path = fs.join(self.remote.home_directory, "qt.tar.gz")

        # Download Qt
        self.remote.download_from_url_to(url, path)

        # Determine commands
        configure_command = "./configure " + " ".join(qt_configure_options)
        make_command = "make"
        install_command = "make install"

        # Execute the commands
        self.remote.execute_lines(configure_command, make_command, install_command, show_output=True)

        # FOR HPC: module load Qt/5.2.1-intel-2015a

        # CHECK IF QMAKE CAN BE FOUND

    # -----------------------------------------------------------------

    def get_skirt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the SKIRT source code ...")

        # Get repository link
        if self.config.repository is not None:
            url = introspection.skirt_git_remote_url(self.config.repository)
        elif self.config.private:
            url = private_skirt_https_link
        else: url = public_skirt_https_link

        # CONVERT TO HTTPS LINK
        host = url.split("@")[1].split(":")[0]
        user_or_organization = url.split(":")[1].split("/")[0]
        repo_name = url.split("/")[-1].split(".git")[0]
        url = "https://" + host + "/" + user_or_organization + "/" + repo_name + ".git"

        # Set the clone command
        command = "git clone " + url + " " + self.skirt_repo_path

        # Find the account file for the repository host (e.g. github.ugent.be)
        username, password = introspection.get_account(host)

        # Set the command lines
        lines = []
        lines.append(command)
        lines.append(("':", username))
        lines.append(("':", password))

        # Clone the repository
        self.remote.execute_lines(*lines, show_output=True)

        # Set PYTHONPATH
        bashrc_path = fs.join(self.remote.home_directory, ".bashrc")
        lines = []
        export_command = "export PATH=" + fs.join(self.skirt_release_path, "SKIRTmain") + ":" + fs.join(self.skirt_release_path, "FitSKIRTmain") + ":$PATH"
        lines.append("")
        lines.append("# For SKIRT and FitSKIRT, added by PTS (Python Toolkit for SKIRT)")
        lines.append(export_command)
        lines.append("")
        self.remote.append_lines(bashrc_path, lines)

        # Set the path to the main SKIRT executable
        self.skirt_path = fs.join(self.skirt_release_path, "SKIRTmain", "skirt")

        # Run export path in the current shell to make SKIRT command visible
        self.remote.execute(export_command)

        # Load bashrc file
        #self.remote.execute("source " + bashrc_path)

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

        #skirt_make_script = "makeSKIRT.sh"

        # Navigate to the SKIRT repo directory
        #self.remote.change_cwd(self.skirt_repo_path)

        # Execute the script
        #self.remote.execute("./" + skirt_make_script, show_output=True)

        # # Create the make file and perform the build
        # $QMAKEPATH BuildSKIRT.pro -o ../release/Makefile CONFIG+=release
        # make -j ${1:-5} -w -C ../release

        # Navigate to the SKIRT repo directory
        self.remote.change_cwd(self.skirt_repo_path)

        # Create command strings
        make_make_command = self.qmake_path + " BuildSKIRT.pro -o ../release/Makefile CONFIG+=release"
        nthreads = self.remote.cores_per_socket
        make_command = "make -j " + str(nthreads) + " -w -C ../release"

        # Execute the commands
        self.remote.execute(make_make_command, show_output=True)
        self.remote.execute(make_command, show_output=True)

        # Success
        log.success("SKIRT was successfully built")

    # -----------------------------------------------------------------

    def build_skirt_hpc(self):

        """
        This function ...
        :return:
        """

        local_script_path = None

        screen_name = "SKIRT installation"

        # Open the job script file
        script_file = open(local_script_path, 'w')

        # Write a general header to the batch script
        script_file.write("#!/bin/sh\n")
        script_file.write("# Batch script for running SKIRT on a remote system\n")
        script_file.write("# To execute manualy, copy this file to the remote filesystem and enter the following commmand:\n")
        script_file.write("# screen -S " + screen_name + " -L -d -m " + fs.name(local_script_path) + "'\n")
        script_file.write("\n")

        # Load modules
        script_file.write("module load lxml/3.4.2-intel-2015a-Python-2.7.9")
        script_file.write("module load Qt/5.2.1-intel-2015a")

        #
        script_file.write("./makeSKIRT.sh")

        self.remote.start_screen(name, local_script_path, script_destination, screen_output_path=None, keep_remote_script=False)

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

# Private repository links
private_pts_ssh_link = "git@github.ugent.be:SKIRT/PTS.git"
private_pts_https_link = "https://github.ugent.be/SKIRT/PTS.git"

# Public repository links
public_pts_link = "https://github.com/SKIRT/PTS.git"

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

        self.conda_executable_path = None
        self.conda_pip_path = None
        self.conda_python_path = None

        self.pts_path = None

    # -----------------------------------------------------------------

    def create_directories_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_directories_remote(self):

        """
        This function ...
        :return:
        """

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

        # Get a python distribution
        self.get_python_distribution_local()

        # Get PTS
        self.get_pts_local()

    # -----------------------------------------------------------------

    def get_python_distribution_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_pts_local(self):

        """
        This function ...
        :return:
        """

        pass

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

        if self.remote.in_python_virtual_environment():

            self.python_path = self.remote.execute("which python")[0]

        else:

            if self.remote.platform == "MacOS":

                # Add conda path to .profile
                profile_path = fs.join(self.remote.home_directory, ".profile")

                # ...

            elif self.remote.platform == "Linux":

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
                    self.remote.execute(command, show_output=True)

                conda_bin_path = fs.join(conda_installation_path, "bin")
                self.conda_executable_path = fs.join(conda_bin_path, "conda")
                self.conda_pip_path = fs.join(conda_bin_path, "pip")
                self.conda_python_path = fs.join(conda_bin_path, "python")

                # Add conda bin path to bashrc
                bashrc_path = fs.join(self.remote.home_directory, ".bashrc")
                line = 'PATH=' + conda_bin_path + ':$PATH'
                lines = []
                lines.append("")
                lines.append("# For Miniconda, added by PTS (Python Toolkit for SKIRT)")
                lines.append(line)
                lines.append("")

                # Debugging
                log.debug("Adding the conda executables to the PATH ...")
                self.remote.append_lines(bashrc_path, lines)

                # Run commands in current shell, so that the conda commands can be found
                self.remote.execute(line)

                # Debugging
                #log.debug("Sourcing the bashrc file ...")
                #self.remote.execute("source " + bashrc_path)

    # -----------------------------------------------------------------

    def get_pts_remote(self):

        """
        This function ...
        :return:
        """

        if self.config.repository is not None:
            url = introspection.pts_git_remote_url(self.config.repository)
        elif self.config.private:
            url = private_pts_https_link
        else: url = public_pts_link


        # CONVERT TO HTTPS LINK
        # git@github.ugent.be:sjversto/PTS.git
        # to
        # https://github.ugent.be/SKIRT/PTS.git

        host = url.split("@")[1].split(":")[0]
        user_or_organization = url.split(":")[1].split("/")[0]
        repo_name = url.split("/")[-1].split(".git")[0]

        url = "https://" + host + "/" + user_or_organization + "/" + repo_name + ".git"

        # Set the clone command
        command = "git clone " + url + " " + self.pts_package_path

        # Find the account file for the repository host (e.g. github.ugent.be)
        username, password = introspection.get_account(host)

        # Set the command lines
        lines = []
        lines.append(command)
        lines.append(("':", username))
        lines.append(("':", password))

        # Clone the repository
        self.remote.execute_lines(*lines, show_output=True)

        # Set PYTHONPATH
        bashrc_path = fs.join(self.remote.home_directory, ".bashrc")
        lines = []
        lines.append("")
        lines.append("# For PTS, added by PTS (Python Toolkit for SKIRT)")
        lines.append("export PYTHONPATH=" + self.pts_root_path + ":$PYTHONPATH")
        lines.append('alias pts="python -m pts.do"')
        lines.append('alias ipts="python -im pts.do"')
        lines.append("")
        self.remote.append_lines(bashrc_path, lines)

        # Run commands in current shell, so that the pts command can be found
        self.remote.execute("export PYTHONPATH=" + self.pts_root_path + ":$PYTHONPATH")
        self.remote.execute('alias pts="python -m pts.do"')
        self.remote.execute('alias ipts="python -im pts.do"')

        # Load bashrc file
        # self.remote.execute("source " + bashrc_path)

        # Set the path to the main PTS executable
        self.pts_path = fs.join(self.pts_package_path, "do", "__main__.py")

    # -----------------------------------------------------------------

    def get_dependencies_remote(self):

        """
        This function ...
        :return:
        """

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
        self.remote.start_python_session()
        self.remote.import_python_package("introspection", from_name="pts.core.tools")
        dependencies = self.remote.get_simple_python_property("introspection", "get_all_dependencies().keys()")
        packages = self.remote.get_simple_python_property("introspection", "installed_python_packages()")
        #self.remote.end_python_session()
        # Don't end the python session just yet

        # Get installation commands
        installation_commands, installed, not_installed = get_installation_commands(dependencies, packages, already_installed, available_packages, self.remote)
        self.remote.end_python_session()

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

    # -----------------------------------------------------------------

    def test_local(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def test_remote(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def find_real_name(module_name, available_packages, remote):

    """
    This function ...
    :param module_name:
    :param available_packages:
    :param remote:
    :return:
    """

    if module_name in available_packages: return module_name

    # Look for real module name
    try:
        module_url = google.lucky(module_name)
    except Exception:
        if remote is not None:
            # use google on the remote end, because there are strange errors when using it on the client end when it has
            # a VPN connection open
            module_url = remote.get_simple_python_property("google", "lucky('" + module_name + "')")
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

def get_installation_commands(dependencies, packages, already_installed, available_packages, remote):

    """
    This function ...
    :return:
    """

    installed = []
    not_installed = []

    commands = dict()

    # Import google, don't pass the remote is importing the google module failed on the remote
    success = remote.import_python_package("google", from_name="pts.core.tools")
    if not success: remote = None

    # Loop over the dependencies
    for module in dependencies:

        module_name = module

        if module_name in packages: continue

        # Skip packages from the standard library
        if introspection.is_std_lib(module_name): continue

        # Check if already installed
        if module_name in already_installed: continue

        # Find name, check if available
        module_name, via = find_real_name(module_name, available_packages, remote)

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
