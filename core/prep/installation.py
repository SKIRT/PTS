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
import urllib
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..simulation.execute import SkirtExec
from ..simulation.remote import Remote
from ..tools import introspection
from ..tools import filesystem as fs
from ..tools.logging import log

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

    @abstractmethod
    def create_directories_remote(self):

        """
        This fucntion ...
        :return:
        """

        pass

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

        # Install SKIRT
        self.install_skirt_local()

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

    # -----------------------------------------------------------------

    def install_skirt_local(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def install_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Installing SKIRT remotely ...")

        # Check Qt installation
        has_qt = self.check_qt_remote()

        # Install Qt
        if not has_qt: self.install_qt_remote()

        # Get the SKIRT code
        self.get_skirt_remote()

        # Install SKIRT
        self.install_skirt_remote()

    # -----------------------------------------------------------------

    def check_qt_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking for Qt installation on remote ...")

        # Installation modules are defined
        if self.remote.host.installation_modules is None:

            # Unload modules that are already loaded
            self.remote.unload_all_modules()

            # Load modules necessary for installation (this will hopefully bring out the qmake executable)
            self.remote.load_installation_modules()

        # Check if qmake can be found
        has_qt = self.remote.is_executable("qmake")

        # Return
        return has_qt

    # -----------------------------------------------------------------

    def install_qt_remote(self):

        """
        This function ...
        :return:
        """

        link = "http://download.qt.io/official_releases/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"

        # Determine the path for the Qt source code
        path = fs.join(None, "qt.tar.gz")

        # Create a
        urllib.urlretrieve(self.config.qt_link, path)

        # FOR HPC: module load Qt/5.2.1-intel-2015a

        # CHECK IF QMAKE CAN BE FOUND

    # -----------------------------------------------------------------

    def get_skirt_remote(self):

        """
        This function ...
        :return:
        """


    # -----------------------------------------------------------------

    def install_skirt_remote(self):

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

                    command = "wget " + miniconda_linux_url + " -O " + conda_installer_path
                    self.remote.execute(command, show_output=True)

                #conda_installer_path = fs.join(self.remote.home_directory, fs.name(miniconda_linux_url))

                # Run the installer
                #self.remote.execute("sh " + conda_installer_path, show_output=True)

                conda_installation_path = fs.join(self.remote.home_directory, "miniconda")

                if not self.remote.is_directory(conda_installation_path):

                    command = "bash " + conda_installer_path + " -b -p " + conda_installation_path
                    self.remote.execute(command, show_output=True)

                conda_bin_path = fs.join(conda_installation_path, "bin")
                conda_executable_path = fs.join(conda_bin_path, "conda")
                conda_pip_path = fs.join(conda_bin_path, "pip")
                conda_python_path = fs.join(conda_bin_path, "python")

                # Add conda bin path to bashrc
                bashrc_path = fs.join(self.remote.home_directory, ".bashrc")
                line = 'PATH=' + conda_bin_path + ':$PATH'
                self.remote.append_line(bashrc_path, line)
                self.remote.execute("source " + bashrc_path)

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
        lines.append("export PYTHONPATH=" + self.pts_root_path + ":$PYTHONPATH")
        lines.append('alias pts="python -m pts.do"')
        lines.append('alias ipts="python -im pts.do"')
        self.remote.append_lines(bashrc_path, lines)
        self.remote.execute("source " + bashrc_path)

    # -----------------------------------------------------------------

    def get_dependencies_remote(self):

        """
        This function ...
        :return:
        """

        # Execute PTS depends
        output = self.remote.execute("pts depends")

        print(output)

        self.remote.start_python_session()
        self.remote.import_python_package("introspection", from_name="pts.core.tools")

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
