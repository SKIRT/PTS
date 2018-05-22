#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.versionchecker Contains the VersionChecker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.basics.log import log
from .configurable import RemotesConfigurable
from pts.core.tools import conda
from .modules import Modules

# -----------------------------------------------------------------

class VersionChecker(RemotesConfigurable):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(VersionChecker, self).__init__(*args, **kwargs)

        # Versions for each remote host
        self.os_versions = dict()
        self.python_versions = dict()
        self.cpp_versions = dict()
        self.mpi_versions = dict()
        self.qt_versions = dict()
        self.skirt_versions = dict()
        self.pts_versions = dict()
        self.conda_versions = dict()
        self.skirt_issues = dict()
        self.pts_issues = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Check versions
        self.check()

        # 3. Show versions if requested
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(VersionChecker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking on each remote host ...")

        # Loop over the remotes
        for remote in self.remotes:

            # Get host ID
            host_id = remote.host_id

            # Create the modules object
            modules = Modules(remote)

            # Get the operating system
            os_version = remote.operating_system_short

            # Get version of compiler, MPI, qmake, and python
            compiler_version = modules.versions["cpp"] if "cpp" in modules.versions else None
            mpi_compiler_version = modules.versions["mpi"] if "mpi" in modules.versions else None
            qmake_version = modules.versions["qmake"] if "qmake" in modules.versions else None
            python_version = modules.versions["python"] if "python" in modules.versions else None

            # Find conda
            conda_installation_path, conda_main_executable_path = remote.find_conda()

            # Check conda version
            if conda_main_executable_path is not None: conda_version = remote.conda_version_at(conda_main_executable_path)
            else: conda_version = None

            # Get SKIRT and PTS version
            skirt_version = remote.skirt_version
            pts_version = remote.pts_version

            # Get conformity issues
            pts_issues = remote.pts_conformity_issues
            skirt_issues = remote.skirt_conformity_issues

            # Set versions
            self.os_versions[host_id] = os_version
            if compiler_version is not None: self.cpp_versions[host_id] = compiler_version
            if mpi_compiler_version is not None: self.mpi_versions[host_id] = mpi_compiler_version
            if qmake_version is not None: self.qt_versions[host_id] = qmake_version
            if python_version is not None: self.python_versions[host_id] = python_version
            if conda_version is not None: self.conda_versions[host_id] = conda_version
            if skirt_version is not None: self.skirt_versions[host_id] = skirt_version
            if pts_version is not None: self.pts_versions[host_id] = pts_version
            if pts_issues is not None: self.pts_issues[host_id] = pts_issues
            if skirt_issues is not None: self.skirt_issues[host_id] = skirt_issues

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing versions ...")

        # OS
        print("")
        print(fmt.bold + fmt.blue + "Operating system" + fmt.reset + ":")
        print("")
        print(" - local: " + introspection.operating_system_short())
        for host_id in self.host_ids: print(" - " + host_id + ": " + self.os_versions[host_id])
        print("")

        # SKIRT
        print(fmt.bold + fmt.blue + "SKIRT" + fmt.reset + ":")
        print("")
        if introspection.has_skirt(): print(" - local: " + fmt.green + introspection.skirt_version() + fmt.reset)
        else: print(" - local: " + fmt.red + "not found" + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.skirt_versions: print(" - " + host_id + ": " + fmt.green + self.skirt_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # PTS
        print(fmt.bold + fmt.blue + "PTS" + fmt.reset + ":")
        print("")
        print(" - local: " + fmt.green + introspection.pts_version() + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.pts_versions: print(" - " + host_id + ": " + fmt.green + self.pts_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # C++
        print(fmt.bold + fmt.blue + "C++ compiler" + fmt.reset + ":")
        print("")
        if introspection.has_cpp_compiler(): print(" - local: " + fmt.green + introspection.cpp_compiler_version() + fmt.reset)
        else: print(" - local: " + fmt.red + "not found" + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.cpp_versions: print(" - " + host_id + ": " + fmt.green + self.cpp_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # MPI
        print(fmt.bold + fmt.blue + "MPI compiler" + fmt.reset + ":")
        print("")
        if introspection.has_mpi_compiler(): print(" - local: " + fmt.green + introspection.mpi_compiler_version() + fmt.reset)
        else: print(" - local: " + fmt.red + "not found" + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.mpi_versions: print(" - " + host_id + ": " + fmt.green + self.mpi_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # Qmake
        print(fmt.bold + fmt.blue + "Qmake" + fmt.reset + ":")
        print("")
        if introspection.has_qmake(): print(" - local: " + fmt.green + introspection.qmake_version() + fmt.reset)
        else: print(" - local: " + fmt.red + "not found" + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.qt_versions: print(" - " + host_id + ": " + fmt.green + self.qt_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # Python
        print(fmt.bold + fmt.blue + "Python" + fmt.reset + ":")
        print("")
        print(" - local: " + fmt.green + introspection.python_version_long() + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.python_versions: print(" - " + host_id + ": " + fmt.green + self.python_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # Conda
        print(fmt.bold + fmt.blue + "Conda" + fmt.reset + ":")
        print("")
        conda_path = conda.find_conda()
        conda_version = conda.get_conda_version(conda_path=conda_path)
        if introspection.has_conda(): print(" - local: " + fmt.green + conda_version + fmt.reset)
        else: print(" - local: " + fmt.red + "not found" + fmt.reset)
        for host_id in self.host_ids:
            if host_id in self.conda_versions: print(" - " + host_id + ": " + fmt.green + self.conda_versions[host_id] + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.red + "not found" + fmt.reset)
        print("")

        # Conform SKIRT installation
        print(fmt.bold + fmt.blue + "Conformity SKIRT installation" + fmt.reset + ":")
        print("")
        if introspection.has_skirt():
            if introspection.skirt_installation_is_conform(): print(" - local: " + fmt.green + "installation is conform" + fmt.reset)
            else: print(" - local: " + fmt.red + " installation not conform (" + introspection.skirt_conformity_issues() + ")" + fmt.reset)
        for host_id in self.host_ids:
            if host_id not in self.skirt_versions: continue
            if host_id in self.skirt_issues: print(" - " + host_id + ": " + fmt.red + "installation not conform (" + self.skirt_issues[host_id] + ")" + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.green + "installation is conform" + fmt.reset)
        print("")

        # Conform PTS installation
        print(fmt.bold + fmt.blue + "Conformity PTS installation" + fmt.reset + ":")
        print("")
        if introspection.pts_installation_is_conform(): print(" - local: " + fmt.green + "installation is conform" + fmt.reset)
        else: print(" - local: " + fmt.red + " installation not conform (" + introspection.pts_conformity_issues() + ")" + fmt.reset)
        for host_id in self.host_ids:
            if host_id not in self.pts_versions: continue
            if host_id in self.pts_issues: print(" - " + host_id + ": " + fmt.red + "installation not conform (" + self.pts_issues[host_id] + ")" + fmt.reset)
            else: print(" - " + host_id + ": " + fmt.green + "installation is conform" + fmt.reset)
        print("")

# -----------------------------------------------------------------
