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
from pts.core.tools.logging import log
from .configurable import RemotesConfigurable

# -----------------------------------------------------------------

class VersionChecker(RemotesConfigurable):
    
    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(VersionChecker, self).__init__(config)

        # Versions for each remote host
        self.python_versions = dict()
        self.cpp_versions = dict()
        self.mpi_versions = dict()
        self.qt_versions = dict()
        self.skirt_versions = dict()
        self.pts_versions = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

            # Loading compilers
            compiler_path = remote.find_and_load_cpp_compiler()
            mpi_compiler_path = remote.find_and_load_mpi_compiler()

            # print(compiler_path)
            # print(mpi_compiler_path)

            # Get C++ compiler and MPI compiler version
            compiler_version = remote.version_of(compiler_path) if compiler_path is not None else None
            mpi_compiler_version = remote.version_of(mpi_compiler_path) if mpi_compiler_path is not None else None

            # Load Qt module, find the qmake path
            qmake_path = remote.find_and_load_qmake()

            # print(qmake_path)

            # Get qmake version
            qmake_version = remote.version_of(qmake_path) if qmake_path is not None else None

            # Load python
            python_path = remote.find_and_load_python()

            # Get python version
            python_version = remote.version_of(python_path)

            # Get SKIRT and PTS version
            skirt_version = remote.skirt_version
            pts_version = remote.pts_version

            # Set versions
            if compiler_version is not None: self.cpp_versions[host_id] = compiler_version
            if mpi_compiler_version is not None: self.mpi_versions[host_id] = mpi_compiler_version
            if qmake_version is not None: self.qt_versions[host_id] = qmake_version
            if python_version is not None: self.python_versions[host_id] = python_version
            if skirt_version is not None: self.skirt_versions[host_id] = skirt_version
            if pts_version is not None: self.pts_versions[host_id] = pts_version

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print("")
        print(fmt.bold + fmt.green + "SKIRT" + fmt.reset + ":")
        print("")
        print(" - local: " + introspection.skirt_version())
        for host_id in self.host_ids:
            if host_id in self.skirt_versions:
                print(" - " + host_id + ": " + self.skirt_versions[host_id])
            else: print(" - " + host_id + ": not found")
        print("")

        print(fmt.bold + fmt.green + "PTS" + fmt.reset + ":")
        print("")
        print(" - local: " + introspection.pts_version())
        for host_id in self.host_ids:
            if host_id in self.pts_versions:
                print(" - " + host_id + ": " + self.pts_versions[host_id])
            else:
                print(" - " + host_id + ": not found")
        print("")

        print(fmt.bold + fmt.green + "C++ compiler" + fmt.reset + ":")
        print("")
        print(" - local: " + " ...")
        for host_id in self.host_ids:
            if host_id in self.cpp_versions:
                print(" - " + host_id + ": " + self.cpp_versions[host_id])
            else:
                print(" - " + host_id + ": not found")
        print("")

        print(fmt.bold + fmt.green + "MPI compiler" + fmt.reset + ":")
        print("")
        print(" - local: " + "...")
        for host_id in self.host_ids:
            if host_id in self.mpi_versions:
                print(" - " + host_id + ": " + self.mpi_versions[host_id])
            else:
                print(" - " + host_id + ": not found")
        print("")

        print(fmt.bold + fmt.green + "Qmake" + fmt.reset + ":")
        print("")
        print(" - local: " + "...")
        for host_id in self.host_ids:
            if host_id in self.qt_versions:
                print(" - " + host_id + ": " + self.qt_versions[host_id])
            else:
                print(" - " + host_id + ": not found")
        print("")

        print(fmt.bold + fmt.green + "Python" + fmt.reset + ":")
        print("")
        print(" - local: " + introspection.python_version_long())
        for host_id in self.host_ids:
            if host_id in self.python_versions:
                print(" - " + host_id + ": " + self.python_versions[host_id])
            else:
                print(" - " + host_id + ": not found")
        print("")

# -----------------------------------------------------------------
