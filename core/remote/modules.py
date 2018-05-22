#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.modules Contains the Modules class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.log import log

# -----------------------------------------------------------------

class Modules(object):

    """
    This class ...
    """

    def __init__(self, remote):

        """
        This function ...
        :param remote:
        """

        # The remote instance
        self.remote = remote

        # Contains the module names
        self.names = dict()

        # Contains the executable paths
        self.paths = dict()

        # Contains the versions
        self.versions = dict()

        # Check the modules
        self.check_compilers()
        self.check_qt()
        self.check_git()
        self.check_python()

    # -----------------------------------------------------------------

    @property
    def all_modules(self):

        """
        This function ...
        :return:
        """

        return list(set(self.names.values()))

    # -----------------------------------------------------------------

    @property
    def all_paths(self):

        """
        This function ...
        :return:
        """

        return list(set(self.paths.values()))

    # -----------------------------------------------------------------

    @property
    def has_cpp(self):

        """
        Thisf unction ...
        :return:
        """

        return "cpp" in self.paths and self.paths["cpp"] is not None

    # -----------------------------------------------------------------

    @property
    def has_mpi(self):

        """
        This function ...
        :return:
        """

        return "mpi" in self.paths and self.paths["mpi"] is not None

    # -----------------------------------------------------------------

    @property
    def has_qmake(self):

        """
        This function ...
        :return:
        """

        return "qmake" in self.paths and self.paths["qmake"] is not None

    # -----------------------------------------------------------------

    @property
    def has_git(self):

        """
        This function ...
        :return:
        """

        return "git" in self.paths and self.paths["git"] is not None

    # -----------------------------------------------------------------

    @property
    def has_python(self):

        """
        This function ...
        :return:
        """

        return "python" in self.paths and self.paths["python"] is not None

    # -----------------------------------------------------------------

    def check_compilers(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of C++ and MPI compilers ...")

        # Get the compiler paths
        compiler_path, cpp_module = self.remote.find_and_load_cpp_compiler(return_module=True, show_output=log.is_debug)
        mpi_compiler_path, mpi_module = self.remote.find_and_load_mpi_compiler(return_module=True, show_output=log.is_debug)

        # Compiler is found
        if compiler_path is not None:

            # Debugging
            log.debug("The C++ compiler path is '" + compiler_path + "'")

            # CPP
            self.names["cpp"] = cpp_module
            self.paths["cpp"] = compiler_path
            self.versions["cpp"] = self.remote.version_of(compiler_path)

        # Compiler is not found
        else: log.debug("C++ compiler is not found on the remote")

        # MPI is found
        if mpi_compiler_path is not None:

            # Debugging
            log.debug("The MPI compiler path is '" + mpi_compiler_path + "'")

            # MPI
            self.names["mpi"] = mpi_module
            self.paths["mpi"] = mpi_compiler_path
            self.versions["mpi"] = self.remote.version_of(mpi_compiler_path)

        # MPI is not found
        else: log.debug("MPI compiler is not found on the remote")

        # Unload all modules to avoid conflicts with the other modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def check_qt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking for Qt installation on remote ...")

        # Load Qt module, find the qmake path
        qmake_path, module = self.remote.find_and_load_qmake(return_module=True, show_output=log.is_debug)

        # Qmake is found
        if qmake_path is not None:

            # Get the version
            version = self.remote.version_of(qmake_path)

            # Debugging
            log.debug("The qmake path is '" + qmake_path + "'")
            log.debug("The qmake version is '" + version + "'")

            # Add module info
            self.names["qmake"] = module
            self.paths["qmake"] = qmake_path
            self.versions["qmake"] = version

        # Qmake is not found
        else: log.debug("Qt installation is not found on the remote")

        # Unload all modules to avoid conflicts with the other modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def check_git(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of git ...")

        # Find and load git
        path, version, module = self.remote.find_and_load_git(return_module=True, show_output=log.is_debug)

        # Git is found
        if path is not None:

            # Debugging
            log.debug("The path of the git installation is '" + path + "'")
            log.debug("The version of git is '" + version + "'")

            # Add module version
            self.names["git"] = module
            self.paths["git"] = path
            self.versions["git"] = version

        # Git is not found
        else: log.debug("Git installation is not found on the remote")

        # Unload all modules to avoid conflicts with the other modules
        self.remote.unload_all_modules()

    # -----------------------------------------------------------------

    def check_python(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the presence of Python ...")

        # Find and load python
        path, module = self.remote.find_and_load_python(return_module=True)

        # Python is found
        if path is not None:

            # Get the version
            version = self.remote.custom_python_version_long(path, show_output=True)
            if version is None: log.warning("Could not load the python module '" + module + "'")
            else:

                # Debugging
                log.debug("The path of the Python installation is '" + path + "'")
                log.debug("The version of Python is '" + version + "'")

                # Add module path and version
                self.names["python"] = module
                self.paths["python"] = path
                self.versions["python"] = version

        # Python is not found
        else: log.debug("Python installation is not found on the remote")

        # Unload all modules to avoid conflicts
        self.remote.unload_all_modules()

# -----------------------------------------------------------------
