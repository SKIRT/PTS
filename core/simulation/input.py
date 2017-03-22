#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.input Contains the SimulationInput class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import introspection, time

# -----------------------------------------------------------------

class SimulationInput(object):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        self.paths = dict()

        # Add args
        for arg in args:

            if fs.is_file(arg): self.add_file(arg)
            elif fs.is_directory(arg): self.add_directory(arg)
            else: raise IOError("The file or directory '" + arg + "' does not exist")

        # Add kwargs
        for name in kwargs:

            path = kwargs[name]
            if not fs.is_file(path): raise IOError("The file '" + path + "' does not exist")
            self.add_file(path, name)

    # -----------------------------------------------------------------

    @classmethod
    def from_any(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

        if isinstance(argument, list): return cls(*argument)
        elif isinstance(argument, basestring): return cls(argument)
        elif isinstance(argument, cls): return argument
        else: raise ValueError("Invalid input specification: should be list, string or SimulationInput object")

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return str(self.paths)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for name in self.paths: yield name, self.paths[name]

    # -----------------------------------------------------------------

    def __contains__(self, name):

        """
        This funtion ...
        :param name:
        :return:
        """

        return name in self.paths

    # -----------------------------------------------------------------

    def __setitem__(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        self.paths[name] = path

    # -----------------------------------------------------------------

    def __getitem__(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.paths[name]

    # -----------------------------------------------------------------

    def add_file(self, path, name=None):

        """
        This function ...
        :param path:
        :param name:
        :return:
        """

        # Determine filename
        if name is None: name = fs.name(path)

        # Check
        if name in self.paths: raise ValueError("Already a file in the simulation input with the name '" + name + "'")

        # Add to dictionary
        self.paths[name] = path

    # -----------------------------------------------------------------

    def add_directory(self, path):

        """
        This function ...
        :return:
        """

        # Add each file in the directory
        for path, name in fs.files_in_path(path, returns=["path", "name"], extensions=True): self.add_file(path, name)

    # -----------------------------------------------------------------

    @property
    def is_scattered(self):

        """
        This function ...
        :return:
        """

        # Determine the single directory where the input is placed
        dir_path = None
        scattered_input = False

        for name in self.paths:

            path = self.paths[name]

            this_dir_path = fs.directory_of(path)

            if dir_path is None: dir_path = this_dir_path

            elif dir_path != this_dir_path:
                # raise RuntimeError("Cannot convert this SkirtArguments instance to a command: input files should be placed in the same directory!")
                scattered_input = True
                break

        return scattered_input

    # -----------------------------------------------------------------

    @property
    def matching_filenames(self):

        """
        This function ...
        :return:
        """

        for name in self.paths:
            if name != fs.name(self.paths[name]): return False

        return True

    # -----------------------------------------------------------------

    def to_single_directory(self):

        """
        This function ...
        :return:
        """

        # Not scattered and filenames match their paths
        if not self.is_scattered and self.matching_filenames:

            first = self.paths[self.paths.keys()[0]]
            dir_path = fs.directory_of(first)

            return dir_path

        # Scatted and/or not matching filenames
        else:

            # Create temporary directory
            temp_path = fs.create_directory_in(introspection.pts_temp_dir, time.unique_name("SKIRT_input"))

            for name in self.paths:

                path = self.paths[name]
                fs.copy_file(path, temp_path, new_name=name)

            dir_path = temp_path

            return dir_path

# -----------------------------------------------------------------
