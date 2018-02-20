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
from ..tools import types

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

        # Dictionary of filepaths
        self.paths = dict()

        # Relative directory names
        self.relative_directories = []

        # Add args
        for arg in args:

            if fs.is_file(arg): self.add_file(arg)
            elif fs.is_directory(arg): self.add_directory(arg)
            else: self.add_relative_directory(arg) # probably a relative directory (e.g. 'in')  #raise IOError("The file or directory '" + arg + "' does not exist")

        # Add kwargs
        for name in kwargs:

            path = kwargs[name]
            if not fs.is_file(path): raise IOError("The file '" + path + "' does not exist")
            self.add_file(path, name)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, paths):

        """
        This function ...
        :param paths:
        :return:
        """

        return cls(*paths)

    # -----------------------------------------------------------------

    @classmethod
    def from_any(cls, argument):

        """
        This function ...
        :param argument:
        :return:
        """

        if isinstance(argument, cls): return argument
        elif types.is_sequence(argument): return cls(*argument)
        elif types.is_string_type(argument): return cls(argument)
        elif types.is_dictionary(argument): return cls(**argument)
        else: raise ValueError("Invalid input specification: should be list, dictionary, string or SimulationInput object")

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, prefix=None):

        """
        This function ...
        :param path:
        :param prefix:
        :return:
        """

        return cls.from_paths(fs.files_in_path(path, startswith=prefix, not_extension="ski"))

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_directory(cls, path, remote, prefix=None):

        """
        This function ...
        :param path:
        :param remote:
        :param prefix:
        :return:
        """

        return cls.from_paths(remote.files_in_path(path, startswith=prefix, not_extension="ski"))

    # -----------------------------------------------------------------

    @classmethod
    def unchecked(cls, **kwargs):

        """
        This function ....
        :param args:
        :param kwargs:
        :return:
        """

        # Create the simulation input object
        input = cls()

        # Add files from kwargs
        for name in kwargs:
            path = kwargs[name]
            input.add_file(path, name)

        # Return the input object
        return input

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.paths.keys()

    # -----------------------------------------------------------------

    def get_filepath(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.paths[name]

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

    def add_relative_directory(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.relative_directories.append(name)

    # -----------------------------------------------------------------

    @property
    def is_scattered(self):

        """
        This function ...
        :return:
        """

        return not self.has_single_directory

    # -----------------------------------------------------------------

    @property
    def has_single_directory(self):

        """
        This function ...
        :return:
        """

        return self.single_directory_path is not None

    # -----------------------------------------------------------------

    @property
    def single_directory_path(self):

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

        #return scattered_input
        if scattered_input: return None
        else: return dir_path

    # -----------------------------------------------------------------

    @property
    def directory_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()

        for name in self.paths:
            path = self.paths[name]
            dir_path = fs.directory_of(path)
            paths[name] = dir_path

        return paths

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

        # If just a relative directory name is given, e.g. 'in'
        if len(self.relative_directories) > 0 and len(self.paths) == 0:
            assert len(self.relative_directories) == 1
            return self.relative_directories[0]

        # Not scattered and filenames match their paths
        if not self.is_scattered and self.matching_filenames:

            first = self.paths[self.paths.keys()[0]]
            dir_path = fs.directory_of(first)

            return dir_path

        # Scatted and/or not matching filenames
        else:

            # Create temporary directory
            temp_path = introspection.create_temp_dir(time.unique_name("SKIRT_input"))

            for name in self.paths:

                path = self.paths[name]
                fs.copy_file(path, temp_path, new_name=name)

            dir_path = temp_path

            return dir_path

    # -----------------------------------------------------------------

    def to_string(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        from ..tools import formatting as fmt

        lines = []

        # Loop over the files
        for name in self.names:

            filepath = self.paths[name]
            line = line_prefix + " - " + fmt.bold + name + fmt.reset + ": " + filepath
            lines.append(line)

        # Add new line
        #lines.append(line_prefix)

        # Return
        return "\n".join(lines)

    # -----------------------------------------------------------------

    def show(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        print(self.to_string(line_prefix=line_prefix))

# -----------------------------------------------------------------

def find_input_filepath(filename, input_path):

    """
    This function ...
    :param filename:
    :param input_path:
    :return:
    """

    # List of file paths
    if types.is_sequence(input_path):

        for path in input_path:
            if fs.name(path) == filename:
                filepath = path
                break
        else: raise ValueError("The list of input paths does not contain the path to the file")

    # Directory path
    elif types.is_string_type(input_path): filepath = fs.join(input_path, filename)

    # Simulation input object
    elif isinstance(input_path, SimulationInput):

        # Check whether present in simulation input
        if filename not in input_path: raise ValueError("The file '" + filename + "' could not be found within the simulation input specification")

        # Otherwise, set the path
        filepath = input_path[filename]

    # Dictionary
    elif types.is_dictionary(input_path):

        # input_path = SimulationInput(**input_path)

        # Check whether present in simulation input
        if filename not in input_path: raise ValueError("The file '" + filename + "' could not be found within the simulation input specification")

        # Otherwise, set the path
        filepath = input_path[filename]

    # Invalid
    else: raise ValueError("Invalid value for 'input_path': '" + str(input_path) + "'")

    # Return the filepath
    return filepath

# -----------------------------------------------------------------
