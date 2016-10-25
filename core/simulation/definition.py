#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.definition Contains the SingleSimulationDefinition and MultiSimulationDefinition classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from .skifile import SkiFile

# -----------------------------------------------------------------

def create_definitions(path, output_path, input_path, recursive=False):

    """
    This function ...
    :param path:
    :param output_path:
    :param input_path:
    :return:
    """

    definitions = []

    # If ski files don't have to be found recursively (in seperate subdirectories)
    if not self.config.recursive:

        # Create an 'out' directory if the output directory is not specified
        if self.config.output is None: output_path = fs.create_directory_in(self.config.path, "out")
        else: output_path = fs.absolute(self.config.output)

    else: output_path = None

    # Keep track of the directories where ski files were found
    ski_dir_paths = []

    # Loop over all files in the current working directory
    for ski_path, prefix in fs.files_in_path(self.config.path, extension="ski", returns=["path", "name"], recursive=self.config.recursive):

        # Determine the path of the directory in which the ski file is found
        dir_path = fs.directory_of(ski_path)

        # Open the ski file and check whether input is required
        ski = SkiFile(ski_path)
        if ski.needs_input:
            input_paths = ski.input_paths(self.config.input, self.config.path)
        else:
            input_paths = None

        # Determine output directory
        if self.config.recursive:

            # Check if the ski file directory is not yet encountered (multiple ski files in a directory)
            if dir_path in ski_dir_paths: raise RuntimeError("There can't be multiple ski files in a directory in recursive mode")

            # Determine output directory
            simulation_output_path = fs.create_directory_in(dir_path, "out")

            # Add the ski directory path
            ski_dir_paths.append(dir_path)

        # Create a directory in the output directory
        else:
            simulation_output_path = fs.create_directory_in(output_path, prefix)

        # Create the simulation definition
        definition = SingleSimulationDefinition(ski_path, simulation_output_path, input_paths)

        # Add the definition
        definitions.append(definition)

    # Return the list of definitions
    return definitions

# -----------------------------------------------------------------

class SingleSimulationDefinition(object):

    """
    This class ...
    """

    def __init__(self, ski_path, output_path, input_path=None, name=None):

        """
        The constructor ...
        :param ski_path:
        :param output_path:
        :param input_path: can be path to input directory or list of file paths
        :param name:
        :return:
        """
        
        # Options for the ski file pattern
        self.ski_path = ski_path

        # The input and output paths
        self.input_path = input_path
        self.output_path = output_path

        # A name for this simulation
        self.name = name
            
    # -----------------------------------------------------------------

    @property
    def base_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.ski_path)

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        return fs.strip_extension(fs.name(self.ski_path))

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        """

        properties = []
        properties.append("ski path: " + self.ski_path)
        if self.input_path is not None: properties.append("input path(s): " + str(self.input_path))
        properties.append("output path: " + str(self.output_path))
    
        return_str = self.__class__.__name__ + ":\n"
        for property in properties: return_str += " -" + property + "\n"
        return return_str

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        """

        return '<' + self.__class__.__name__ + " ski path: '" + self.ski_path + "'>"

# -----------------------------------------------------------------

class MultiSimulationDefinition(object):

    """
    This class ...
    """

    def __init__(self, base_path, pattern, input_name, output_name, recursive=True):

        """
        The constructor ...
        :param base_path: the directory with the ski files
        :param pattern: the ski pattern
        :param input_name: the name of the input directory (or a list of file paths)
        :param output_name: the name of the output directory
        :param recursive:
        """

        self.base_path = base_path
        self.pattern = pattern
        self.input = input_name
        self.output = output_name
        self.recursive = recursive

# -----------------------------------------------------------------
