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

# Import standard modules
import copy

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from .skifile import SkiFile
from ..tools.utils import lazyproperty
from .tree import DustGridTree, get_nleaves

# -----------------------------------------------------------------

def create_definitions(path, output_path, input_path, recursive=False):

    """
    This function ...
    :param path:
    :param output_path:
    :param input_path:
    :param recursive:
    :return:
    """

    definitions = []

    # If ski files don't have to be found recursively (in seperate subdirectories)
    if not recursive:
        # Create an 'out' directory if the output directory is not specified
        if output_path is None: output_path = fs.create_directory_in(path, "out")
        else: output_path = fs.absolute_path(output_path)
    else: output_path = None

    # Keep track of the directories where ski files were found
    ski_dir_paths = []

    # Loop over all files in the current working directory
    for ski_path, prefix in fs.files_in_path(path, extension="ski", returns=["path", "name"], recursive=recursive):

        # Determine the path of the directory in which the ski file is found
        dir_path = fs.directory_of(ski_path)

        # Open the ski file and check whether input is required
        ski = SkiFile(ski_path)
        if ski.needs_input: input_paths = ski.input_paths(input_path, path)
        else: input_paths = None

        # Determine output directory
        if recursive:

            # Check if the ski file directory is not yet encountered (multiple ski files in a directory)
            if dir_path in ski_dir_paths: raise RuntimeError("There can't be multiple ski files in a directory in recursive mode")

            # Determine output directory
            simulation_output_path = fs.create_directory_in(dir_path, "out")

            # Add the ski directory path
            ski_dir_paths.append(dir_path)

        # Create a directory in the output directory
        else: simulation_output_path = fs.create_directory_in(output_path, prefix)

        # Create the simulation definition
        definition = SingleSimulationDefinition(ski_path, simulation_output_path, input_paths)

        # Add the definition
        definitions.append(definition)

    # Return the list of definitions
    return definitions

# -----------------------------------------------------------------

class SimulationDefinition(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This constructor does not implement anything ...
        """

        pass

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------

class SingleSimulationDefinition(SimulationDefinition):

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

        # Call the constructor of the base class
        super(SingleSimulationDefinition, self).__init__()
        
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
        return fs.directory_of(self.ski_path)

    # -----------------------------------------------------------------

    @property
    def prefix(self):
        return fs.strip_extension(fs.name(self.ski_path))

    # -----------------------------------------------------------------

    @lazyproperty
    def ski(self):
        return SkiFile(self.ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):
        return self.ski.nwavelengthsfile(self.input_path) if self.ski.wavelengthsfile() else self.ski.nwavelengths()

    # -----------------------------------------------------------------

    @lazyproperty
    def dustlib_dimension(self):
        return self.ski.dustlib_dimension()

    # -----------------------------------------------------------------

    @lazyproperty
    def uses_dust_grid_tree(self):
        return self.ski.filetreegrid()

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree_path(self):
        if not self.uses_dust_grid_tree: return None
        else: return self.ski.treegridfile(self.input_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):
        if not self.uses_dust_grid_tree: return None
        else: return DustGridTree.from_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_cells(self):
        if self.ski.treegrid_notfile(): return None # number of cells cannot be defined
        elif self.ski.filetreegrid(): return get_nleaves(self.dust_grid_tree_path) #return self.dust_grid_tree.nleaves
        else: return self.ski.ncells()

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

class MultiSimulationDefinition(SimulationDefinition):

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

        # Call the constructor of the base class
        super(MultiSimulationDefinition, self).__init__()

        # Set properties
        self.base_path = base_path
        self.pattern = pattern
        self.input = input_name
        self.output = output_name
        self.recursive = recursive

# -----------------------------------------------------------------
