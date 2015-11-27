#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import fnmatch

# Import the relevant PTS classes and modules
from ..tools import configuration
from .simulation import SkirtSimulation

# -----------------------------------------------------------------

class SkirtParameters(object):

    """
    This class ...
    """

    def __init__(self, parameters=None):

        """
        The constructor ...
        :return:
        """

        # Load a configuration object according to the default template and the passed parameters
        config = configuration.set("skirt", parameters)

        # Set the configuration entries as attributes of the object
        for entry in config: setattr(self, entry, config[entry])

    # -----------------------------------------------------------------

    def simulations(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the simulation objects
        simulations = []

        # Loop over the seperate ski files defined in the ski pattern
        pattern = [self.ski_pattern] if isinstance(self.ski_pattern, basestring) else self.ski_pattern
        for skifile in pattern:

            # Determine the directory path and the actual file descriptor
            root, name = os.path.split(skifile)
            root = os.path.realpath(root)

            # Construct the 'dirlist' variable; this is a list of 3-tuples (dirpath, dirnames, filenames)
            if self.recursive: dirlist = os.walk(root)
            else: dirlist = [(root, [], filter(lambda fn: os.path.isfile(os.path.join(root,fn)), os.listdir(root)))]

            # Search for ski files matching the pattern and construct SkirtSimulation objects
            for dirpath, dirnames, filenames in dirlist:
                for filename in fnmatch.filter(filenames, name):

                    inp = os.path.join(dirpath, self.input_path)  if (self.relative) else self.input_path
                    out = os.path.join(dirpath, self.output_path) if (self.relative) else self.output_path

                    # Create the simulation and add it to the list
                    simulations.append(SkirtSimulation(filename, inpath=inp, outpath=out))

        # Check whether the ski pattern is ought to represent only one particular simulation
        if self.single:

            # If multiple matching ski files are found, raise an error
            if len(simulations) > 1: raise ValueError("The specified ski pattern defines multiple simulations")
            else: return simulations[0]

        # Else, just return the list of simulations (even when containing only one item)
        else: return simulations

    # -----------------------------------------------------------------

    def to_command(self, skirt_path):

        """
        This function ...
        :return:
        """

        # Create the argument list
        arguments = skirt_command(self.parallel.processes)

        ## Parallelization

        if self.parallel.threads > 0: arguments += ["-t", str(self.parallel.threads)]
        if self.parallel.simulations > 1 and self.parallel.processes <= 1: arguments += ["-s", str(self.parallel.simulations)]

        ## Logging

        if self.logging.brief: arguments += ["-b"]
        if self.logging.verbose: arguments += ["-v"]
        if self.logging.memory: arguments += ["-m"]
        if self.logging.allocation: arguments += ["-l ", str(self.logging.allocation_limit)]

        ## Input and output

        if self.input_path is not None: arguments += ["-i", self.input_path]
        if self.output_path is not None: arguments += ["-o", self.output_path]

        ## Other options

        if self.emulate: arguments += ["-e"]

        ## Ski file pattern

        if self.relative: arguments += ["-k"]
        if self.recursive: arguments += ["-r"]
        if isinstance(self.ski_pattern, basestring): arguments += [self.ski_pattern]
        elif isinstance(self.ski_pattern, list): arguments += self.ski_pattern
        else: raise ValueError("The ski pattern must consist of either a string or a list of strings")

        # Create the final command string for this simulation
        command = " ".join(arguments)

        # Return the command string
        return command

# -----------------------------------------------------------------

def skirt_command(skirt_path, mpi_command, processes, scheduler):

    """
    This function ...
    :param processes:
    :return:
    """

    # Multiprocessing mode
    if processes > 1:

        # Determine the command based on whether or not a scheduling system is used
        if scheduler: return [mpi_command, skirt_path]
        else: return [mpi_command, "-np", str(processes), skirt_path]

    # Singleprocessing mode
    else: return [skirt_path]

# -----------------------------------------------------------------

class FitSkirtParameters(object):

    """
    This class ...
    """

    def __init__(self, parameters):

        """
        The constructor ...
        :param parameters:
        :return:
        """

        # Load a configuration object according to the default template and the passed parameters
        config = configuration.set("fitskirt", parameters)

        # Set the configuration entries as attributes of the object
        for entry in config: setattr(self, entry, config[entry])

    # -----------------------------------------------------------------

    def simulations(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------