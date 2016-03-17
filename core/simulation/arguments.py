#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.arguments Contains the SkirtArguments class, used for representing the set of
#  command-line arguments that can be passed to SKIRT.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import fnmatch

# Import the relevant PTS classes and modules
from .simulation import SkirtSimulation
from ..basics.map import Map

# -----------------------------------------------------------------

class SkirtArguments(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # If the set of arguments is passed to the constructor as a configuration mapping, set the
        # attributes of this object according to the configuration settings
        if config is not None:

            # Loop over the entries in the configuration and create an attribute with the same name
            for entry in config: setattr(self, entry, config[entry])

        # If no configuration mapping is passed to this object, set all attributes to default values
        else:

            ## Set the attributes of the object to default values

            # Options for the ski file pattern
            self.ski_pattern = None
            self.recursive = None
            self.relative = None

            # The input and output paths
            self.input_path = None
            self.output_path = None

            # Other options
            self.emulate = False    # Run in emulation mode
            self.single = False     # True if only a single simulation is expected

            # Options for logging
            self.logging = Map()
            self.logging.brief = False            # Brief console logging
            self.logging.verbose = False          # Verbose logging
            self.logging.memory = False           # State the amount of used memory with each log message
            self.logging.allocation = False       # Write log messages with the amount of (de)allocated memory
            self.logging.allocation_limit = 1e-5  # The lower limit for the amount of (de)allocated memory to be logged

            # Options for parallelization
            self.parallel = Map()
            self.parallel.simulations = None  # The number of parallel simulations
            self.parallel.threads = None      # The number of parallel threads per simulation
            self.parallel.processes = None    # The number of parallel processes per simulation

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

    def to_command(self, skirt_path, mpi_command, scheduler, bind_to_cores=False, threads_per_core=1, to_string=False):

        """
        This function ...
        :param skirt_path:
        :param mpi_command:
        :param scheduler:
        :param bind_to_cores:
        :param threads_per_core:
        :param to_string:
        :return:
        """

        # Create the argument list
        arguments = skirt_command(skirt_path, mpi_command, bind_to_cores, self.parallel.processes, self.parallel.threads, threads_per_core, scheduler)

        # Parallelization options
        if self.parallel.threads > 0: arguments += ["-t", str(self.parallel.threads)]
        if self.parallel.simulations > 1 and self.parallel.processes <= 1: arguments += ["-s", str(self.parallel.simulations)]

        # Logging options
        if self.logging.brief: arguments += ["-b"]
        if self.logging.verbose: arguments += ["-v"]
        if self.logging.memory: arguments += ["-m"]
        if self.logging.allocation: arguments += ["-l", str(self.logging.allocation_limit)]

        # Options for input and output
        if self.input_path is not None: arguments += ["-i", self.input_path]
        if self.output_path is not None: arguments += ["-o", self.output_path]

        # Other options
        if self.emulate: arguments += ["-e"]

        # Ski file pattern
        if self.relative: arguments += ["-k"]
        if self.recursive: arguments += ["-r"]
        if isinstance(self.ski_pattern, basestring): arguments += [self.ski_pattern]
        elif isinstance(self.ski_pattern, list): arguments += self.ski_pattern
        else: raise ValueError("The ski pattern must consist of either a string or a list of strings")

        # If requested, convert the argument list into a string
        if to_string:

            # Create the final command string for this simulation
            command = " ".join(arguments)

            # Return the command string
            return command

        # Otherwise, return the list of argument values
        else: return arguments

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function creates a copy of this SkirtArguments object
        :return:
        """

        # Create a new SkirtArguments object
        arguments = SkirtArguments()

        ## Set options identical to this instance

        # Options for the ski file pattern
        arguments.ski_pattern = self.ski_pattern
        arguments.recursive = self.recursive
        arguments.relative = self.relative

        # The input and output paths
        arguments.input_path = self.input_path
        arguments.output_path = self.output_path

        # Other options
        arguments.emulate = self.emulate    # Run in emulation mode
        arguments.single = self.single     # True if only a single simulation is expected

        # Options for logging
        arguments.logging.brief = self.logging.brief            # Brief console logging
        arguments.logging.verbose = self.logging.verbose        # Verbose logging
        arguments.logging.memory = self.logging.memory          # State the amount of used memory with each log message
        arguments.logging.allocation = self.logging.allocation  # Write log messages with the amount of (de)allocated memory
        arguments.logging.allocation_limit = self.logging.allocation_limit  # The lower limit for the amount of (de)allocated memory to be logged

        # Options for parallelization
        arguments.parallel.simulations = self.parallel.simulations  # The number of parallel simulations
        arguments.parallel.threads = self.parallel.threads          # The number of parallel threads per simulation
        arguments.parallel.processes = self.parallel.processes      # The number of parallel processes per simulation

        # Return the new object
        return arguments

# -----------------------------------------------------------------

def skirt_command(skirt_path, mpi_command, bind_to_cores, processes, threads, threads_per_core, scheduler):

    """
    This function ...
    :param skirt_path:
    :param mpi_command:
    :param bind_to_cores:
    :param processes:
    :param threads:
    :param threads_per_core:
    :param scheduler:
    :return:
    """

    # Multiprocessing mode
    if processes > 1:

        # Determine the command based on whether or not a scheduling system is used
        if scheduler: command = mpi_command.split()
        else: command = mpi_command.split() + ["-np", str(processes)]

        # If 'process to core' binding must be enabled, add the 'cpus-per-proc' option
        # (see https://www.open-mpi.org/faq/?category=tuning)
        if bind_to_cores:
            # Hyperthreading: threads_per_core will be > 1
            # No hyperthreading: threads_per_core will be 1
            # cores / process = (cores / thread) * (threads / process)
            cores_per_process = threads / threads_per_core
            command += ["--cpus-per-proc", str(cores_per_process)] # "CPU'S per process" means "core per process" in our definitions

        # Add the SKIRT path and return the final command list
        command += [skirt_path]
        return command

    # Singleprocessing mode, no MPI command or options
    else: return [skirt_path]

# -----------------------------------------------------------------
