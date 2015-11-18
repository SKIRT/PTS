#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skirtexec Executing the SKIRT command line application
#
# An instance of the SkirtExec class in this module represents a particular SKIRT executable,
# and allows invoking it with given command line arguments.

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import os.path
import subprocess

# Import the relevant PTS modules
from .log import Log
from ..core.parameters import SkirtParameters

# -----------------------------------------------------------------
#  SkirtExec class
# -----------------------------------------------------------------

## An instance of the SkirtExec class represents a particular SKIRT executable, and allows invoking it with
# given command line arguments. The current implementation only supports local execution; remote execution
# on another host will be added in the future.
class SkirtExec:

    ## The constructor accepts a single argument specifying the path to the <tt>SKIRT</tt> executable to be used.
    # The path may or may not include the "skirt" filename and if not, it may or may not include the final "/" of a
    # directory path. The path may be absolute, relative to a user's home folder, or relative to the current
    # working directory.
    # The default value is the empty path, which means <tt>SKIRT</tt> is looked for in the standard system path \c $PATH.
    def __init__(self, path="", log="", mpi_style="generic"):

        # Set the SKIRT path
        self._path = path
        if not self._path.endswith("skirt"): self._path = os.path.join(self._path, "skirt")
        if self._path != "skirt": self._path = os.path.realpath(os.path.expanduser(self._path))

        # Indicate no simulations are running yet
        self._process = None

        # Set the log mechanism
        self._log = log if log else Log()

        # Check whether MPI is present on this system
        self.mpi = mpi_present()

        # Set the MPI style
        self.mpi_style = mpi_style.lower()

    ## This function returns the command for calling SKIRT, the mpirun command is added if MPI is present on the
    #  system and if the number of processes is greater than one.
    def skirt_command(self, processes):

        # Multiprocessing mode
        if processes > 1:

            # Check whether MPI is installed on the system
            if not self.mpi:
                self._log.warning("No mpirun executable: skipping simulations")
                return None

            # Determine the command based on the MPI style
            # For 'lsf', the number of processes and hosts is derived from environment variables
            if self.mpi_style == 'lsf': return ["mpirun", "-lsf", self._path]
            else: return ["mpirun", "-np", str(processes), self._path]

        # Singleprocessing mode
        else: return [self._path]

    ## This function invokes the <tt>SKIRT</tt> executable with the simulations and command line options corresponding to the
    #  values of the function arguments:
    #
    # - simulations: the list of simulations (instances of SkirtSimulation) to be executed
    # - recursive: if one or more simulations passed to the execute function contain wildcards in their \em ski file
    #   path, the recursive argument can be used to specify whether all directories recursively nested within the base
    #   path are searched as well, using the same filename pattern.
    # - inpath: a string specifying the absolute or relative path for simulation input files.
    # - outpath: a string specifying the absolute or relative path for simulation output files.
    # - skirel: if \c True, the simulation input/output paths are relative to the \em ski file being processed;
    #   if \c False or missing, they are relative to the current directory.
    # - threads: a positive integer specifying the number of parallel threads for each simulation; if zero or missing
    #   the number of logical cores on the computer is used.
    # - parallel: a positive integer specifying the number of simulations to be executed in parallel (by means of
    #   multithreading); if missing the default value is one.
    # - processes: the number of parallel MPI processes to be launched. The default value is one, which means MPI
    #   is not used. If the specified number of processes is larger than one, the value of \em parallel argument is
    #   ignored (i.e. you can't run multiple simulations in parallel when using MPI).
    # - mpistyle: the method or style to invoke the mpirun command; currently supported values are 'generic' (the
    #   default), which uses the standard -np switch useful for launching a number of MPI processes on a single
    #   computing node, and 'lsf', which uses the -lsf switch supported by platform MPI under the LSF cluster queueing
    #   system.
    # - verbose: This option has effect only if the number of processes is larger than one. If set to \c True, each
    #   process creates its own complete log file. If missing or set to \em False, only the root process creates a
    #   full log file, and the other processes only create a log file when there are errors or warnings.
    # - wait: if \c True or missing, the function waits until <tt>SKIRT</tt> execution completes and sends the brief
    #   SKIRT log to the standard console; if \c False the function returns immediately without waiting for SKIRT,
    #   and SKIRT's log messages are sent to the null device.
    #
    # The function returns a list of SkirtSimulation instances corresponding to the simulations to be performed
    # (after processing any wildcards in the ski filenames), in arbitrary order.
    #
    def execute(self, skipattern, recursive=False, inpath=None, outpath=None, skirel=False, threads=0, parallel=1,
                processes=1, mpistyle='generic', brief=False, verbose=False, memory=False, allocation=False,
                emulate=False, single=False, wait=True, silent=False):

        self.mpi_style = mpistyle

        # Create a SkirtParameters object
        parameters = SkirtParameters()

        # The ski file pattern
        parameters.ski_pattern = skipattern
        parameters.recursive = recursive
        parameters.relative = skirel

        # Input and output
        parameters.input_path = inpath
        parameters.output_path = outpath

        # Parallelization settings
        parameters.parallel.threads = threads
        parameters.parallel.processes = processes
        parameters.parallel.simulations = parallel

        # Logging settings
        parameters.logging.brief = brief
        parameters.logging.verbose = verbose
        parameters.logging.memory = memory
        parameters.logging.allocation = allocation

        # Other settings
        parameters.emulate = emulate
        parameters.single = single

        # Run SKIRT with the specified parameters
        self.run(parameters, wait, silent)

    ## This function does the same as the execute function, but obtains its arguments from a SkirtParameters object
    def run(self, parameters, wait=True, silent=False):

        # Create the argument list
        arguments = self.skirt_command(parameters.parallel.processes)
        if arguments is None: return []

        ## Parallelization

        if parameters.parallel.threads > 0: arguments += ["-t", str(parameters.parallel.threads)]
        if parameters.parallel.simulations > 1 and parameters.parallel.processes <= 1: arguments += ["-s", str(parameters.parallel.simulations)]

        ## Logging

        if parameters.logging.brief: arguments += ["-b"]
        if parameters.logging.verbose: arguments += ["-v"]
        if parameters.logging.memory: arguments += ["-m"]
        if parameters.logging.allocation: arguments += ["-l ", str(parameters.logging.allocation_limit)]

        ## Input and output

        if parameters.input_path is not None: arguments += ["-i", parameters.input_path]
        if parameters.output_path is not None: arguments += ["-o", parameters.output_path]

        ## Other options

        if parameters.emulate: arguments += ["-e"]

        ## Ski file pattern

        if parameters.relative: arguments += ["-k"]
        if parameters.recursive: arguments += ["-r"]
        if isinstance(parameters.ski_pattern, basestring): arguments += [parameters.ski_pattern]
        elif isinstance(parameters.ski_pattern, list): arguments += parameters.ski_pattern
        else: raise ValueError("The ski pattern must consist of either a string or a list of strings")

        # Launch the SKIRT command
        if wait:
            self._process = None
            if silent: subprocess.call(" ".join(arguments), shell=True, stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w'))
            else: subprocess.call(" ".join(arguments), shell=True)
        else: self._process = subprocess.Popen(arguments, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # Return the list of simulations so that their results can be followed up
        return parameters.simulations()

    ## This function returns True if the previously started SKIRT process is still running, False otherwise
    def isrunning(self):
        return (self._process != None and self._process.poll() == None)

    ## This function waits for the previously started SKIRT process to complete, if needed
    def wait(self):
        if self.isrunning(): self._process.wait()

    ## This function returns a string with version information on the SKIRT executable represented by this
    # object. The function invokes SKIRT with an incorrect command line argument to obtain this information.
    def version(self):
        # execute skirt with incorrect argument list and get its output
        process = subprocess.Popen([self._path, "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = process.communicate()[0];

        # return the relevant portion of the output
        return "SKIRT" + output.splitlines()[0].partition("SKIRT")[2]

    ## This function returns a string with the path of the SKIRT directory that is used, for example: HOME/SKIRT
    def directory(self):
        return os.path.dirname(os.path.dirname(os.path.dirname(self._path)))

# -----------------------------------------------------------------

## An instance of the FitSkirtExec class represents a particular FitSKIRT executable, and allows invoking it
#  with given command line arguments.
class FitSkirtExec(object):

    ## The constructor ...
    def __init__(self, path=None):

        # Set the FitSKIRT path
        self._path = path
        if not self._path.endswith("fitskirt"): self._path = os.path.join(self._path, "fitskirt")
        if self._path != "fitskirt": self._path = os.path.realpath(os.path.expanduser(self._path))

        # Indicate no simulations are running yet
        self._process = None

    ## This function ...
    def execute(self, fskipattern, recursive=False, inpath="", outpath="", skirel=False, threads=0, parallel=1,
                processes=1, mpistyle='generic', brief=False, verbose=False, memory=False, wait=True):

        pass

    ## This function ...
    def run(self, config):

        pass

# *****************************************************************

def mpi_present():

    """
    This function checks whether the MPI executable is installed on the system
    :return:
    """

    try:
        devnull = open(os.devnull)
        subprocess.Popen("mpirun", stdout=devnull, stderr=devnull).communicate()
        return True
    except: return False

# -----------------------------------------------------------------
