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

# -----------------------------------------------------------------

import fnmatch
import os
import os.path
import subprocess

from pts.log import Log
from pts.skirtsimulation import SkirtSimulation

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
    def __init__(self, path="", log=""):

        # Set the SKIRT path
        self._path = path
        if not self._path.endswith("skirt"): self._path = os.path.join(self._path, "skirt")
        if self._path != "skirt": self._path = os.path.realpath(os.path.expanduser(self._path))

        # Indicate no simulations are running yet
        self._process = None

        # Set the log mechanism
        self._log = log if log else Log()

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
    #   and SKIRT's log messages are sent to the null device. If \c True, SKIRT's log messages are logged to the console
    #   except when the 'console' argument is set to \c False.
    # - console: this argument only applies when the 'wait' argument is set to \c True. It specifies whether the SKIRT
    #   should log the console output or whether SKIRT's log messages are sent to the null device.
    #
    # The function returns a list of SkirtSimulation instances corresponding to the simulations to be performed
    # (after processing any wildcards in the ski filenames), in arbitrary order.
    #
    def execute(self, skipattern, recursive=False, inpath="", outpath="", skirel=False,
                threads=0, parallel=1, processes=1, mpistyle='generic', brief=False, verbose=False, wait=True, console=True):

        # In multiprocessing mode, check whether MPI is installed on the system
        if processes > 1 and not self._MPIinstalled(): return []

        # Create the argument list, starting with the SKIRT path in singleprocessing mode and "mpirun" in multiprocessing mode
        args = [self._path]
        if processes > 1:
            if mpistyle.lower()=='lsf':
                args = ["mpirun", "-lsf"] + args  # the number of processes and hosts is derived from environment variables
            else:
                args = ["mpirun", "-np", str(processes)] + args

        # Set general command line options
        if brief: args += ["-b"]
        if verbose: args += ["-v"]
        if skirel: args += ["-k"]
        if recursive: args += ["-r"]
        if threads > 0: args += ["-t", str(threads)]
        if inpath != "": args += ["-i", inpath]
        if outpath != "": args += ["-o", outpath]
        if parallel > 1 and processes<=1: args += ["-s", str(parallel)]
        if isinstance(skipattern, basestring): skipattern = [skipattern]
        args += skipattern

        # Execute SKIRT
        if wait:
            self._process = None
            subprocess.call(args) if console else subprocess.call(args, stdout=open(os.path.devnull, 'w'))
        else:
            self._process = subprocess.Popen(args, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # Create a list of the simulations that are executed
        simulations = []
        for skifile in skipattern:
            root, pattern = os.path.split(skifile)
            root = os.path.realpath(root)
            # The "dirlist" variable becomes a list of 3-tuples (dirpath, dirnames, filenames) where dirnames is never used
            if recursive:
                dirlist = os.walk(root)
            else:
                dirlist = [( root, [ ], filter(lambda fn: os.path.isfile(os.path.join(root,fn)), os.listdir(root)) )]
            # Process "dirlist"
            for dirpath, dirnames, filenames in dirlist:
                for filename in fnmatch.filter(filenames, pattern):
                    inp = os.path.join(dirpath, inpath)  if (skirel) else inpath
                    out = os.path.join(dirpath, outpath) if (skirel) else outpath
                    simulations.append(SkirtSimulation(filename, inpath=inp, outpath=out))

        # Return the list of simulations so that their results can be followed up
        return simulations

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

    ## This function checks whether the MPI executable is installed on the system.
    def _MPIinstalled(self):
        MPI = True
        try:
            devnull = open(os.devnull)
            subprocess.Popen("mpirun", stdout=devnull, stderr=devnull).communicate()
        except:
            self._log.warning("No mpirun executable: skipping MPI simulations!")
            MPI = False
        return MPI

# -----------------------------------------------------------------
