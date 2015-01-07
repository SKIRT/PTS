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
    # - recursive: if one or more simulations passed to the execute function contain wildcards in their \em ski file path,
    #   the recursive argument can be used to specify whether all directories recursively nested within the base path are searched 
    #   as well, using the same filename pattern.
    # - inpath: a string specifying the absolute or relative path for simulation input files.
    # - outpath: a string specifying the absolute or relative path for simulation output files.
    # - skirel: if \c True, the simulation input/output paths are relative to the \em ski file being processed;
    #   if \c False or missing, they are relative to the current directory.
    # - threads: a positive integer specifying the number of parallel threads for each simulation; if zero or missing
    #   the number of logical cores on the computer is used.
    # - parallel: a positive integer specifying the number of simulations to be executed in parallel (by means of multithreading); 
    #   if missing the default value is one.
    # - processes: the number of (MPI) processes to be used for the execution of simulations that are to be executed with MPI
    #   (this can be checked with the MPI() function of SkirtSimulation). This argument is ignored for simulation that should not
    #   be executed with MPI, whereas the 'parallel' argument is ignored for the MPI simulations (thus multithreading is not used for
    #   running different MPI simulations in parallel). The default argument is one, in case the particular simulations are still
    #   executed with the MPI executable but with only one process. 
    # - wait: if \c True or missing, the function waits until <tt>SKIRT</tt> execution completes and sends SKIRT's brief log
    #   to the standard console; if \c False the function returns immediately without waiting for <tt>SKIRT</tt>, and
    #   <tt>SKIRT</tt>'s log messages are sent to the null device.
    #
    # The function returns a list of SkirtSimulation instances corresponding to the simulations to be performed
    # (after processing any wildcards in the ski filenames), in arbitrary order.
    #
    def execute(self, skipattern, recursive=False, inpath="", outpath="", skirel=False, threads=0, parallel=1, processes=1, wait=True):

        # In multiprocessing mode, check whether MPI is installed on the system
        if processes > 1 and not self._MPIinstalled(): return []
        
        # Create the argument list, starting with the SKIRT path in singleprocessing mode and "mpirun" in multiprocessing mode
        args = ["mpirun", "-np", str(processes), self._path, "-b"] if processes > 1 else [self._path, "-b"]
        
        # Set general command line options 
        if skirel: args += ["-k"]
        if recursive: args += ["-r"]
        if threads > 0: args += ["-t", str(threads)]
        if inpath != "": args += ["-i", inpath]
        if outpath != "": args += ["-o", outpath]
        if parallel > 1: args += ["-s", str(parallel)]
        if isinstance(skipattern, basestring): skipattern = [skipattern]
        args += skipattern

        # Execute SKIRT
        if wait:
            self._process = None
            subprocess.call(args)
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
