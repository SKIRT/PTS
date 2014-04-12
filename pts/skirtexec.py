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

from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------
#  SkirtExec class
# -----------------------------------------------------------------

## An instance of the SkirtExec class represents a particular SKIRT executable, and allows invoking it with
# given command line arguments. The current implementation only supports local execution; remote execution
# on another host will be added in the future.
class SkirtExec:

    ## The constructor accepts a single argument specifying the path to the SKIRT executable to be used.
    # The path may or may not include the "skirt" filename and if not, it may or may not include the final "/" of a
    # directory path. The path may be absolute, relative to a user's home folder, or relative to the current
    # working directory.
    # The default value is the empty path, which means SKIRT is looked for in the standard system path $PATH.
    def __init__(self, path=""):
        self._path = path
        if not self._path.endswith("skirt"): self._path = os.path.join(self._path, "skirt")
        if self._path != "skirt": self._path = os.path.realpath(os.path.expanduser(self._path))
        self._process = None

    ## This function invokes the SKIRT executable with the ski file(s) and command line options corresponding to the
    # values of the function arguments:
    # - skifiles: a single string or a sequence of strings, each specifying a relative or absolute file path for
    #   a ski file, including the \c .ski filename extension. The filename (not the base path) may also contain
    #   ? and * wildcards forming a pattern to match multiple files.
    # - recursive: if \c True, all directories recursively nested within the base path of each \em skifiles item
    #   are searched as well, using the same filename pattern; if \c False or missing, this does not happen.
    # - inpath: a string specifying the absolute or relative path for simulation input files.
    # - outpath: a string specifying the absolute or relative path for simulation output files.
    # - skirel: if \c True, the simulation input/output paths are relative to the ski file being processed;
    #   if \c False or missing, they are relative to the current directory.
    # - threads: a positive integer specifying the number of parallel threads for each simulation; if zero or missing
    #   the number of logical cores on the computer is used.
    # - simulations: a positive integer specifying the number of simulations to be executed in parallel; if missing
    #   the default value is one.
    # - wait: if \c True or missing, the function waits until SKIRT execution completes and sends SKIRT's brief log
    #   to the standard console; if \c False the function returns immediately without waiting for SKIRT, and
    #   SKIRT's log messages are sent to the null device.
    #
    # The function returns a list of SkirtSimulation instances corresponding to the simulations to be performed
    # (after processing any wildcards in the ski filenames), in arbitrary order.
    #
    def execute(self, skifiles, recursive=False, inpath="", outpath="", skirel=False,  \
                threads=0, simulations=1, wait=True):
        # construct argument list
        args = [self._path, "-b"]
        if simulations > 1: args += ["-s", str(simulations)]
        if threads > 0: args += ["-t", str(threads)]
        if skirel: args += ["-k"]
        if recursive: args += ["-r"]
        if inpath != "": args += ["-i", inpath]
        if outpath != "": args += ["-o", outpath]
        if isinstance(skifiles, basestring): skifiles = [skifiles]
        args += skifiles

        # execute skirt
        if wait:
            self._process = None
            subprocess.call(args)
        else:
            self._process = subprocess.Popen(args, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # create list of simulations
        simulations = [ ]
        for skifile in skifiles:
            root, pattern = os.path.split(skifile)
            root = os.path.realpath(root)
            # dirlist becomes a list of 3-tuples (dirpath, dirnames, filenames) where dirnames is never used
            if recursive:
                dirlist = os.walk(root)
            else:
                dirlist = [( root, [ ],  \
                          filter(lambda fn: os.path.isfile(os.path.join(root,fn)), os.listdir(root)) )]
            # process dirlist
            for dirpath, dirnames, filenames in dirlist:
                for filename in fnmatch.filter(filenames, pattern):
                    inp = os.path.join(dirpath, inpath)  if (skirel) else inpath
                    out = os.path.join(dirpath, outpath) if (skirel) else outpath
                    simulations.append( SkirtSimulation(filename, inpath=inp, outpath=out) )

        return simulations

    ## This function returns True if the previously started SKIRT process is still running, False otherwise
    def isrunning(self):
        # return true if there is a process and it has no return code yet
        return self._process != None and self._process.poll() == None

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

# -----------------------------------------------------------------
