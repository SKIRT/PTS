#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.skirtrun Managing the files related to a particular SKIRT run for EAGLE.
#
# An instance of the SkirtRun class in this module manages the files related to a particular SKIRT run in
# the context of the EAGLE project. It also provides access to corresponding SkirtSimulation and SkirtExec objects.
#

# -----------------------------------------------------------------

import os
import os.path

from . import config as config
from ..core.simulation.simulation import SkirtSimulation
from ..core.simulation.execute import SkirtExec

# -----------------------------------------------------------------

## This function returns a list of run-ids (integer numbers) corresponding to the given range specification,
# or None in case of syntax error. The range specification must be a comma-seperated list of run-ids and/or
# run-id ranges expressed as two run-ids with a dash in between.
def runids_in_range(runidspec):
    try:
        runids = []
        for segment in runidspec.split(","):
            if "-" in segment:
                first,last = map(int,segment.split("-"))
                runids += [ id for id in range(first,last+1) ]
            else:
                if segment!="": runids += [ int(segment) ]
        return runids
    except Exception:
        return None

# -----------------------------------------------------------------
#  SkirtRun class
# -----------------------------------------------------------------

## An instance of the SkirtRun class manages the files related to a particular SKIRT run in the context of
# the EAGLE project, and provides access to corresponding SkirtSimulation and SkirtExec objects.
#
# Many file systems have performance problems with directories that contain a large number of files or subdirectories.
# To avoid these issues, the skirt results are placed in directory hierarchy based on the unique SKIRT run identifier
# defined in the SKIRT-runs database.
# The current implementation assumes that the run-id is a positive integer below 10 million and uses two levels with
# at total of 10000 x 1000 directories to organize the SKIRT run data. If this implementation limitation is exceeded
# we can switch to a more complex hierarchy without too much hassle.
#
# The bottom-level directories are named according to the pattern "r-nnnn-nnn" where the literal \em r stands for
# "run" and the \em n stand for the consecutive digits of the run-id. The dashes are included for readability.
# The top-level directories are named according to the pattern "g-nnnn" where the literal \em r stands for "group"
# and the \em n stand for the most significant digits of the run-id. Each top-level directory contains the
# 1000 bottom-level directories with corresponding most-significant digits.
#
# Inside each run directory the data is organized as follows:
#  - the ski file is placed immediately inside the run directory
#  - input data for SKIRT is placed in the "in" subdirectory
#  - output data from SKIRT is placed in the "out" subdirectory
#  - data derived from the SKIRT results (e.g. for visualization) are placed in the "vis" subdirectory
#
# This structure facilitates managing the data files. For example to save disk space one could remove
# the contents of the input and visualization directories (since this information is easily regenerated
# respectively from the EAGLE hdf5 files and from the SKIRT results), and one could compress the contents
# of the SKIRT output directory into a ZIP file.
#
class SkirtRun:

    ## The constructor takes a SKIRT run-id and determines the corresponding absolute directory path.
    # By default the SKIRT-run directory will be located in the results path defined in the current configuration;
    # this can be overriden by specifying a value for the \em alternate_results_path argument.
    # If the \em create flag is true, the relevant directories are created (including the input, output
    # and visualization subdirectories).
    def __init__(self, runid, create=False, alternate_results_path=None):
        self._runid = int(runid)
        if self._runid<=0 or self._runid>9999999:
            raise ValueError("SKIRT run-id exceeds implementation limits: " + str(runid))
        topdir = "g-{0:04}".format(self._runid//1000)
        botdir = "r-{0:04}-{1:03}".format(self._runid//1000, self._runid%1000)
        resultspath = config.results_path if alternate_results_path==None else alternate_results_path
        self._runpath = os.path.join(resultspath, topdir, botdir)

        if create:
            _createdir(self.inpath())
            _createdir(self.outpath())
            _createdir(self.vispath())

    ## This function returns the SKIRT run identifier that was passed to the constructor
    def runid(self):
        return self._runid

    ## This function returns an absolute path to the directory where the ski file should be placed for this SKIRT run.
    def runpath(self):
        return self._runpath

    ## This function returns an absolute path to the input subdirectory for this SKIRT run.
    def inpath(self):
        return os.path.join(self._runpath, "in")

    ## This function returns an absolute path to the output subdirectory for this SKIRT run.
    def outpath(self):
        return os.path.join(self._runpath, "out")

    ## This function returns an absolute path to the visualization subdirectory for this SKIRT run.
    def vispath(self):
        return os.path.join(self._runpath, "vis")

    ## This function returns true if the run directory contains exactly one ski file, false if not
    def hasskifile(self):
        skifiles = filter(lambda fn: fn.endswith(".ski"), os.listdir(self._runpath))
        return len(skifiles)==1

    ## This function returns the name of the ski file stored in the run directory for this SKIRT run.
    # (without the ".ski" filename extension)
    # An error is raised if the run directory does not contain exactly one ski file.
    def prefix(self):
        skifiles = filter(lambda fn: fn.endswith(".ski"), os.listdir(self._runpath))
        if len(skifiles)!=1:
            raise ValueError("SKIRT run directory does not contain exactly one ski file: " + self._runpath)
        skifile = skifiles[0]
        return skifile[0:-len(".ski")]

    ## This function returns a SkirtSimulation object for the SKIRT results stored in the output subdirectory
    # for this SKIRT run. An error is raised if the run directory does not contain exactly one ski file.
    def simulation(self):
        return SkirtSimulation(self.prefix(), self.inpath(), self.outpath())

    ## This function invokes the SKIRT executable with the ski file for this SKIRT run using the appropriate
    # input and output paths, waits until SKIRT execution completes, and returns a SkirtSimulation object
    # for the SKIRT results. An error is raised if the run directory does not contain exactly one ski file.
    #
    # The \em processes argument specifies the number of parallel MPI processes to be launched; to disable MPI
    # specify a value of one, or omit the argument.
    # The \em threads argument specifies the number of parallel threads for each simulation; if zero or missing
    # the number of logical cores on the computer is used.
    # The \em mpistyle argument specifies the style to invoke the mpirun command; default is 'generic'.
    def execute(self, processes=1, threads=0, mpistyle='generic'):
        skirt = SkirtExec(config.skirt_path)
        skifile = os.path.join(self._runpath, self.prefix()+".ski")
        simulations = skirt.execute(skifile, inpath=self.inpath(), outpath=self.outpath(),
                      processes=processes, threads=threads, mpistyle=mpistyle, brief=True, verbose=True)
        return simulations[0]

# -----------------------------------------------------------------

## This private helper function creates directories along the specified path, and raises an error if there is a problem
def _createdir(path):
    # makedirs() fails if the directory already exists, so we catch exceptions and verify the result
    try:
        os.makedirs(path)
    except Exception:
        if (not os.path.isdir(path)): raise ValueError("Can't create directory " + path)

# -----------------------------------------------------------------
