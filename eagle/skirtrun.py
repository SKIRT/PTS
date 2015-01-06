#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.skirtrun Managing the files related to a particular SKIRT run for EAGLE.
#
# An instance of the SkirtRun class in this module manages the files related to a particular SKIRT run in
# the context of the EAGLE project. It also provides access to corresponding SkirtSimulation and SkirtExec objects.
#

# -----------------------------------------------------------------

import os
import os.path
import eagle.config as config
from pts.skirtsimulation import SkirtSimulation
from pts.skirtexec import SkirtExec

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
    # If the \em create flag is true, the relevant directories are created (including the input, output
    # and visualization subdirectories).
    def __init__(self, runid, create=False):
        self._runid = int(runid)
        if self._runid<=0 or self._runid>9999999:
            raise ValueError("SKIRT run-id exceeds implementation limits: " + str(runid))
        topdir = "g-{0:04}".format(self._runid//1000)
        botdir = "r-{0:04}-{1:03}".format(self._runid//1000, self._runid%1000)
        self._runpath = os.path.join(config.results_path, topdir, botdir)

        if create:
            try:
                os.makedirs(self.inpath())
                os.makedirs(self.outpath())
                os.makedirs(self.vispath())
            except:
                if (not os.path.isdir(self.inpath())): raise ValueError("Can't create directory " + self.inpath())
                if (not os.path.isdir(self.outpath())): raise ValueError("Can't create directory " + self.outpath())
                if (not os.path.isdir(self.vispath())): raise ValueError("Can't create directory " + self.vispath())

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

    ## This function returns a SkirtSimulation object for the <tt>SKIRT</tt> results stored in the output subdirectory
    # for this <tt>SKIRT</tt> run. An error is raised if the run directory does not contain exactly one \em ski file.
    def simulation(self):
        return SkirtSimulation(self.prefix(), self.inpath(), self.outpath())

    ## This function invokes the <tt>SKIRT</tt> executable with the \em ski file for this <tt>SKIRT</tt> run using 
    # the appropriate input and output paths, waits until <tt>SKIRT</tt> execution completes, and returns a 
    # SkirtSimulation object for the <tt>SKIRT</tt> results. An error is raised if the run directory does not contain 
    # exactly one \em ski file. The \em threads argument specifies the number of parallel threads for each simulation; 
    # if zero or missing the number of logical cores on the computer is used.
    def execute(self, threads=0):
        
        # Create the SKIRT execution environment
        skirt = SkirtExec(config.skirt_path)
        
        # Create the simulation
        skifile = os.path.join(self._runpath, self.prefix()+".ski")
        simulation = SkirtSimulation(self.prefix(), inpath=self.inpath(), outpath=self.outpath(), skifile=skifile)
        
        # Execute the simulation
        skirt.execute([simulation], inpath=self.inpath(), outpath=self.outpath(), threads=threads)
        return simulation

# -----------------------------------------------------------------
