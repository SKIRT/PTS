#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.scalingtest Performing a scaling test for SKIRT
#
# An instance of the ScalingTest class represents a SKIRT scaling benchmark test.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import multiprocessing
import datetime

# Import the relevant PTS class
try:
    from pts.skirtexec import SkirtExec
except ImportError:
    print "the pts.skirtexec module could not be loaded (this is not necessary on the cluster)"
from pts.log import Log
from do.extractscaling import extract
from pts.jobscript import JobScript

# -----------------------------------------------------------------

## An instance of the ScalingTest class represents a SKIRT scaling benchmark test for a particular ski file.
#
class ScalingTest:

    ## The constructor accepts the following arguments:
    #
    #  - path: the path of the directory used for SKIRT scaling benchmark tests. This directory should contain
    #          different folders, each containing a particular ski file.
    #  - simulation: this string indicates which of the simulations or subfolders of the path (see above) should
    #                be used for this scaling test. This string should match the name of one of the subdirectories of
    #                the path.
    #  - system: a string identifying the system on which this scaling test is carried out. This system name is used
    #            in the name of the scaling test results file. If the system name has a particular value, for example
    #            "delcatty" or "cosma", the scaling test will automatically switch to a mode where simulations are
    #            scheduled instead of executed immediately.
    #  - mode: the mode in which to run the scaling test. This can be either "mpi", "threads" or "hybrid".
    #
    def __init__(self, path, simulation, system, mode):
        
        # Set the path
        self._path = path

        # Set the simulation name (subdirectory of self._path) and the simulation path
        self._simulationname = simulation
        self._simulationpath = os.path.join(self._path, self._simulationname)

        # Set the output and result paths
        self._outpath = os.path.join(self._simulationpath, "out")
        self._respath = os.path.join(self._simulationpath, "res")
        
        # Create the output directories if they do not already exist
        try: os.mkdir(self._outpath)
        except OSError: pass
        try: os.mkdir(self._respath)
        except OSError: pass
        
        # Set the system name
        self._system = system

        # Determine whether we are dealing with a scheduling system or we can launch simulations right away
        if self._system == "delcatty" or self._system == "cosma":
            self._do = self._schedule
        else:
            self._do = self._run

        # Set the mode
        self._mode = mode

        # Determine the number of cores (per node) on this system
        self._cores = multiprocessing.cpu_count()

        # TODO: Search this file automatically; support other names
        self._skifilename = "scaling"

        # The path of the ski file
        self._skifilepath = os.path.join(self._simulationpath, self._skifilename + ".ski")

        # Create the logging mechanism
        self._log = Log()
        
        # Create skirt execution context
        try:
            self._skirt = SkirtExec(log=self._log)
        except NameError:
            self._skirt = None

    ## When this function is invoked, the scaling test is started. This function takes the following arguments:
    #
    #  - maxnodes: the maximum number of 'nodes' to be used for this scaling test. On a desktop system, a node is
    #              considered to be the entire set of processors. This number can range from zero to infinity, and
    #              can be a fractional number. To use 2 processors on a system with 4 processors (per node) for example,
    #              a value of 0.5 can be used for maxnodes.
    #  - keepoutput: this optional argument indicates whether the output of the SKIRT simulations has to be kept
    #                or can be deleted from the disk.
    #
    def run(self, maxnodes, minnodes, keepoutput=False):

        # Log the system name, the test mode and the version of SKIRT used for this test
        self._log.info("Starting parallel scaling benchmark for " + self._system + " in " + self._mode + " mode.")
        if self._skirt is not None:
            self._log.info("Using " + self._skirt.version())

        # Create a file containing the results of the scaling test
        resultsfilepath = self._createresultsfile(maxnodes, minnodes)

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        maxprocessors = int(maxnodes * self._cores)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one.
        minprocessors = int(minnodes * self._cores)
        processors = minprocessors
        if processors == 0: processors = 1

        # Perform the simulations with increasing number of processors
        while processors <= maxprocessors:

            # Perform this run
            self._do(processors, resultsfilepath, keepoutput)

            # The next run will be performed with double the amount of processors
            processors *= 2

        # End with some log messages
        self._log.success("Finished parallel scaling benchmark script")
        self._log.info("The results are / will be written to " + resultsfilepath)

    ## This functions schedules a simulation on the cluster. This function takes the following arguments:
    #
    #   - processors: the total number of processors to be used for this run
    #   - resultsfilepath: the path of the file that should contain the timings from this run
    #   - keepoutput: a flag indicating whether the SKIRT output should be kept or deleted
    #
    def _schedule(self, processors, resultsfilepath, keepoutput):

        # Inform the user about the number of processes and threads used for this run
        self._log.info("Scheduling simulation with " + str(processors) + " processors")

        # Calculate the necessary amount of nodes
        nodes = processors/self._cores + (processors % self._cores > 0)

        # Determine the number of processors per node
        ppn = processors if nodes == 1 else self._cores

        # Scaling benchmark name (="SKIRTscaling")
        name = os.path.basename(os.path.normpath(self._path))

        # The path of the output directory to be created
        dataoutputpath = os.path.join(os.getenv("VSC_DATA"), name, self._simulationname, "out_" + self._mode + "_" + str(processors))

        # Create a seperate output directory for this run (different runs can be executed simultaneously)
        self._log.info("The output of this run will be placed in " + dataoutputpath)
        os.makedirs(dataoutputpath)

        # The path of the log file for this simulation run
        logfilepath = os.path.join(dataoutputpath, self._skifilename + "_log.txt")

        # Get the timings from a serial run of the simulation
        seriallogfilepath = os.path.join(self._outpath, self._skifilename + "_log.txt")
        timings = extract(seriallogfilepath)
        serialtime = timings[0] + timings[2]     # setuptime + writingtime in seconds
        paralleltime = timings[1]                # in seconds

        # Set the expected walltime
        walltime = int((serialtime + paralleltime / processors)*1.5 + 100)  # in seconds

        # Create the job script
        jobscriptpath = os.path.join(self._outpath, "job_" + self._mode + "_" + str(processors) + ".sh")
        jobscript = JobScript(jobscriptpath, self._skifilepath, nodes, ppn, self._mode == "hybrid", dataoutputpath, walltime)

        # Add the command to go the the PTS do directory
        jobscript.addcommand("cd $VSC_HOME/PTS/git/do", comment="Navigate to the PTS do directory")

        # Add the PTS command to extract the timings of this run
        processes = processors if self._mode == "mpi" else nodes
        threads = 1 if self._mode == "mpi" else ppn

        command = "python extractscaling.py " + logfilepath + " " + str(processes) + " " + str(threads) + " " + resultsfilepath
        jobscript.addcommand(command, comment="Extract the results")

        # Add the command to remove the output directory of this run
        if not keepoutput:
            command = "cd; rm -rf " + dataoutputpath
            jobscript.addcommand(command, comment="Remove the temporary output directory")

        # Change the directory to the output directory for this simulation for the output.txt and error.txt files
        os.chdir(self._outpath)

        # Submit the job script to the cluster scheduler
        jobscript.submit()

        # Remove this job script (it has been submitted)
        if not keepoutput:
            jobscript.remove()

    ## This function runs the simulation once with the specified number of threads,
    #  and writes the timing results to the specified file object. This function takes the following arguments:
    #
    #  - processors: the number of processors to be used for this run
    #  - resultsfilepath: the path of the file that should contain the timings from this run
    #   - keepoutput: a flag indicating whether the SKIRT output should be kept or deleted
    #
    def _run(self, processors, resultsfilepath, keepoutput):

        # Calculate the number of processes and the number of threads
        processes = processors if self._mode == "mpi" else 1
        threads = 1 if self._mode == "mpi" else processors

        # Inform the user about the number of processes and threads used for this run
        self._log.info("Running simulation with " + str(processes) + " process(es) consisting of " + str(threads) + " thread(s)")

        # Run the simulation
        simulation = self._skirt.execute(skipattern=self._skifilepath, outpath=self._outpath, threads=threads, processes=processes)[0]

        # Check whether the simulation finished
        if simulation.status() != "Finished": raise ValueError("Simulation " + simulation.status())

        # Determine the path to the simulation log file
        logfilepath = os.path.join(self._outpath, self._simulationname + "_log.txt")

        # Extract the timings from the log file and place them in the results file
        extract(logfilepath, processes, threads, resultsfilepath)

        # Remove the contents of the output directory, if requested
        if not keepoutput:
            fileList = os.listdir(self._outpath)
            for filename in fileList:

                filepath=os.path.join(self._outpath,filename)

                if os.path.isfile(filepath):
                    os.remove(filepath)

    ## This function creates the file containing the results of the scaling benchmark test. It takes the maximum
    #  number of nodes, as given to the command line, as an argument. This number is used in the name of the results file,
    #  to identify this particular scaling test.
    def _createresultsfile(self, maxnodes, minnodes):

        # Generate a timestamp identifying this particular run for the ski file
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # Create a new file, whose name includes the system identifier, the scaling test mode, the maximum and
        # minium number of nodes and a timestamp.
        filepath = os.path.join(self._respath, self._system + "_" + self._mode + "_" + str(maxnodes)
                                + "_" + str(minnodes) + "_" + timestamp + ".dat")
        resultsfile = open(filepath, "w")

        # Write a header containing useful information about this test to the results file
        resultsfile.write("# Parallel scaling benchmark results for " + self._system + " in " + self._mode + " mode\n")
        if self._skirt is not None:
            resultsfile.write("# Using " + self._skirt.version() + "\n")
        resultsfile.write("# Column 1: Number of processes p\n")
        resultsfile.write("# Column 2: Number of threads per process t\n")
        resultsfile.write("# Column 3: Total number of threads (t*p)\n")
        resultsfile.write("# Column 4: Execution time for the setup (s)\n")
        resultsfile.write("# Column 5: Execution time for the stellar emission phase (s)\n")
        resultsfile.write("# Column 6: Execution time for the writing phase (s)\n")
        resultsfile.write("# Column 7: Execution time for the simulation (s)\n")

        # Close the results file (results will be appended!)
        resultsfile.close()

        # Return the path of the newly created results file
        return filepath

# -----------------------------------------------------------------