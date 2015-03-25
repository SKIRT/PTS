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
import subprocess
import multiprocessing
import datetime
import numpy as np
import shutil

# Import the relevant PTS class
from pts.skirtexec import SkirtExec
from pts.log import Log
from do.extractscaling import extract
from pts.jobscript import JobScript

# -----------------------------------------------------------------

# Ignore warnings, otherwise Canopy would give a UserWarning on top of the error encountered when a scaling
# test results file does not contain any data (an error which is catched an produces an error message).
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------

cores = {'delcatty': 16, 'gulpin': 32, 'cosma': 16}

# -----------------------------------------------------------------

## An instance of the ScalingTest class represents a SKIRT scaling benchmark test for a particular ski file.
#
class ScalingTest(object):

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

        # Create the logging mechanism
        self._log = Log()

        # Set the path
        self._path = path

        # Scaling benchmark name (="SKIRTscaling")
        self._name = os.path.basename(os.path.normpath(self._path))

        # Set the simulation name (subdirectory of self._path) and the simulation path
        self._simulationname = simulation
        self._simulationpath = os.path.join(self._path, self._simulationname)

        # Check whether the simulation path exists
        if not os.path.exists(self._simulationpath):
            self._log.error("No directory called " + self._simulationname + " was found in " + self._path)
            exit()

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
        if self._system in cores.keys():

            # Set the appropriate function for scheduling simulations
            self._do = self._schedule

            # Get the number of cores (per node) on this system from a pre-defined dictionary
            self._cores = cores[self._system]

            # Set the environment to launch jobs to the requested cluster
            clusterstring = "cluster/" + self._system
            subprocess.call(["module", "swap", clusterstring])

        else:

            # Set the appropriate function for running simulations
            self._do = self._run

            # Determine the number of cores (per node) on this system
            self._cores = multiprocessing.cpu_count()

        # Set the mode
        self._mode = mode

        # Search the ski file that is in the specified simulation path
        self._skifilename = ""
        for filename in os.listdir(self._simulationpath):
            if filename.endswith(".ski") and not filename.startswith("."):
                self._skifilename = filename[:-4]
                break

        # Check whether a ski file is found
        if not self._skifilename:
            self._log.error("No ski file was found in the simulation path " + self._simulationpath)
            exit()

        # The path of the ski file
        self._skifilepath = os.path.join(self._simulationpath, self._skifilename + ".ski")

        # Create skirt execution context
        self._skirt = SkirtExec(log=self._log)

    ## When this function is invoked, the scaling test is started. This function takes the following arguments:
    #
    #  - maxnodes: the maximum number of 'nodes' to be used for this scaling test. On a desktop system, a node is
    #              considered to be the entire set of processors. This number can range from zero to infinity, and
    #              can be a fractional number. To use 2 processors on a system with 4 processors (per node) for example,
    #              a value of 0.5 can be used for maxnodes.
    #  - minnodes: the minimum number of 'nodes' to be used for this scaling test. The usage is similar as with
    #              maxnodes. The default value of minnodes is zero, denoting no lower limit on the number of nodes.
    #              The minimum number of processors used for the scaling test will then equal one.
    #              In hybrid mode, minnodes also defines the number of parallel threads per process.
    #  - manual: a flag indicating whether the job script should be submitted from within this script, or just saved
    #            so that the user can inspect it and launch it manually.
    #  - keepoutput: this optional argument indicates whether the output of the SKIRT simulations has to be kept
    #                or can be deleted from the disk.
    #
    def run(self, maxnodes, minnodes, manual=False, keepoutput=False):

        # Log the system name, the test mode and the version of SKIRT used for this test
        self._log.info("Starting parallel scaling benchmark for " + self._system + " in " + self._mode + " mode.")
        self._log.info("Using " + self._skirt.version())

        # Generate a timestamp identifying this particular run for the scaling test
        self._timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        maxprocessors = int(maxnodes * self._cores)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one.
        minprocessors = int(minnodes * self._cores)
        if minprocessors == 0: minprocessors = 1
        processors = minprocessors

        # In hybrid mode, the minimum number of processors also represents the number of threads per process
        self._threadspp = 1
        if self._mode == "hybrid":
            self._threadspp = minprocessors

        # Create a file containing the results of the scaling test
        resultsfilepath = self._createresultsfile(maxnodes, minnodes)

        # Perform the simulations with increasing number of processors
        while processors <= maxprocessors:

            # Perform this run
            self._do(processors, resultsfilepath, manual, keepoutput)

            # The next run will be performed with double the amount of processors
            processors *= 2

        # End with some log messages
        self._log.success("Finished parallel scaling benchmark script")
        self._log.info("The results are / will be written to " + resultsfilepath)

    ## This functions schedules a simulation on the cluster. This function takes the following arguments:
    #
    #  - processors: the total number of processors to be used for this run
    #  - resultsfilepath: the path of the file that should contain the timings from this run
    #  - manual: a flag indicating whether the job script should be submitted from within this script, or just saved
    #            so that the user can inspect it and launch it manually
    #  - keepoutput: a flag indicating whether the SKIRT output should be kept or deleted
    #
    def _schedule(self, processors, resultsfilepath, manual, keepoutput):

        # Determine the number of processes and threads per process
        processes, threads = self._getmapping(processors)

        # Determine the number of nodes and processors per node
        nodes, ppn = self._getrequirements(processors)

        # In threads mode, show a warning message if the number of threads > the number of cores per node
        # (we can't use multiple nodes in threads mode)
        if self._mode == "threads" and threads > self._cores:

            # Show a warning and return immediately
            self._log.warning("The number of threads " + str(threads) + " exceeds the number of cores on this system: skipping")
            return

        # Inform the user about the number of processors, processes, threads per process, nodes and processors per node
        self._log.info("Scheduling simulation with:")
        self._log.info(" - total number of processors = " + str(processors))
        self._log.info(" - number of parallel processes = " + str(processes))
        self._log.info(" - number of parallel threads per process = " + str(threads))
        self._log.info(" - number of nodes = " + str(nodes))
        self._log.info(" - number of requested processors per node = " + str(ppn))

        # Check that we are indeed on the UGent HPC system: look for the VSC_DATA directory (we need it)
        vscdatapath = os.getenv("VSC_DATA")
        if vscdatapath is None:
            self._log.error("Can't find $VSC_DATA (Are you on indeed on the HPC infrastructure?)")
            exit()

        # All the output generated by SKIRT is placed in a directory within VSC_DATA, identified by the scaling
        # benchmark name (="SKIRTScaling") and the simulation name
        outputpath = os.path.join(vscdatapath, self._name, self._simulationname)

        # Create a seperate output directory for this run (different runs can be executed simultaneously)
        dataoutputpath = self._createdatadir(processors, outputpath)

        # The path of the log file for this simulation run
        logfilepath = os.path.join(dataoutputpath, self._skifilename + "_log.txt")

        # Calculate the expected walltime for this number of processors
        walltime = self._estimatewalltime(processors)

        # Create the job script. The name of the script indicates the mode in which we run this scaling test and
        # the current number of processors used. We enable the SKIRT verbose logging mode to be able to compare
        # the progress of the different parallel processes afterwards. Because for scaling tests, we don't want
        # processes to end up on different nods or the SKIRT processes sensing interference from other programs,
        # we set the 'fullnode' flag to True, which makes sure we always request at least one full node, even when
        # the current number of processors is less than the number of cores per node.
        jobscriptpath = os.path.join(self._outpath, "job_" + self._mode + "_" + str(processors) + ".sh")
        jobscript = JobScript(jobscriptpath, self._skifilepath, nodes, ppn, threads, dataoutputpath, walltime, verbose=True, fullnode=True)

        # Add the command to go the the PTS do directory
        jobscript.addcommand("cd $VSC_HOME/PTS/git/do", comment="Navigate to the PTS do directory")

        # Add the PTS command to extract the timings of this run
        command = "python extractscaling.py " + logfilepath + " " + str(processes) + " " + str(threads) + " " + resultsfilepath
        jobscript.addcommand(command, comment="Extract the results")

        # Add the command to remove the output directory of this run
        if not keepoutput:
            command = "cd; rm -rf " + dataoutputpath
            jobscript.addcommand(command, comment="Remove the temporary output directory")

        # Change the directory to the output directory for this simulation for the output.txt and error.txt files
        os.chdir(self._outpath)

        # Submit the job script to the cluster scheduler
        if not manual: jobscript.submit()

        # Remove this job script (it has been submitted)
        if not manual and not keepoutput: jobscript.remove()

    ## This function runs the simulation once with the specified number of threads,
    #  and writes the timing results to the specified file object. This function takes the following arguments:
    #
    #  - processors: the number of processors to be used for this run
    #  - resultsfilepath: the path of the file that should contain the timings from this run
    #  - manual: a flag indicating whether the job script should be submitted from within this script, or just saved
    #            so that the user can inspect it and launch it manually
    #  - keepoutput: a flag indicating whether the SKIRT output should be kept or deleted
    #
    def _run(self, processors, resultsfilepath, manual, keepoutput):

        # Determine the number of processes and threads per process
        processes, threads = self._getmapping(processors)

        # Inform the user about the number of processes and threads used for this run
        self._log.info("Running simulation with " + str(processes) + " process(es), each consisting of " + str(threads) + " thread(s)")

        # Create a seperate output directory for this run (different runs can be executed simultaneously)
        dataoutputpath = self._createdatadir(processors)

        # Run the simulation
        simulation = self._skirt.execute(skipattern=self._skifilepath, outpath=dataoutputpath, threads=threads, processes=processes, verbose=True)[0]

        # Check whether the simulation finished
        if simulation.status() != "Finished": raise ValueError("Simulation " + simulation.status())

        # Extract the timings from the simulation's log file and place them in the results file
        extract(simulation.logfilepath(), processes, threads, resultsfilepath)

        # Remove the contents of the output directory, if requested
        if not keepoutput: shutil.rmtree(dataoutputpath)

    ## This function creates the file containing the results of the scaling benchmark test. It takes the maximum
    #  and minimum number of nodes, as given to the command line, as arguments. These numbers are used in the name
    #  of the results file, to identify this particular scaling test.
    def _createresultsfile(self, maxnodes, minnodes):

        # If hybrid mode is selected, add the number of threads per process to the name of the results file
        hybridinfo = ""
        if self._mode == "hybrid":
            hybridinfo = str(self._threadspp)

        # Create a new file, whose name includes the system identifier, the scaling test mode, the maximum and
        # minium number of nodes and a timestamp.
        filepath = os.path.join(self._respath, self._system + "_" + self._mode + hybridinfo + "_" + str(maxnodes)
                                + "_" + str(minnodes) + "_" + self._timestamp + ".dat")
        resultsfile = open(filepath, "w")

        # Write a header containing useful information about this test to the results file
        resultsfile.write("# Parallel scaling benchmark results for " + self._system + " in " + self._mode + " mode\n")
        resultsfile.write("# Using " + self._skirt.version() + "\n")
        resultsfile.write("# Column 1: Number of processes p\n")
        resultsfile.write("# Column 2: Number of threads per process t\n")
        resultsfile.write("# Column 3: Total number of threads (t*p)\n")
        resultsfile.write("# Column 4: Execution time for the setup (s)\n")
        resultsfile.write("# Column 5: Execution time for the stellar emission phase (s)\n")
        resultsfile.write("# Column 6: Execution time for the dust self-absorption phase (s)\n")
        resultsfile.write("# Column 7: Execution time for the dust emission phase (s)\n")
        resultsfile.write("# Column 8: Execution time for the writing phase (s)\n")
        resultsfile.write("# Column 9: Execution time for the simulation (s)\n")

        # Close the results file (results will be appended!)
        resultsfile.close()

        # Return the path of the newly created results file
        return filepath

    ## This function calculates the number of processes and the number of threads (per process) for
    #  a certain number of processors, depending on the mode in which this scaling test is run.
    #  In other words, this function determines the 'mapping' from a set of processors to an appropriate
    #  set of threads and processes. This function takes the number of processors as the sole argument.
    def _getmapping(self, processors):

        threads = 1
        processes = 1
        if self._mode == "mpi":

            # In mpi mode, each processor runs a different process
            processes = processors

        if self._mode == "threads":

            # In threads mode, each processor runs a seperate thread within the same process
            threads = processors

        if self._mode == "hybrid":

            # In hybrid mode, the number of processes depends on how many threads are requested per process
            # and the current number of processors
            threads = self._threadspp
            processes = processors / self._threadspp

        return processes, threads

    ## This function calculates the required amount of nodes and processors per node, given a certain number
    #  of processors. This function is only used when on a cluster with a scheduling system.
    def _getrequirements(self, processors):

        # Calculate the necessary amount of nodes
        nodes = processors/self._cores + (processors % self._cores > 0)

        # Determine the number of processors per node
        ppn = processors if nodes == 1 else self._cores

        # Return the number of nodes and processors per node
        return nodes, ppn

    ## This function creates a directory to contain the output of a certain run during the scaling test.
    #  The name of the directory includes the mode in which the scaling test was run, the used number of
    #  processors and the timestamp identifying this test. This function takes the following arguments:
    #
    #  - processors: the number of processors for this run of the scaling test
    #  - outputpath: this optional argument determines the directory where the results from each scaling test will
    #                be placed (in a seperate directory identifying the mode, number of processors and a timestamp.
    #
    def _createdatadir(self, processors, outputpath=None):

        # If no outputpath is given, use the default "out" directory next to the ski file
        if outputpath is None: outputpath = self._outpath

        # The path of the output directory to be created
        dataoutputpath = os.path.join(outputpath, "out_" + self._mode + "_" + str(processors) + "_" + self._timestamp)

        # Create a seperate output directory for this run (different runs can be executed simultaneously)
        self._log.info("The output of this run will be placed in " + dataoutputpath)
        os.makedirs(dataoutputpath)

        # Return the path to the output directory
        return dataoutputpath

    ## This function estimates the total runtime (or walltime) for the current simulation, number of processors,
    #  system and scaling test mode. This function takes the following arguments:
    #
    #  - processors: the number of processors (or total threads) used for this run of the simulation
    #  - factor: this optional argument determines how much the upper limit on the walltime should deviate from
    #            a previous run of this simulation.
    #
    def _estimatewalltime(self, processors, factor=1.5):

        # Try to get the runtimes for this number of processors
        runtimes = self._getruntimes(processors)

        # If these runtimes could be found, use the total runtime times a surplus of 1.5 as an upper limit to the
        # walltime for this run
        if runtimes is not None:

            # Return the estimated walltime (as an integer number in seconds)
            return int(runtimes["total"]*factor)

        # If runtimes for this number of processors could not be found, we estimate the walltime by getting the
        # runtimes from a serial run of the simulation
        else:

            # Get the runtimes for this simulation, ran on 1 processor (with one thread), and don't look at
            # the system name or scaling test mode.
            runtimes = self._getruntimes(1, anysystem=True, anymode=True)

        # Check if finding the serial runtimes was successfull or not
        if runtimes is not None:

            # Calculate the portion of the total runtime spent in serial and parallel parts of the code
            serialtime = runtimes['setup'] + runtimes['writing']
            paralleltime = runtimes['stellar'] + runtimes['dustselfabs'] + runtimes['dustem']

            # Estimate the total runtime for this number of processors, by taking an overhead of 1 percent per
            # parallel process
            totaltime = (serialtime + paralleltime / processors) * (1.0 + 0.01*processors)

            # Calculate and return the expected walltime (as an integer number in seconds)
            return int(totaltime*factor)

        # As a last resort, look for a log file that was placed next to the ski file of this scaling test
        else:

            # The path of the serial log file
            seriallogfilepath = os.path.join(self._simulationpath, self._skifilename + "_log.txt")

            # Check whether such a file exists
            if os.path.exists(seriallogfilepath):

                # TODO: check whether the log file comes from a simulation with 1 process and 1 thread

                # If such a log file is present, extract the timings from it
                runtimes = extract(seriallogfilepath)

                # Calculate the portion of the total runtime spent in serial and parallel parts of the code
                serialtime = runtimes['setup'] + runtimes['writing']
                paralleltime = runtimes['stellar'] + runtimes['dustselfabs'] + runtimes['dustem']

                # Estimate the total runtime for this number of processors, by taking an overhead of 1 percent per
                # parallel process
                totaltime = (serialtime + paralleltime / processors) * (1.0 + 0.01*processors)

                # Calculate and return the expected walltime (as an integer number in seconds)
                return int(totaltime*factor)

            # If not, exit with an error
            else:

                self._log.error("The walltime could not be estimated for a run with " + str(processors) + " processors.")
                exit()

    ## This function extracts the timings for the current simulation from a scaling test results file that was
    #  created earlier for this simulation. This function takes the following arguments:
    #
    #  - processors: the number of processors (or total threads) for which to look up the runtimes
    #  - anysystem: this flag tells whether we must look for results of any system, or only from the system
    #               we are currently running on
    #   - anymode: this flag tells whether we must look for results created with any mode, or only with the mode
    #              in which we are running the current scaling test
    #
    def _getruntimes(self, processors, anysystem=False, anymode=False):

        # Search for a results file corresponding with this system and the current mode
        for filename in os.listdir(self._respath):

            # Check whether this file is a data file and not hidden
            if not filename.endswith(".dat") or filename.startswith("."): continue

            # Determine the full path to this results file
            filepath = os.path.join(self._respath, filename)

            # Split the file name into its segments
            segments = filename.split("_")

            # Get the system name in which the scaling test was run for this results file
            systemname = segments[0]

            # Get the mode in which the scaling test was run for this results file
            mode = segments[1]

            # Check whether this results file corresponds to a test on this system and the current mode
            if (anysystem or systemname == self._system) and (anymode or mode.startswith(self._mode)):

                # Try extracting the columns from the data file
                try:
                    threads, setup, stellar, dustselfabs, dustem, writing, total = np.loadtxt(filepath, usecols=(2,3,4,5,6,7,8), unpack=True)
                except (IndexError, ValueError):
                    # Try the next file, this one is probably empty
                    continue

                # Look for an entry corresponding to the current number of processors (the total number of threads)
                try:
                    index = [int(nthreads) for nthreads in threads].index(processors)

                    # Create a dictionary specifying the serial runtime of each of the different simulation phases
                    runtimes = dict()
                    runtimes['setup'] = setup[index]
                    runtimes['stellar'] = stellar[index]
                    runtimes['dustselfabs'] = dustselfabs[index]
                    runtimes['dustem'] = dustem[index]
                    runtimes['writing'] = writing[index]
                    runtimes['total'] = total[index]

                    # Return the runtimes
                    return runtimes

                except ValueError:
                    # Try the next file, no entry for a serial run could be found in this one
                    pass

        # Return None if no file was found for this system and/or mode, or with the timings for this number of processors
        return None

# -----------------------------------------------------------------
