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
import numpy as np
import shutil

# Import the relevant PTS class
from pts.skifile import SkiFile
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

        # Set the input, output, result and visualisation paths
        self._inpath = os.path.join(self._simulationpath, "in")
        self._outpath = os.path.join(self._simulationpath, "out")
        self._respath = os.path.join(self._simulationpath, "res")
        self._vispath = os.path.join(self._simulationpath, "vis")

        # Create the output directories if they do not already exist
        try: os.mkdir(self._outpath)
        except OSError: pass
        try: os.mkdir(self._respath)
        except OSError: pass
        
        # Set the system name
        self._system = system

        # Determine whether we are dealing with a scheduling system or we can launch simulations right away
        if self._system in cores.keys():

            # Indicate that we are working with a scheduling system
            self._scheduler = True

            # Get the number of cores (per node) on this system from a pre-defined dictionary
            self._cores = cores[self._system]

        else:

            # Indicate that we are not dealing with a scheduling system
            self._scheduler = False

            # Determine the number of cores (per node) on this system
            self._cores = multiprocessing.cpu_count()

        # Set the dataoutputpath
        if self._scheduler:

            # Check that we are indeed on the UGent HPC system: look for the VSC_SCRATCH or alternatively the VSC_DATA
            # directory (we need it to store the SKIRT output)
            vscdatapath = os.getenv("VSC_SCRATCH_" + self._system.upper())

            # If no VSC_SCRATCH directory could be found for this cluster (we are not using delcatty or gulpin), then use
            # the VSC_DATA directory
            if vscdatapath is None: vscdatapath = os.getenv("VSC_DATA")

            # If VSC_DATA is not found, we are not on the HPC infrastructure
            if vscdatapath is None:

                self._log.error("Can't find $VSC_DATA (Are you on indeed on the HPC infrastructure?)")
                exit()

            # All the output generated by SKIRT is placed in a directory within VSC_SCRATCH or VSC_DATA, identified by
            # the scaling benchmark name (="SKIRTScaling") and the simulation name
            self._dataoutputpath = os.path.join(vscdatapath, self._name, self._simulationname)

        else:

            self._dataoutputpath = self._outpath

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
    #  - extractprogress: a flag indicating whether the progress of the different SKIRT processes should be extracted
    #                     from the simulation's log files or not
    #  - extracttimeline: a flag indicating whether timeline information of the different SKIRT processes should be
    #                     extracted from the simulation's log files or not
    #  - weak: a flag indicating whether a weak or a strong scaling test should be performed
    #
    def run(self, maxnodes, minnodes, manual=False, keepoutput=False, extractprogr=False, extracttimeline=False, weak=False):

        # Set the properties of the scaling test
        self._manual = manual
        self._keepoutput = keepoutput
        self._extractprogr = extractprogr
        self._extracttimeline = extracttimeline
        self._weak = weak

        # Set a string specifying whether this scaling test is 'strong' or 'weak'
        scalingtype = "weak" if self._weak else "strong"

        # Log the system name, the test mode and the version of SKIRT used for this test
        self._log.info("Starting " + scalingtype + " scaling benchmark for " + self._system + " in " + self._mode + " mode")
        self._log.info("Using " + self._skirt.version())

        # Generate a timestamp identifying this particular run for the scaling test
        self._timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        maxprocessors = int(maxnodes * self._cores)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one
        minprocessors = int(minnodes * self._cores)
        if minprocessors == 0: minprocessors = 1
        processors = minprocessors

        # In hybrid mode, the minimum number of processors also represents the number of threads per process
        self._threadspp = 1
        if self._mode == "hybrid": self._threadspp = minprocessors

        # If hybrid mode is selected, add the number of threads per process to the name of the results directory
        hybridinfo = str(self._threadspp) if self._mode == "hybrid" else ""

        # Define a name identifying this scaling test run
        self._scalingrunname = self._system + "_" + self._mode + hybridinfo + "_" + str(maxnodes) + "_" + str(minnodes) + "_" + self._timestamp

        # Create a directory that will contain the output of the simulations that are part of this scaling test run
        self._runoutputdir = self._createrunoutputdir()

        # Create a directory to contain the results (scaling, progress and timeline data) for this scaling test
        self._resultsdirpath = self._createresultsdir()

        # Inside this directory, create a file named 'scaling.dat' to contain the runtimes from which the scaling
        # behaviour can be inferred
        self._scalingfilepath = self._createscalingfile()

        # Create a file which gives useful information about this scaling test run
        self._infofilepath = self._createinfofile(maxnodes, minnodes)

        # Perform the simulations with increasing number of processors
        while processors <= maxprocessors:

            # Perform this run
            self._do(processors)

            # The next run will be performed with double the amount of processors
            processors *= 2

        # End with some log messages
        self._log.success("Finished parallel scaling benchmark script")
        self._log.info("The results are / will be written to " + self._scalingfilepath)

    ## This function performs one run of the scaling test for a particular number of processors
    def _do(self, processors):

        # Open the info file
        infofile = open(self._infofilepath, 'a')

        # Determine the number of processes and threads per process
        processes, threads = self._getmapping(processors)

        # Set the path of the ski file used for this run
        skifilepath = self._skifilepath

        # Create a seperate output directory for this simulation (different simulations can be executed simultaneously)
        simulationoutputpath = self._createskirtoutputdir(processors)

        # If a 'weak' scaling test is performed, create a ski file that is adjusted to the current number of processors
        if self._weak and processors > 1:

            # Open the original ski file and adjust the number of photon packages and the number of dust cells
            skifile = SkiFile(self._skifilepath)
            skifile.increasepackages(processors)
            skifile.increasedustcells(processors)

            # Save the adjusted ski file in the output path for this simulation
            skifilepath = os.path.join(simulationoutputpath, self._skifilename + ".ski")
            skifile.saveto(skifilepath)

        # Write some information about this simulation to the info file
        infofile.write("Simulation performed on " + str(processors) + " processors\n")
        infofile.write(" - ski file: " + skifilepath + "\n")
        infofile.write(" - output directory: " + simulationoutputpath + "\n")
        infofile.write(" - number of processes: " + str(processes) + "\n")
        infofile.write(" - number of threads per processes: " + str(threads) + "\n")

        # Schedule or launch the simulation
        if self._scheduler: self._schedule(processors, processes, threads, skifilepath, simulationoutputpath, infofile)
        else: self._launch(processors, processes, threads, skifilepath, simulationoutputpath, infofile)

        # Close the info file
        infofile.write("\n")
        infofile.close()

    ## This functions schedules a simulation on the cluster. This function takes the following arguments:
    #
    #  - processors: the total number of processors to be used for this run
    #
    def _schedule(self, processors, processes, threads, skifilepath, dataoutputpath, infofile):

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

        # Write the number of nodes and processors per node to the info file
        infofile.write(" - number of used nodes: " + str(nodes) + "\n")
        infofile.write(" - number of requested processors per node: " + str(ppn) + "\n")

        # The path of the log file for this simulation run
        logfilepath = os.path.join(dataoutputpath, self._skifilename + "_log.txt")

        # Calculate the expected walltime for this number of processors
        walltime = self._estimatewalltime(processors)

        # Create the job script. The name of the script indicates the mode in which we run this scaling test and
        # the current number of processors used. We enable the SKIRT verbose logging mode to be able to compare
        # the progress of the different parallel processes afterwards. Because for scaling tests, we don't want
        # processes to end up on different nodes or the SKIRT processes sensing interference from other programs,
        # we set the 'fullnode' flag to True, which makes sure we always request at least one full node, even when
        # the current number of processors is less than the number of cores per node.
        jobscriptpath = os.path.join(self._outpath, "job_" + self._mode + "_" + str(processors) + ".sh")
        jobscript = JobScript(jobscriptpath, skifilepath, self._system, nodes, ppn, threads, self._inpath, dataoutputpath, walltime, brief=True, verbose=True, fullnode=True)

        # Add the command to go the the PTS do directory
        jobscript.addcommand("cd $VSC_HOME/PTS/git/do", comment="Navigate to the PTS do directory")

        # Add the PTS command to extract the timings of this run
        command = "python extractscaling.py " + logfilepath + " " + str(processes) + " " + str(threads) + " " + self._scalingfilepath
        jobscript.addcommand(command, comment="Extract the results")

        # Add the command to extract the progress information after the job finished, if requested
        if self._extractprogr and processes > 1:

            # Create the file to contain the progress information of this run
            progressfilepath = self._createprogressfile(processes)

            # Add the command to the jobscript to extract the progress
            command = "python extractprogress.py " + self._skifilename + " " + dataoutputpath + " " + progressfilepath
            jobscript.addcommand(command, comment="Extract the progress of the different processes")

            # Write the path of the progress file to the info file
            infofile.write(" - progress information extracted to: " + progressfilepath + "\n")

        # Add the command to extract the timeline information after the job finished, if requested
        if self._extracttimeline:

            # Create the file to contain the timeline information for this run
            timelinefilepath = self._createtimelinefile(processes)

            # Add the command to the jobscript to extract the timeline data
            command = "python extracttimeline.py " + self._skifilename + " " + dataoutputpath + " " + timelinefilepath
            jobscript.addcommand(command, comment="Extract the timeline information for the different processes")

            # Write the path of the timeline file to the info file
            infofile.write(" - timeline information extracted to: " + timelinefilepath + "\n")

        # Add the command to remove the output directory of this run
        if not self._keepoutput:

            # Add the command to the jobscript
            command = "cd; rm -rf " + dataoutputpath
            jobscript.addcommand(command, comment="Remove the temporary output directory")

        # Change the directory to the output directory for this simulation for the 'output.txt' and 'error.txt' files
        os.chdir(self._outpath)

        # Submit the job script to the cluster scheduler
        if not self._manual: jobscript.submit()

        # Remove this job script (it has been submitted)
        if not self._manual and not self._keepoutput: jobscript.remove()

    ## This function launches the simulation once with the specified number of threads,
    #  and writes the timing results to the specified file object. This function takes the following arguments:
    #
    #  - processors: the number of processors to be used for this run
    #
    def _launch(self, processors, processes, threads, skifilepath, dataoutputpath, infofile):

        # Inform the user about the number of processes and threads used for this run
        self._log.info("Launching simulation with " + str(processes) + " process(es), each consisting of " + str(threads) + " thread(s)")

        # Run the simulation
        simulation = self._skirt.execute(skipattern=skifilepath, inpath=self._inpath, outpath=dataoutputpath, threads=threads, processes=processes, brief=True, verbose=True)[0]

        # Check whether the simulation finished
        if simulation.status() != "Finished": raise ValueError("Simulation " + simulation.status())

        # Extract the timings from the simulation's log file and place them in the results file
        extract(simulation.logfilepath(), processes, threads, self._scalingfilepath)

        # Extract the progress of the different processes, if requested
        if self._extractprogr and processes > 1:

            # Create the file to contain the progress information of this run
            progressfilepath = self._createprogressfile(processes)

            # Load the extractprogress module
            import do.extractprogress

            # Extract the progress information
            do.extractprogress.extract(self._skifilename, dataoutputpath, progressfilepath)

            # Write the path of the progress file to the info file
            infofile.write(" - progress information extracted to: " + progressfilepath + "\n")

        # Add the command to extract the timeline information after the job finished, if requested
        if self._extracttimeline:

            # Create the file to contain the timeline information for this run
            timelinefilepath = self._createtimelinefile(processes)

            # Load the extracttimeline module
            import do.extracttimeline

            # Extract the timeline information
            do.extracttimeline.extract(self._skifilename, dataoutputpath, timelinefilepath)

            # Write the path of the timeline file to the info file
            infofile.write(" - timeline information extracted to: " + timelinefilepath + "\n")

        # Remove the contents of the output directory, if requested
        if not self._keepoutput: shutil.rmtree(dataoutputpath)

    ## This function creates the directory containing the results of the scaling benchmark test and creates the file
    #  'scaling.dat' inside this directory that will contain the runtimes of the test. This function takes the
    #  maximum and minimum number of nodes, as given to the command line, as arguments. These numbers are used in
    #  the name of the results directory, to identify this particular scaling test.
    def _createresultsdir(self):

        # Create a new directory, whose name includes the system identifier, the scaling test mode, the maximum
        # and minium number of nodes and a timestamp
        resultsdirpath = os.path.join(self._respath, self._scalingrunname)
        os.mkdir(resultsdirpath)

        # Return the path to the newly created directory
        return resultsdirpath

    ## This function creates the directory containing the output of this scaling test run
    def _createrunoutputdir(self):

        # Create the new directory
        runoutputpath = os.path.join(self._dataoutputpath, self._scalingrunname)
        os.makedirs(runoutputpath)

        # Return the path to the run output directory
        return runoutputpath

    ## This function creates a directory to contain the output of a certain simulation during a scaling test run.
    #  The name of the directory includes the mode in which the scaling test was run and the used number of
    #  processors. This function takes the following arguments:
    #
    #  - processors: the number of processors for this run of the scaling test
    #
    def _createskirtoutputdir(self, processors):

        # The path of the output directory to be created
        skirtoutputpath = os.path.join(self._runoutputdir, "out_" + self._mode + "_" + str(processors))

        # Create a seperate output directory for this simulation (different simulations can be executed simultaneously)
        self._log.info("The output of this simulation will be placed in " + skirtoutputpath)
        os.makedirs(skirtoutputpath)

        # Return the path to the SKIRT output directory
        return skirtoutputpath

    ## This function creates a file containing general information about the current scaling test run
    def _createinfofile(self, maxnodes, minnodes):

        # Create the file
        infofilepath = os.path.join(self._resultsdirpath, "info.txt")
        infofile = open(infofilepath, "w")

        # Set a string specifying whether this scaling test is 'strong' or 'weak'
        scalingtype = "weak" if self._weak else "strong"

        # If hybrid mode is selected, add the number of threads per process to the name of the results directory
        hybridinfo = " with " + str(self._threadspp) + " threads per process" if self._mode == "hybrid" else ""

        # Write some useful information to the file
        infofile.write("Scaling benchmark test " + self._scalingrunname + "\n")
        infofile.write("Scaling type: " + scalingtype + "\n")
        infofile.write("System: " + self._system + "\n")
        infofile.write("SKIRT version: " + self._skirt.version() + "\n")
        infofile.write("Mode: " + self._mode + hybridinfo + "\n")
        infofile.write("Maximum number of nodes: " + str(maxnodes) + "\n")
        infofile.write("Minimum number of nodes: " + str(minnodes) + "\n")
        infofile.write("\n")

        # Close the info file (information on specific simulations will be appended)
        infofile.close()

        # Return the path of the info file
        return infofilepath

    ## This function creates the file containing the scaling information
    def _createscalingfile(self):

        # Create the file
        scalingfilepath = os.path.join(self._resultsdirpath, "scaling.dat")
        scalingfile = open(scalingfilepath, "w")

        # Write a header to this new file which contains some general info about its contents
        scalingfile.write("# Timing results for scaling benchmark test " + self._scalingrunname + "\n")
        scalingfile.write("# Column 1: Number of processes p\n")
        scalingfile.write("# Column 2: Number of threads per process t\n")
        scalingfile.write("# Column 3: Total number of threads (t*p)\n")
        scalingfile.write("# Column 4: Execution time for the setup (s)\n")
        scalingfile.write("# Column 5: Execution time for the stellar emission phase (s)\n")
        scalingfile.write("# Column 6: Execution time for the dust self-absorption phase (s)\n")
        scalingfile.write("# Column 7: Execution time for the dust emission phase (s)\n")
        scalingfile.write("# Column 8: Execution time for the writing phase (s)\n")
        scalingfile.write("# Column 9: Execution time for the simulation (s)\n")

        # Close the scaling results file (results will be appended)
        scalingfile.close()

        # Return the path of the scaling results file
        return scalingfilepath

    ## This function creates the file containing the progress information
    def _createprogressfile(self, processes):

        # Create a new file whose name includes the current number of processors
        filepath = os.path.join(self._resultsdirpath, "progress_" + str(processes) + ".dat")
        progressfile = open(filepath, 'w')

        # Write a header to this new file which contains some general info about its contents
        progressfile.write("# Progress results for scaling benchmark test " + self._scalingrunname + "\n")
        progressfile.write("# Column 1: Simulation phase (0=stellar, 1=spectra, 2=dust)\n")
        progressfile.write("# Column 2: Process rank\n")
        progressfile.write("# Column 3: Execution time (s)\n")
        progressfile.write("# Column 4: Progress (%)\n")

        # Close the progress file (results will be appended)
        progressfile.close()

        # Return the path of the newly created progress file
        return filepath

    ## This function creates the file containing the timeline information
    def _createtimelinefile(self, processes):

        # Create a new file whose name includes the current number of processors
        filepath = os.path.join(self._resultsdirpath, "timeline_" + str(processes) + ".dat")
        timelinefile = open(filepath, 'w')

        # Write a header to this new file which contains some general info about its contents
        timelinefile.write("# Timeline results for scaling benchmark test " + self._scalingrunname + "\n")

        # Close the timeline file (results will be appended)
        timelinefile.close()

        # Return the path of the newly created timeline data file
        return filepath

    ## This function calculates the number of processes and the number of threads (per process) for
    #  a certain number of processors, depending on the mode in which this scaling test is run.
    #  In other words, this function determines the 'mapping' from a set of processors to an appropriate
    #  set of threads and processes. This function takes the number of processors as the sole argument.
    def _getmapping(self, processors):

        # Set default values for the number of threads and processes
        threads = 1
        processes = 1

        # In mpi mode, each processor runs a different process
        if self._mode == "mpi":

            processes = processors

        # In threads mode, each processor runs a seperate thread within the same process
        if self._mode == "threads":

            threads = processors

        # In hybrid mode, the number of processes depends on how many threads are requested per process
        # and the current number of processors
        if self._mode == "hybrid":

            threads = self._threadspp
            processes = processors / self._threadspp

        # Return the number of processes and the number of threads
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

    ## This function estimates the total runtime (or walltime) for the current simulation, number of processors,
    #  system and scaling test mode. This function takes the following arguments:
    #
    #  - processors: the number of processors (or total threads) used for this run of the simulation
    #  - factor: this optional argument determines how much the upper limit on the walltime should deviate from
    #            a previous run of this simulation.
    #
    def _estimatewalltime(self, processors, factor=1.5):

        # If the scaling test is weak, use the total runtime on one processor as the reference
        if self._weak:

            # Get the runtime for one processor, add a surplus of one percent per processor and apply the extra factor
            runtimes = self._getruntimes(1)

            if runtimes is not None:

                totaltime = runtimes["total"] * (1.0 + 0.01*processors)
                return int(totaltime*factor)

            else:

                # The path of the log file
                logfilepath = os.path.join(self._simulationpath, self._skifilename + "_log.txt")

                # TODO: don't assume the logfile was created by a simulation with only 1 process!

                # Check whether such a file exists
                if os.path.exists(logfilepath):

                    # If such a log file is present, extract the timings from it
                    runtimes = extract(logfilepath)

                    totaltime = runtimes["total"] * (1.0 + 0.01*processors)
                    return int(totaltime*factor)

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

            # The path of the log file
            logfilepath = os.path.join(self._simulationpath, self._skifilename + "_log.txt")

            # Check whether such a file exists
            if os.path.exists(logfilepath):

                # Initially, set the number of processes and threads from the log file to one
                logfileprocesses = 1
                logfilethreads = 1

                # Check with how many processes and threads this simulation was run, by reading each line of the
                # specified log file and searching for indications of multiple processes and/or multiple threads
                for line in open(logfilepath):

                    if 'Starting simulation ' + self._skifilename + ' with' in line:

                        logfileprocesses = int(line.split(' with ')[1].split()[0])

                    elif 'Initializing random number generator for thread number' in line:

                        # The last such line that is found states the rank of the last (highest-ranked) thread
                        logfilethreads = int(line.split(' for thread number ')[1].split()[0]) + 1

                # Calculate the total number of used processors used to create the log file
                logfileprocessors = logfileprocesses * logfilethreads

                # If such a log file is present, extract the timings from it
                runtimes = extract(logfilepath)

                # Calculate the portion of the total runtime spent in serial and parallel parts of the code
                serialtime = runtimes['setup'] + runtimes['writing']
                paralleltime = runtimes['stellar']*logfileprocessors \
                             + runtimes['dustselfabs']*logfileprocessors \
                             + runtimes['dustem']*logfileprocessors

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
