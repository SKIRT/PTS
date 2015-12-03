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

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import the relevant PTS classes and modules
from ..simulation.arguments import SkirtArguments
from .scalinganalyser import ScalingAnalyser
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..tools import time

# -----------------------------------------------------------------

class ScalingTest(Configurable):

    """
    An instance of the ScalingTest class represents a SKIRT scaling benchmark test for a particular ski file.
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ScalingTest, self).__init__(config)

        ## Attributes

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        # Create a ScalingAnalyser instance
        self.analyser = ScalingAnalyser()

        # Initialize a list to contain the retreived simulations
        self.simulations = []

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ScalingTest instance
        test = cls()

        ## Adjust the configuration settings according to the command-line arguments

        # Logging
        if arguments.debug: test.config.logging.level = "DEBUG"

        # Ski file
        test.config.ski_path = arguments.filepath

        # Remote host and cluster (if applicable)
        test.config.remote = arguments.remote
        test.config.cluster = arguments.cluster

        # Parallelization mode
        test.config.mode = arguments.mode

        # Maximum and minimum number of nodes
        test.config.min_nodes = arguments.minnodes
        test.config.max_nodes = arguments.maxnodes

        # Other options
        test.config.manual = arguments.manual
        test.config.keep_output = arguments.keep

        # Extraction
        test.config.extraction.progress = arguments.progress
        test.config.extraction.timeline = arguments.timeline
        test.config.extraction.memory = arguments.memory

        # Plotting
        test.config.plotting.progress = arguments.progress
        test.config.plotting.timeline = arguments.timeline
        test.config.plotting.memory = arguments.memory

        # Return the new scaling test
        return test

    # -----------------------------------------------------------------

    def run(self):

        """
        When this function is called, the scaling test is started.
        """

        # 1. Call the setup function
        self.setup()

        # 2. Perform the simulations as a part of the scaling test
        self.perform()

        # 3. Retreive the output of finished simulations (of other scaling test runs)
        self.retreive()

        # 3. Analyse the output of the retreived simulations
        self.analyse()

        # End with some log messages
        self.log.success("Finished scaling test run")
        self.log.info("The results are / will be written to " + self.scaling_file_path)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ScalingTest, self).setup()

        # Set the input, output, result, plot and temp paths
        self.base_path = os.path.dirname(self.config.ski_path) if "/" in self.config.ski_path else os.getcwd()
        self.input_path = os.path.join(self.base_path, "in")
        self.output_path = os.path.join(self.base_path, "out")
        self.result_path = os.path.join(self.base_path, "res")
        self.plot_path = os.path.join(self.base_path, "plot")
        self.temp_path = os.path.join(self.base_path, "temp")

        # Check if an input directory exists
        if not os.path.isdir(self.input_path): self.input_path = None

        # Create the output, result, plot and temp directories if necessary
        if not os.path.isdir(self.output_path): os.makedirs(self.output_path)
        if not os.path.isdir(self.result_path): os.makedirs(self.result_path)
        if not os.path.isdir(self.plot_path): os.makedirs(self.plot_path)
        if not os.path.isdir(self.temp_path): os.makedirs(self.temp_path)

        # Setup the remote execution context
        self.remote.setup(self.config.remote, self.config.cluster)

        # Determine whether we are dealing with a scheduling system or we can launch simulations right away
        self.scheduler = self.remote.config.scheduler

        # Get the number of cores (per node) on this system from a pre-defined dictionary
        self.cores = self.remote.cores

        # Determine the simulation prefix
        self.prefix = os.path.basename(self.config.ski_path).split(".")[0]

        ## Set the minimum and maximum number of processors and the number of threads per process (for hybrid mode)

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        self.max_processors = int(self.config.max_nodes * self.cores)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one
        self.min_processors = int(self.config.min_nodes * self.cores)
        if self.min_processors == 0: self.min_processors = 1

        # In hybrid mode, the minimum number of processors also represents the number of threads per process
        self.threads_per_process = 1
        if self.config.mode == "hybrid": self.threads_per_process = self.min_processors

        # If hybrid mode is selected, add the number of threads per process to the name of the results directory
        hybridinfo = str(self.threads_per_process) if self.config.mode == "hybrid" else ""

        ## Names, directories

        # Define a name identifying this scaling test run
        self.scaling_run_name = time.unique_name(self.config.remote + "_" + self.config.mode + hybridinfo + "_" + str(self.config.max_nodes) + "_" + str(self.config.min_nodes))

        # Determine the paths to the directories that will contain the output, results, plots and temporary files of this particular scaling test run
        self.output_path_run = os.path.join(self.output_path, self.scaling_run_name)
        self.result_path_run = os.path.join(self.result_path, self.scaling_run_name)
        self.plot_path_run = os.path.join(self.plot_path, self.scaling_run_name)
        self.temp_path_run = os.path.join(self.temp_path, self.scaling_run_name)

        # Create the output, result, plot and temp directories for this run if necessary
        if not os.path.isdir(self.output_path_run): os.makedirs(self.output_path_run)
        if not os.path.isdir(self.result_path_run): os.makedirs(self.result_path_run)
        if not os.path.isdir(self.plot_path_run): os.makedirs(self.plot_path_run)
        if not os.path.isdir(self.temp_path_run): os.makedirs(self.temp_path_run)

        # Inside the results directory of this run, create a file named 'scaling.dat' to contain the runtimes
        # from which the scaling behaviour can be inferred
        self.create_scaling_file()

        # Inide the results directory of this run, create a file which gives useful information about this run
        self.create_info_file()

        ### SKIRT arguments

        self.arguments = SkirtArguments()

        # Ski file options
        self.arguments.ski_pattern = self.config.ski_path
        self.arguments.recursive = False
        self.arguments.relative = False

        # Input path for the simulation
        self.arguments.input_path = self.input_path

        # The output path is adjusted seperately for each simulation

        # Other options
        self.arguments.emulate = False
        self.arguments.single = True

        # Options for logging
        self.arguments.logging.brief = False
        self.arguments.logging.verbose = True
        if self.config.extraction.memory:
            self.arguments.logging.memory = True
            self.arguments.logging.allocation = True
            self.arguments.logging.allocation_limit = 1e-5

        # Options for parallelization are adjusted seperately for each simulation

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing for the next scaling test run")

        # Set default values for attributes
        self.simulations = []

        # Clear the scaling analyser
        self.analyser.clear()

        # Other ...

    # -----------------------------------------------------------------

    def perform(self):

        """
        This function ...
        :return:
        """

        # Log the remote host name, the parallelization mode and the version of SKIRT used for this test
        self.log.info("Starting scaling test run " + self.scaling_run_name + ":")
        self.log.info("  - remote host: " + self.config.remote)
        if self.config.scheduler: self.log.info("  - cluster: " + self.config.cluster)
        self.log.info("  - parallelization mode: " + self.config.mode)
        self.log.info("Using " + self.remote.skirt_version)

        # Perform the simulations with increasing number of processors
        processors = self.min_processors
        while processors <= self.max_processors:

            # Perform this run
            self.launch(processors)

            # The next run will be performed with double the amount of processors
            processors *= 2

        if not self.scheduler:
            shell_script_path = os.path.join(self.temp_path_run, "simulations.sh")
            self.remote.start_queue(self.scaling_run_name, shell_script_path)

    # -----------------------------------------------------------------

    def retreive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Retreiving finished simulations...")

        # Get a list of the simulations that have been succesfully retreived
        self.simulations = self.remote.retreive()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Analysing retreived simulations...")

        # Run the scaling analyser for the current ski file
        self.analyser.run(self.config.ski_path)

    # -----------------------------------------------------------------

    def launch(self, processors):

        """
        This function performs one simulation of the scaling test with a particular number of processors
        :param processors:
        :return:
        """

        # Open the info file
        infofile = open(self.info_file_path, 'a')

        # Determine the number of processes and threads per process
        processes, threads = self.get_mapping(processors)

        # Determine the paths to the simulation's output, result and plot directories
        self.output_path_simulation = os.path.join(self.output_path_run, str(processors))
        self.result_path_simulation = os.path.join(self.result_path_run, str(processors))
        self.plot_path_simulation = os.path.join(self.plot_path_run, str(processors))

        # Create the output, result and plot directories for this run if necessary
        if not os.path.isdir(self.output_path_simulation): os.makedirs(self.output_path_simulation)
        if not os.path.isdir(self.result_path_simulation): os.makedirs(self.result_path_simulation)
        if not os.path.isdir(self.plot_path_simulation): os.makedirs(self.plot_path_simulation)

        # Write some information about this simulation to the info file
        infofile.write("Simulation performed on " + str(processors) + " processors\n")
        infofile.write(" - ski file: " + self.config.ski_path + "\n")
        infofile.write(" - output directory: " + self.output_path_simulation + "\n")
        infofile.write(" - number of processes: " + str(processes) + "\n")
        infofile.write(" - number of threads per processes: " + str(threads) + "\n")

        # Schedule or launch the simulation
        #if self.config.scheduler: self.schedule(processors, processes, threads, infofile)
        #else: self.execute(processors, processes, threads, infofile)

        self.run_simulation(processors, processes, threads, infofile)

        # Close the info file
        infofile.write("\n")
        infofile.close()

    # -----------------------------------------------------------------

    def run_simulation(self, processors, processes, threads, infofile):

        """
        This function schedules a simulation on the cluster. This function takes the following arguments:
        :param processors: the total number of processors to be used for this run
        :param processes:
        :param threads:
        :param skifilepath:
        :param dataoutputpath:
        :param infofile:
        :return:
        """

        # Determine the number of nodes and processors per node
        nodes, ppn = self.get_requirements(processors)

        # In threads mode, show a warning message if the number of threads > the number of cores per node
        # (we can't use multiple nodes in threads mode)
        if self.config.mode == "threads" and threads > self.cores:

            # Show a warning and return immediately
            self.log.warning("The number of threads " + str(threads) + " exceeds the number of cores per node: skipping")
            return

        # Inform the user about the number of processors, processes, threads per process, nodes and processors per node
        self.log.info("Scheduling simulation with:")
        self.log.info(" - total number of processors = " + str(processors))
        self.log.info(" - number of parallel processes = " + str(processes))
        self.log.info(" - number of parallel threads per process = " + str(threads))
        self.log.info(" - number of nodes = " + str(nodes))
        self.log.info(" - number of requested processors per node = " + str(ppn))

        # Write the number of nodes and processors per node to the info file
        infofile.write(" - number of used nodes: " + str(nodes) + "\n")
        infofile.write(" - number of requested processors per node: " + str(ppn) + "\n")

        # The path of the log file for this simulation run
        #logfilepath = os.path.join(dataoutputpath, self._skifilename + "_log.txt")

        # Calculate the expected walltime for this number of processors if a scheduling system is used
        if self.scheduler: walltime = self.estimate_walltime(processors)
        else: walltime = None

        ###

        # Adjust the SKIRT command-line arguments
        self.arguments.parallel.processes = processes
        self.arguments.parallel.threads = threads

        # The local output path
        self.arguments.output_path = self.output_path_simulation

        # Create a unique name for this simulation, based on the scaling run name and the current number of processors
        simulation_name = self.scaling_run_name + "_" + str(processors)

        # Run the simulation
        if self.scheduler:
            # Create the job script. The name of the script indicates the mode in which we run this scaling test and
            # the current number of processors used. We enable the SKIRT verbose logging mode to be able to compare
            # the progress of the different parallel processes afterwards. Because for scaling tests, we don't want
            # processes to end up on different nodes or the SKIRT processes sensing interference from other programs,
            # we set the 'fullnode' flag to True, which makes sure we always request at least one full node, even when
            # the current number of processors is less than the number of cores per node
            jobscript_path = os.path.join(self.temp_path_run, "job_" + str(processors) + ".sh")
            simulation_file_path = self.remote.run(self.arguments, walltime, simulation_name, jobscript_path, mail=False, full_node=True)
        else:
            simulation_file_path = self.remote.add_to_queue(self.arguments)

        # Add additional information to the simulation file
        simulation_file = open(simulation_file_path, 'a')
        simulation_file.write("extract progress: " + str(self.config.extraction.progress) + "\n")
        simulation_file.write("extract timeline: " + str(self.config.extraction.timeline) + "\n")
        simulation_file.write("extract memory: " + str(self.config.extraction.memory) + "\n")
        simulation_file.write("plot progress: " + str(self.config.plotting.progress) + "\n")
        simulation_file.write("plot timeline: " + str(self.config.plotting.timeline) + "\n")
        simulation_file.write("plot memory: " + str(self.config.plotting.memory) + "\n")
        simulation_file.write("remove remote input: " + str(False) + "\n")
        simulation_file.write("remove remote output: " + str(False) + "\n")
        simulation_file.write("extraction directory: " + self.result_path_simulation + "\n")
        simulation_file.write("plotting directory: " + self.plot_path_simulation + "\n")
        simulation_file.write("part of scaling test " + self.scaling_run_name + "\n")

        # Close the file
        simulation_file.close()

        # Add information about the path to the directory where the extracted data will be placed
        if self.config.extraction.progress: infofile.write(" - progress information will be extracted to: " + self.result_path_simulation + "\n")
        if self.config.extraction.timeline: infofile.write(" - timeline information will be extracted to: " + self.result_path_simulation + "\n")
        if self.config.extraction.memory: infofile.write(" - memory information will be extract to: " + self.result_path_simulation)

        # Add information about the path to the directory where the plots will be placed
        if self.config.plotting.progress: infofile.write(" - progress data will be plotted to: " + self.plot_path_simulation + "\n")
        if self.config.plotting.timeline: infofile.write(" - timeline data will be plotted to: " + self.plot_path_simulation + "\n")
        if self.config.plotting.memory: infofile.write(" - memory data will be plotted to: " + self.plot_path_simulation + "\n")

    # -----------------------------------------------------------------

    def create_info_file(self):

        """
        This function creates a file containing general information about the current scaling test run
        :return:
        """

        # Create the file and set the path
        self.info_file_path = os.path.join(self.result_path_run, "info.txt")
        infofile = open(self.info_file_path, "w")

        # If hybrid mode is selected, add the number of threads per process to the name of the results directory
        hybridinfo = " with " + str(self.threads_per_process) + " threads per process" if self.config.mode == "hybrid" else ""

        # Write some useful information to the file
        infofile.write("Scaling benchmark test " + self.scaling_run_name + "\n")
        infofile.write("Remote host: " + self.config.remote + "\n")
        infofile.write("SKIRT version: " + self.remote.skirt_version + "\n")
        infofile.write("Parallelization mode: " + self.config.mode + hybridinfo + "\n")
        infofile.write("Maximum number of nodes: " + str(self.config.max_nodes) + " (" + str(self.max_processors) + " processors)\n")
        infofile.write("Minimum number of nodes: " + str(self.config.min_nodes) + " (" + str(self.min_processors) + " processors)\n")
        infofile.write("\n")

        # Close the info file (information on specific simulations will be appended)
        infofile.close()

    # -----------------------------------------------------------------

    def create_scaling_file(self):

        """
        This function creates the file to contain the information about the scaling behaviour
        :return:
        """

        # Create the file and set the path
        self.scaling_file_path = os.path.join(self.result_path_run, "scaling.dat")
        scalingfile = open(self.scaling_file_path, "w")

        # Write a header to this new file which contains some general info about its contents
        scalingfile.write("# Timing results for scaling benchmark test " + self.scaling_run_name + "\n")
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

    # -----------------------------------------------------------------

    def get_mapping(self, processors):

        """
        This function calculates the number of processes and the number of threads (per process) for
        a certain number of processors, depending on the mode in which this scaling test is run.
        In other words, this function determines the 'mapping' from a set of processors to an appropriate
        set of threads and processes. This function takes the number of processors as the sole argument.
        :param processors:
        :return:
        """

        # Set default values for the number of threads and processes
        threads = 1
        processes = 1

        # In mpi mode, each processor runs a different process
        if self.config.mode == "mpi": processes = processors

        # In threads mode, each processor runs a seperate thread within the same process
        if self.config.mode == "threads": threads = processors

        # In hybrid mode, the number of processes depends on how many threads are requested per process
        # and the current number of processors
        if self.config.mode == "hybrid":

            threads = self.threads_per_process
            processes = processors // self.threads_per_process

        # Return the number of processes and the number of threads
        return processes, threads

    # -----------------------------------------------------------------

    def get_requirements(self, processors):

        """
        This function calculates the required amount of nodes and processors per node, given a certain number of
        processors.
        :param processors:
        :return:
        """

        # Calculate the necessary amount of nodes
        nodes = processors // self.cores + (processors % self.cores > 0)

        # Determine the number of processors per node
        ppn = processors if nodes == 1 else self.cores

        # Return the number of nodes and processors per node
        return nodes, ppn

    # -----------------------------------------------------------------

    def estimate_walltime(self, processors, factor=1.5):

        """
        This function estimates the total runtime (walltime) for the current simulation, number of processors,
        system and parallelization mode. This function takes the following arguments:
        :param processors: the number of processors (or total threads) used for this run of the simulation
        :param factor: this optional argument determines how much the upper limit on the walltime should deviate from
            a previous run of the same simulation.
        :return:
        """

        # Try to get the runtimes for this number of processors
        runtimes = self.get_runtimes(processors)

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
            runtimes = self.get_runtimes(1, anysystem=True, anymode=True)

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
            logfilepath = os.path.join(self.base_path, self.config.ski_path + "_log.txt")

            # Check whether such a file exists
            if os.path.exists(logfilepath):

                # Initially, set the number of processes and threads from the log file to one
                logfileprocesses = 1
                logfilethreads = 1

                # Check with how many processes and threads this simulation was run, by reading each line of the
                # specified log file and searching for indications of multiple processes and/or multiple threads
                for line in open(logfilepath):

                    if 'Starting simulation ' + os.path.basename(self.config.ski_path).split(".")[0] + ' with' in line:

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

                self.log.error("The walltime could not be estimated for a run with " + str(processors) + " processors.")
                exit()

    # -----------------------------------------------------------------

    def get_runtimes(self, processors, anysystem=False, anymode=False):

        """
        This function extracts the timings for the current simulation from a scaling test results file that was
        created earlier for the current simulation. This function takes the following arguments:
        :param processors: the number of processors (or total threads) for which to look up the runtimes
        :param anysystem: this flag tells whether we must look for results of any system, or only from the system
                we are currently running on
        :param anymode: this flag tells whether we must look for results created with any mode, or only with the mode
                in which we are running the current scaling test
        :return:
        """

        # Search for a scaling results file corresponding with this system and the current mode
        for itemname in os.listdir(self.result_path):

            # Define the full path to this item
            itempath = os.path.join(self.result_path, itemname)

            # Check whether this item is a directory and it is not hidden
            if not os.path.isdir(itempath) or itemname.startswith("."): continue

            # Define the name of the scaling results file inside this directory
            filepath = os.path.join(itempath, "scaling.dat")

            # Split the directory name into its segments
            segments = itemname.split("_")

            # Get the system name in which the scaling test was run for this results file
            systemname = segments[0]

            # Get the mode in which the scaling test was run for this results file
            mode = segments[1]

            # Check whether this results file corresponds to a test on this system and the current mode
            if (anysystem or systemname == self.config.remote) and (anymode or mode.startswith(self.config.mode)):

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
