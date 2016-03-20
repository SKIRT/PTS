#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.scaling Contains the ScalingTest class, used for launching a SKIRT scaling test for a
#  particular ski file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..simulation.simulation import SkirtSimulation
from ..simulation.arguments import SkirtArguments
from ..launch.analyser import SimulationAnalyser
from .resources import ResourceEstimator
from ..basics.configurable import Configurable
from ..simulation.remote import SkirtRemote
from ..extract.timeline import TimeLineExtractor
from ..tools import time, filesystem
from ..tools.logging import log

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
        super(ScalingTest, self).__init__(config, "core")

        # -- Attributes --

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # Initialize a list to contain the retrieved simulations
        self.simulations = []

        # The number of cores on the remote system
        self.cores = None

        # The number of (hyper)threads to be used per core
        self.threads_per_core = None

        # The minimum and maximum number of processors to use for the scaling test
        self.min_processors = None
        self.max_processors = None

        # Information about the parallelization mode used for the scaling test
        self.mode_info = None
        self.mode_info_long = None

        # The simulation prefix
        self.prefix = None

        # The base path
        self.base_path = None

        # The name of the scaling run
        self.scaling_run_name = None
        self.long_scaling_run_name = None

        # The SKIRT arguments object
        self.arguments = None

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
        test.config.keep = arguments.keep

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

        # 3. Retrieve the output of finished simulations (of other scaling test runs)
        self.retrieve()

        # 3. Analyse the output of the retrieved simulations
        self.analyse()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ScalingTest, self).setup()

        # -- Remote execution environment --

        # Setup the remote SKIRT execution context
        self.remote.setup(self.config.remote, self.config.cluster)

        # Get the number of cores (per node) on the remote system
        # (cache it because it requires a external call each time)
        self.cores = self.remote.cores

        # Determine how many threads that we want to use per core, depending on the number of hyperthreads per core
        # on the remote system and whether hyperthreading is enabled in the remote host configuration file
        if self.remote.use_hyperthreading:

            # Inform the user
            log.info("Hybrid and pure-threading scaling tests will be performed with hyperthreading enabled")

            # Set the number of threads per core equal to the number of hyperthreads on the system
            self.threads_per_core = self.remote.threads_per_core

        # If hyperthreading is not enabled for the remote host, set the number of threads per core to one
        else: self.threads_per_core = 1

        # -- The minimum and maximum number of processors --

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        self.max_processors = int(self.config.max_nodes * self.cores)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one
        self.min_processors = int(self.config.min_nodes * self.cores)
        if self.min_processors == 0: self.min_processors = 1

        # In hybrid mode, the minimum number of processors also represents the number of threads per process
        self.threads_per_process = 1
        if self.config.mode == "hybrid": self.threads_per_process = self.min_processors * self.threads_per_core

        # Create a string that identifies the parallelization mode, where the number of threads per process has been
        # appended for hybrid mode
        self.mode_info = self.config.mode + "-" + str(self.threads_per_process) if self.config.mode == "hybrid" else self.config.mode
        self.mode_info_long = self.config.mode + " mode with " + str(self.threads_per_process) + " threads per process" if self.config.mode == "hybrid" else self.config.mode + " mode"

        # -- Directory structure --

        # Determine the simulation prefix
        self.prefix = filesystem.name(self.config.ski_path).split(".")[0]

        # Set the base path (the ski file directory)
        self.base_path = filesystem.directory_of(self.config.ski_path) if "/" in self.config.ski_path else filesystem.cwd()

        # Define a name identifying this scaling test run
        self.scaling_run_name = time.unique_name(self.mode_info, separator="__")
        self.long_scaling_run_name = "SKIRT__scaling__" + self.prefix + "__" + self.remote.system_name + "__" + self.scaling_run_name

        # Create the input, output, result, plot and temp directories
        self.create_directory_structure()

        # Inside the results directory of this run, create a file named 'scaling.dat' to contain the runtimes
        # from which the scaling behaviour can be inferred
        self.create_scaling_file()

        # Inide the results directory of this run, create a file which gives useful information about this run
        self.create_info_file()

        # -- SKIRT arguments --

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
        self.arguments.logging.memory = True
        #self.arguments.logging.allocation = True
        #self.arguments.logging.allocation_limit = 1e-6

        # Options for parallelization are adjusted seperately for each simulation

    # -----------------------------------------------------------------

    def create_directory_structure(self):

        """
        This function ...
        :return:
        """

        # -- Top - level directories --

        # Set the paths to the input, output, result, plot and temp directories
        self.input_path = filesystem.join(self.base_path, "in")
        self.output_path = filesystem.join(self.base_path, "out")
        self.result_path = filesystem.join(self.base_path, "res")
        self.plot_path = filesystem.join(self.base_path, "plot")
        self.temp_path = filesystem.join(self.base_path, "temp")

        # Check if an input directory exists
        if not filesystem.is_directory(self.input_path): self.input_path = None

        # Create the output, result, plot and temp directories if necessary
        filesystem.create_directories([self.output_path, self.result_path, self.plot_path, self.temp_path])

        # -- System - level directories --

        # Set the input, output, result, plot and temp paths for the system we are running this scaling test on
        self.output_path_system = filesystem.join(self.output_path, self.remote.system_name)
        self.result_path_system = filesystem.join(self.result_path, self.remote.system_name)
        self.plot_path_system = filesystem.join(self.plot_path, self.remote.system_name)
        self.temp_path_system = filesystem.join(self.temp_path, self.remote.system_name)

        # Create the output, result, plot and temp directories for the system if necessary
        filesystem.create_directories([self.output_path_system, self.result_path_system, self.plot_path_system, self.temp_path_system])

        # -- Scaling run - level directories --

        # Determine the paths to the directories that will contain the output, results, plots and temporary files of this particular scaling test run
        self.output_path_run = filesystem.join(self.output_path_system, self.scaling_run_name)
        self.result_path_run = filesystem.join(self.result_path_system, self.scaling_run_name)
        self.plot_path_run = filesystem.join(self.plot_path_system, self.scaling_run_name)
        self.temp_path_run = filesystem.join(self.temp_path_system, self.scaling_run_name)

        # Create the output, result, plot and temp directories for this run if necessary
        filesystem.create_directories([self.output_path_run, self.result_path_run, self.plot_path_run, self.temp_path_run])

    # -----------------------------------------------------------------

    def create_simulation_directories(self, processors):

        """
        This function ...
        :param processors:
        :return:
        """

        # Determine the paths to the simulation's output, result and plot directories
        self.output_path_simulation = filesystem.join(self.output_path_run, str(processors))
        self.result_path_simulation = filesystem.join(self.result_path_run, str(processors))
        self.plot_path_simulation = filesystem.join(self.plot_path_run, str(processors))

        # Create the output, result and plot directories for this simulation if necessary
        filesystem.create_directories([self.output_path_simulation, self.result_path_simulation, self.plot_path_simulation])

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing for the next scaling test run ...")

        # Set default values for attributes
        self.simulations = []

        # Clear the simulation analyser
        self.analyser.clear()

    # -----------------------------------------------------------------

    def perform(self):

        """
        This function ...
        :return:
        """

        # Log the remote host name, the parallelization mode and the version of SKIRT used for this test
        log.info("Starting scaling test run " + self.scaling_run_name + ":")
        log.info("  - remote host: " + self.remote.system_name)
        log.info("  - parallelization mode: " + self.mode_info_long)
        log.info("Using " + self.remote.skirt_version)

        # Perform the simulations with increasing number of processors
        processors = self.min_processors
        while processors <= self.max_processors:

            # Perform this run
            self.launch(processors)

            # The next run will be performed with double the amount of processors
            processors *= 2

        # If the remote host does not use a scheduling system, manually start the queued simulations
        if not self.remote.scheduler:

            # Determine a local path for the batch script for manual inspection
            shell_script_path = filesystem.join(self.temp_path_run, "simulations.sh")
            self.remote.start_queue(self.long_scaling_run_name, shell_script_path)

        # End with some log messages
        log.info("Finished scaling test run")
        log.info("The results are / will be written to " + self.scaling_file_path)

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving finished simulations ...")

        # Get a list of the simulations that have been succesfully retrieved
        self.simulations = self.remote.retrieve()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing retrieved simulations ...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation)

            # Clear the analyser
            self.analyser.clear()

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

        # Create the directories
        self.create_simulation_directories(processors)

        # Write some information about this simulation to the info file
        infofile.write("Simulation performed on " + str(processors) + " processor(s)\n")
        infofile.write(" - ski file: " + self.config.ski_path + "\n")
        infofile.write(" - output directory: " + self.output_path_simulation + "\n")
        infofile.write(" - number of processes: " + str(processes) + "\n")
        infofile.write(" - number of threads per processes: " + str(threads) + "\n")
        infofile.write(" - number of (hyper)threads per core: " + str(self.threads_per_core) + "\n")

        # Schedule or launch the simulation
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
        :param infofile:
        :return:
        """

        # Determine the number of nodes and processors per node
        nodes, ppn = self.get_requirements(processors)

        # In threads mode, show a warning message if the number of threads > the number of cores per node
        # (we can't use multiple nodes in threads mode)
        if self.config.mode == "threads" and threads > self.cores:

            # Show a warning and return immediately
            log.warning("The number of threads " + str(threads) + " exceeds the number of cores per node: skipping")
            return

        # Inform the user about the number of processors, processes, threads per process, nodes and processors per node
        log.info("Scheduling simulation with:")
        log.info(" - total number of processors = " + str(processors))
        log.info(" - number of parallel processes = " + str(processes))
        log.info(" - number of parallel threads per process = " + str(threads))
        log.info(" - number of nodes = " + str(nodes))
        log.info(" - number of requested processors per node = " + str(ppn))

        # Write the number of nodes and processors per node to the info file
        infofile.write(" - number of used nodes: " + str(nodes) + "\n")
        infofile.write(" - number of requested processors per node: " + str(ppn) + "\n")

        # Calculate the expected walltime for this number of processors if a scheduling system is used
        if self.remote.scheduler:
            walltime = self.estimate_walltime(processes, threads)
            log.info(" - expected walltime: " + str(walltime) + " seconds")
        else: walltime = None

        ###

        # Adjust the SKIRT command-line arguments
        self.arguments.parallel.processes = processes
        self.arguments.parallel.threads = threads

        # The local output path
        self.arguments.output_path = self.output_path_simulation

        # Create a unique name for this simulation, based on the scaling run name and the current number of processors
        simulation_name = self.long_scaling_run_name + "_" + str(processors)

        # Run the simulation
        scheduling_options = None
        if self.remote.scheduler:

            # Create the job script. The name of the script indicates the mode in which we run this scaling test and
            # the current number of processors used. We enable the SKIRT verbose logging mode to be able to compare
            # the progress of the different parallel processes afterwards. Because for scaling tests, we don't want
            # processes to end up on different nodes or the SKIRT processes sensing interference from other programs,
            # we set the 'fullnode' flag to True, which makes sure we always request at least one full node, even when
            # the current number of processors is less than the number of cores per node
            jobscript_path = filesystem.join(self.temp_path_run, "job_" + str(processors) + ".sh")

            # Adjust the scheduling options
            scheduling_options = {}
            scheduling_options["nodes"] = nodes
            scheduling_options["ppn"] = ppn
            scheduling_options["walltime"] = walltime
            scheduling_options["jobscript_path"] = jobscript_path
            scheduling_options["mail"] = False
            scheduling_options["full_node"] = True

        # Add the simulation to the remote queue
        simulation = self.remote.add_to_queue(self.arguments, simulation_name, scheduling_options)

        # Add additional information to the simulation object
        simulation.extract_progress = True
        simulation.extract_timeline = True
        simulation.extract_memory = True
        simulation.plot_progress = True
        simulation.plot_timeline = True
        simulation.plot_memory = True
        simulation.remove_remote_input = not self.config.keep
        simulation.remove_remote_output = not self.config.keep
        simulation.remove_remote_simulation_directory = not self.config.keep
        simulation.retrieve_types = ["log"]
        simulation.extraction_path = self.result_path_simulation
        simulation.plotting_path = self.plot_path_simulation
        simulation.scaling_run_name = self.long_scaling_run_name
        simulation.scaling_data_file = self.scaling_file_path
        simulation.scaling_plot_path = self.plot_path_system
        if not self.remote.scheduler: simulation.screen_name = self.long_scaling_run_name

        # Save the simulation file
        simulation.save()

        # Add information about the path to the directory where the extracted data will be placed
        infofile.write(" - progress information will be extracted to: " + self.result_path_simulation + "\n")
        infofile.write(" - timeline information will be extracted to: " + self.result_path_simulation + "\n")
        infofile.write(" - memory information will be extract to: " + self.result_path_simulation + "\n")

        # Add information about the path to the directory where the plots will be placed
        infofile.write(" - progress data will be plotted to: " + self.plot_path_simulation + "\n")
        infofile.write(" - timeline data will be plotted to: " + self.plot_path_simulation + "\n")
        infofile.write(" - memory data will be plotted to: " + self.plot_path_simulation + "\n")

    # -----------------------------------------------------------------

    def create_info_file(self):

        """
        This function creates a file containing general information about the current scaling test run
        :return:
        """

        # Create the file and set the path
        self.info_file_path = filesystem.join(self.result_path_run, "info.txt")
        infofile = open(self.info_file_path, "w")

        # Write some useful information to the file
        infofile.write("Scaling benchmark test " + self.scaling_run_name + "\n")
        infofile.write("Remote host: " + self.remote.system_name + "\n")
        infofile.write("SKIRT version: " + self.remote.skirt_version + "\n")
        infofile.write("Parallelization mode: " + self.mode_info_long + "\n")
        infofile.write("Maximum number of nodes: " + str(self.config.max_nodes) + " (" + str(self.max_processors) + " processors)\n")
        infofile.write("Minimum number of nodes: " + str(self.config.min_nodes) + " (" + str(self.min_processors) + " processor(s))\n")
        infofile.write("\n")

        # Close the info file (information on specific simulations will be appended)
        infofile.close()

    # -----------------------------------------------------------------

    def create_scaling_file(self):

        """
        This function creates the file to contain the information about the scaling behaviour
        :return:
        """

        # Set the path to the scaling file for the current system (remote host - cluster)
        self.scaling_file_path = filesystem.join(self.result_path_system, "scaling.dat")

        # If the file has already been created, skip the rest of the function
        if filesystem.is_file(self.scaling_file_path): return

        names = []
        names.append("Parallelization mode")        # Parallelization mode
        names.append("Processes")                   # Number of processes p
        names.append("Threads")                     # Number of threads (per process) t
        names.append("Setup time")                  # Time spent in simulation setup (s)
        names.append("Stellar emission time")       # Time spent shooting stellar photon packages (s)
        names.append("Spectra calculation time")    # Time spent in calculation of dust emission spectra (s)
        names.append("Dust emission time")          # Time spent shooting dust emission photon packages (s)
        names.append("Writing time")                # Time spent writing to disk (s)
        names.append("Waiting time")                # Time spent waiting for other processes (s)
        names.append("Communication time")          # Time spent in inter-process communication (s)
        names.append("Total time")                  # Total simulation time (s)
        names.append("Peak memory usage")           # Peak memory usage (GB)

        # Create an empty table (it just has the header)
        table = Table(names=names)

        table["Setup time"].unit = "s"
        table["Stellar emission time"].unit = "s"
        table["Spectra calculation time"].unit = "s"
        table["Dust emission time"].unit = "s"
        table["Writing time"].unit = "s"
        table["Waiting time"].unit = "s"
        table["Communication time"].unit = "s"
        table["Total time"].unit = "s"
        table["Peak memory usage"].unit = "GB"

        # Write the table to file
        table.write(self.scaling_file_path, format="ascii.commented_header")

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
        if self.config.mode == "threads": threads = processors * self.threads_per_core

        # In hybrid mode, the number of processes depends on how many threads are requested per process
        # and the current number of processors
        if self.config.mode == "hybrid":

            threads = self.threads_per_process
            cores_per_process = self.self.threads_per_process / self.threads_per_core
            processes = processors // cores_per_process

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

    def estimate_walltime(self, processes, threads, factor=1.2):

        """
        This function estimates the total runtime (walltime) for the current simulation, number of processors,
        system and parallelization mode. This function takes the following arguments:
        :param processes: the number of processes used for this run of the simulation
        :param threads: the number of threads used for this run of the simulation
        :param factor: this optional argument determines how much the upper limit on the walltime should deviate from
            a previous run of the same simulation.
        :return:
        """

        processors = processes * threads

        # Create a dictionary to contain the paths to timeline data files found for the ski file,
        # indexed on (system_name, mode, processors)
        timeline_paths = defaultdict(list)

        # Recursively search for files contained in the result directory
        for file_path, file_name in filesystem.files_in_path(self.result_path, contains="timeline", extension="dat", returns=["path", "name"], recursive=True):

            # Get the path to the directory where this timeline.dat file is in
            dir_path = filesystem.directory_of(file_path)
            dir_of_dir_path = filesystem.directory_of(dir_path)
            dir_of_dir_of_dir_path = filesystem.directory_of(dir_path)

            # Get properties of scaling test run
            processors = int(filesystem.name(dir_path))
            scaling_run_name = filesystem.name(dir_of_dir_path)
            system_name = filesystem.name(dir_of_dir_of_dir_path)
            mode = scaling_run_name.split("_")[0]

            # Add the timeline file path to the dictionary
            timeline_paths[(system_name, mode, processors)].append(file_path)

        # 1. Try to find an extracted timeline for the system and parallelization mode of this run and the current
        #    number of processors
        key = (self.remote.system_name, self.mode_info, processors)
        if key in timeline_paths:

            # Create a TimeLineExtractor instance (from the first timeline file in the list)
            extractor = TimeLineExtractor.open_table(timeline_paths[key][0])

            # Return the total runtime (multiplied by the specified factor)
            return extractor.total * factor

        # 2. Try to find an extracted timeline for the system of this run and the current number of processors, but
        #    a different parallelization mode
        for key in timeline_paths:

            # Check whether the system name and number of processors in the key correspond to those of the current scaling test run
            if key[0] == self.remote.system_name and key[2] == processors:

                # Create a TimeLineExtractor instance (from the first timeline file in the list)
                extractor = TimeLineExtractor.open_table(timeline_paths[key][0])

                # Return the total runtime (multiplied by the specified factor)
                return extractor.total * factor

        # 3. Try to find an extracted timeline for the current number of processors, but for a different system
        for key in timeline_paths:

            # Check whether the number of processors in the key corresponds to the current one
            if key[2] == processors:

                # Create a TimeLineExtractor instance (from the first timeline file in the list)
                extractor = TimeLineExtractor.open_table(timeline_paths[key][0])

                # Return the total runtime (multiplied by the specified factor)
                return extractor.total * factor

        # 4. Try to find a log file placed next to the ski file used for the scaling test
        log_file_path = filesystem.join(self.base_path, self.prefix + "_log.txt")
        if filesystem.is_file(log_file_path):

            # Create a SkirtSimulation object
            simulation = SkirtSimulation(self.prefix, self.input_path, self.base_path)

            # Create a new TimeLineExtractor instance
            extractor = TimeLineExtractor()
            extractor.run(simulation)

            # Determine the number of threads and processes for the log file
            log_processes = simulation.processes()
            log_threads = simulation.threads()

            # Determine the number of used processors (assuming each single threads was run on a seperate processor)
            log_processors = log_processes * log_threads

            # Return the estimated total runtime for the current number of processors (assuming the overhead increases linearly with the number of processors)
            return (extractor.serial + extractor.parallel * log_processors / processors + extractor.overhead / log_processors * processors) * factor

        # 5. Try to estimate the runtime by using the ResourceEstimator class
        else:

            # Create and run a ResourceEstimator instance
            estimator = ResourceEstimator()
            estimator.run(self.arguments.ski_pattern, processes, threads)

            # Return the estimated walltime
            return estimator.walltime * factor

# -----------------------------------------------------------------
