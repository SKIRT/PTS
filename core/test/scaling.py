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

# Import the relevant PTS classes and modules
from ..basics.map import Map
from ..simulation.definition import SingleSimulationDefinition
from ..simulation.parallelization import Parallelization
from ..basics.configurable import Configurable
from ..tools import time
from ..tools import filesystem as fs
from ..basics.log import log
from ..launch.options import SchedulingOptions
from ..launch.batch import BatchLauncher
from ..launch.timing import TimingTable
from ..launch.memory import MemoryTable
from ..tools import tables
from ..advanced.runtimeestimator import RuntimeEstimator

# -----------------------------------------------------------------

# To distribute the wavelengths over the processes based on some kind of 'load' curve as a function of wavelength

# from http://stackoverflow.com/questions/4632322/finding-all-possible-combinations-of-numbers-to-reach-a-given-sum

def subset_sum(numbers, target, partial=[]):

    s = sum(partial)

    # check if the partial sum is equals to target
    if s == target:
        print("sum(%s)=%s" % (partial, target))
    if s >= target:
        return  # if we reach the number why bother to continue

    for i in range(len(numbers)):
        n = numbers[i]
        remaining = numbers[i+1:]
        subset_sum(remaining, target, partial + [n])


#if __name__ == "__main__":
#    subset_sum([3,9,8,4,5,7,10],15)

    #Outputs:
    #sum([3, 8, 4])=15
    #sum([3, 5, 7])=15
    #sum([8, 7])=15
    #sum([5, 10])=15

# -----------------------------------------------------------------

class ScalingTest(Configurable):

    """
    An instance of the ScalingTest class represents a SKIRT scaling benchmark test for a particular ski file.
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(ScalingTest, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The system name
        self.system_name = None

        # The SKIRT version being used on the remote system
        self.skirt_version = None

        # Whether or not the remote host uses a scheduling system
        self.scheduler = None

        # The simulation prefix
        self.prefix = None

        # The path to the base directory
        self.base_path = None

        # The main directories
        self.main_paths = Map()

        # System-level directories
        self.system_paths = Map()

        # Run-level directories
        self.run_paths = Map()

        # The number of cores (per node) on the remote system
        self.cores_per_node = None

        # The number of (hyper)threads to be used per core
        self.threads_per_core = None

        # The number of threads per process, only relevant for hybrid scaling test runs
        self.threads_per_process = None

        # The minimum and maximum number of processing cores to use for the scaling test
        self.min_cores = None
        self.max_cores = None

        # Information about the parallelization mode used for the scaling test
        self.mode_info = None
        self.mode_info_long = None

        # The table with information about the different simulations
        self.info_table = None

        # The path to the info table
        self.info_table_path = None

        # The runtime estimator
        self.estimator = None

        # The paths of the timing and memory table
        self.timing_table_path = None
        self.memory_table_path = None

        # The name of the scaling run
        self.scaling_run_name = None
        self.long_scaling_run_name = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        When this function is called, the scaling test is started.
        """

        # 2. Load (and create) the info tables
        self.create_info_table()

        # 3. Set the runtime estimator
        if self.scheduler: self.set_estimator()

        # 4. Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ScalingTest, self).setup(**kwargs)

        # Basic setup
        self.setup_basic()

        # Setup the batch launcher
        self.setup_launcher()

        # Setup properties
        self.setup_properties()

        # Setup the properties of this scaling run
        self.setup_run()

    # -----------------------------------------------------------------

    def setup_basic(self):

        """
        This function ...
        :return:
        """

        # Determine the simulation prefix
        self.prefix = fs.name(self.config.ski_path).split(".")[0]

        # Set the base path (the ski file directory)
        self.base_path = fs.directory_of(self.config.ski_path) if "/" in self.config.ski_path else fs.cwd()

        # Set the path to the timing table and initialize it if necessary
        self.timing_table_path = fs.join(self.base_path, "timing.dat")
        if not fs.is_file(self.timing_table_path):
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        # Set the path to the memory table and initialize it if necessary
        self.memory_table_path = fs.join(self.base_path, "memory.dat")
        if not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

    # -----------------------------------------------------------------

    def setup_launcher(self):

        """
        This function ...
        :return:
        """

        # Set options for the batch launcher: basic options
        self.launcher.config.shared_input = True  # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = [self.config.remote]  # the remote hosts on which to run the simulations
        self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file
        self.launcher.config.retrieve_types = ["log"] # only log files should be retrieved from the simulation output
        self.launcher.config.keep = self.config.keep

        # Options for the batch launcher: extraction analysis options
        self.launcher.config.analysis.extraction.path = self.run_paths.result
        self.launcher.config.analysis.extraction.timeline = True  # extract the simulation timeline
        self.launcher.config.analysis.extraction.progress = True  # extract the simulation phase progress information
        self.launcher.config.analysis.extraction.memory = True  # extract memory information

        # Options for the batch launcher: plotting analysis options
        self.launcher.config.analysis.plotting.path = self.run_paths.plot  # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.plotting.progress = True
        self.launcher.config.analysis.plotting.timeline = True
        self.launcher.config.analysis.plotting.memory = True
        # self.launcher.config.analysis.plotting.format = "png"  # plot in PNG format so that an animation can be made from the fit SEDs

        # Options for the batch launcher: miscellaneous analysis options
        self.launcher.config.analysis.misc.path = self.run_paths.result  # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output

        # Options for the batch launcher: logging options
        self.launcher.config.logging.brief = False
        self.launcher.config.logging.verbose = True
        self.launcher.config.logging.memory = True
        self.launcher.config.logging.allocation = True
        self.launcher.config.logging.allocation_limit = 1e-6

        # Set the cluster name
        self.launcher.set_cluster_for_host(self.config.remote, self.config.cluster)

        # Setup the remotes (so that we can get properties such as the number of cores)
        self.launcher.setup_remotes()

    # -----------------------------------------------------------------

    def setup_properties(self):

        """
        This function ...
        :return:
        """

        # Get the remote instance from the batch launcher
        remote = self.launcher.single_remote

        # Get the number of cores (per node) on the remote system
        self.cores_per_node = remote.cores_per_node

        # Check whether hyperthreading can be enabled for the remote
        hyperthreading = remote.use_hyperthreading
        if hyperthreading:

            # Inform the user
            log.info("Hybrid and pure-threading scaling tests will be performed with hyperthreading enabled")

            # Set the number of threads per core equal to the number of hyperthreads on the system
            self.threads_per_core = remote.threads_per_core

        # If hyperthreading is not enabled for the remote host, set the number of threads per core to one
        else: self.threads_per_core = 1

        # Get the name of the remote system
        self.system_name = remote.system_name

        # Get the SKIRT version
        self.skirt_version = remote.skirt_version

        # Check whether the remote host uses a scheduling system
        self.scheduler = remote.scheduler

    # -----------------------------------------------------------------

    def setup_run(self):

        """
        This function ...
        :return:
        """

        # Calculate the maximum number of processors to use for the scaling test (maxnodes can be a decimal number)
        self.max_cores = int(self.config.max_nodes * self.cores_per_node)

        # Calculate the minimum number of processors to use for the scaling test and set this as the value
        # to start the loop below with. The default setting is minnodes and minprocessors both equal to zero. If
        # this is the case, the starting value for the loop is set to one
        self.min_cores = int(self.config.min_nodes * self.cores_per_node)
        if self.min_cores == 0: self.min_cores = 1

        # In hybrid mode, the minimum number of processors also represents the number of threads per process
        self.threads_per_process = 1
        if self.config.mode == "hybrid": self.threads_per_process = self.min_cores * self.threads_per_core

        # Create a string that identifies the parallelization mode, where the number of threads per process has been
        # appended for hybrid mode
        self.mode_info = self.config.mode + "-" + str(self.threads_per_process) if self.config.mode == "hybrid" else self.config.mode
        self.mode_info_long = self.config.mode + " mode with " + str(self.threads_per_process) + " threads per process" if self.config.mode == "hybrid" else self.config.mode + " mode"

        # -- Directory structure --

        # Define a name identifying this scaling test run
        self.scaling_run_name = time.unique_name(self.mode_info, separator="__")
        self.long_scaling_run_name = "SKIRT__scaling__" + self.prefix + "__" + self.system_name + "__" + self.scaling_run_name

        # Create the directory structure
        self.create_directories()

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Create the main directories
        self.create_main_directories()

        # Create the system-level directories
        self.create_system_directories()

        # Create the scaling run-level directories
        self.create_run_directories()

    # -----------------------------------------------------------------

    def create_main_directories(self):

        """
        This function ...
        :return:
        """

        # Set the paths to the input, output, result, plot and temp directories
        self.main_paths.input = fs.join(self.base_path, "in")
        self.main_paths.output = fs.join(self.base_path, "out")
        self.main_paths.result = fs.join(self.base_path, "res")
        self.main_paths.plot = fs.join(self.base_path, "plot")
        self.main_paths.temp = fs.join(self.base_path, "temp")

        # Check if an input directory exists
        if not fs.is_directory(self.main_paths.input): self.main_paths.input = None

        # Create the output, result, plot and temp directories if necessary
        fs.create_directories(self.main_paths.output, self.main_paths.result, self.main_paths.plot, self.main_paths.temp)

    # -----------------------------------------------------------------

    def create_system_directories(self):

        """
        This function ...
        :return:
        """

        # Set the input, output, result, plot and temp paths for the system we are running this scaling test on
        self.system_paths.output = fs.join(self.main_paths.output, self.system_name)
        self.system_paths.result = fs.join(self.main_paths.result, self.system_name)
        self.system_paths.plot = fs.join(self.main_paths.plot, self.system_name)
        self.system_paths.temp = fs.join(self.main_paths.temp, self.system_name)

        # Create the output, result, plot and temp directories for the system if necessary
        fs.create_directories(self.system_paths.output, self.system_paths.result, self.system_paths.plot, self.system_paths.temp)

    # -----------------------------------------------------------------

    def create_run_directories(self):

        """
        This function ...
        :return:
        """

        # Determine the paths to the directories that will contain the output, results, plots and temporary files of this particular scaling test run
        self.run_paths.output = fs.join(self.system_paths.output, self.scaling_run_name)
        self.run_paths.result = fs.join(self.system_paths.result, self.scaling_run_name)
        self.run_paths.plot = fs.join(self.system_paths.plot, self.scaling_run_name)
        self.run_paths.temp = fs.join(self.system_paths.temp, self.scaling_run_name)

        # Create the output, result, plot and temp directories for this run if necessary
        fs.create_directories(self.run_paths.output, self.run_paths.result, self.run_paths.plot, self.run_paths.temp)

    # -----------------------------------------------------------------

    def create_info_table(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the runs info table
        info_table_runs_path = fs.join(self.system_paths.result, "info.dat")

        # If the table does not yet exist, create a new one
        if not fs.is_file(info_table_runs_path):

            # Initialize the columns
            names = ["Scaling run name", "Remote", "SKIRT version", "Parallelization mode", "Max nodes", "Min nodes", "Max cores", "Min cores"]
            data = [[], [], [], [], [], [], [], []]
            #dtypes = ["S24", "float64", "float64", "float64"]
            dtypes = None

            # Create the table and save it
            info_table_runs = tables.new(data, names, dtypes=dtypes)
            tables.write(info_table_runs, info_table_runs_path, format="ascii.ecsv")

        # Add an entry to the runs info table
        with open(info_table_runs_path, 'a') as info_table_runs:

            info_table_runs.write(self.scaling_run_name + " " + self.system_name + " " + self.skirt_version + " "
                                  + self.mode_info_long + " " + str(self.config.max_nodes) + " "
                                  + str(self.config.min_nodes)) + " " + str(self.max_cores) + " " + str(self.min_cores)

        # Determine the path to the simulations info table
        self.info_table_path = fs.join(self.run_paths.result, "info.dat")

        # Initialize the columns
        names = ["Simulation name", "Cores", "Processes", "Threads per process", "Nodes", "Cores per node"]
        data = [[], [], [], [], [], []]
        #dtypes = []

        # Create the table
        self.info_table = tables.new(data, names)

    # -----------------------------------------------------------------

    def set_estimator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the runtime estimator ...")

        # Load the timing table
        self.estimator = RuntimeEstimator.from_file(self.timing_table_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Log the remote host name, the parallelization mode and the version of SKIRT used for this test
        log.info("Starting scaling test run " + self.scaling_run_name + ":")
        log.info("  - remote host: " + self.system_name)
        log.info("  - parallelization mode: " + self.mode_info_long)
        log.info("Using " + self.skirt_version)

        # Perform the simulations with increasing number of cores
        cores = self.min_cores
        while cores <= self.max_cores:

            # Perform this run
            self.launch_simulation(cores)

            # The next run will be performed with double the amount of cores
            cores *= 2

        # Set the local path for the batch script for manual inspection
        self.launcher.set_script_path(self.config.remote, self.run_paths.temp)

        # Enable screen output to be written out (for debugging)
        if not self.scheduler: self.launcher.enable_screen_output(self.config.remote)

        # Run the launcher, schedules or initiates the simulations
        self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in self.launcher.launched_simulations:

            # Set scaling test characteristics
            simulation.analysis.scaling_path = self.base_path
            simulation.analysis.scaling_run_name = self.long_scaling_run_name

            # Save the simulation object
            simulation.save()

        # End with some log messages
        log.success("Finished scaling test run")

    # -----------------------------------------------------------------

    def launch_simulation(self, cores):

        """
        This function performs one simulation of the scaling test with a particular number of processors
        :param cores:
        :return:
        """

        # Inform the user
        log.info("Adding a simulation to the queue with " + str(cores) + " cores ...")

        # Determine the parallelization scheme for the simulation
        parallelization = Parallelization.from_mode(self.config.mode, cores, self.threads_per_core, self.threads_per_process)

        # Create a unique name for this simulation, based on the scaling run name and the current number of processors
        simulation_name = self.long_scaling_run_name + "_" + str(parallelization.cores)

        # Create the simulation output directory
        # (paths for extraction, plotting and miscellaneous output are created by the batch launcher)
        output_path_simulation = fs.join(self.run_paths.output, simulation_name)
        fs.create_directory(output_path_simulation)

        # Create the SKIRT simulation definition
        definition = SingleSimulationDefinition(self.config.ski_path, output_path_simulation, self.main_paths.input, name=simulation_name)

        # Determine the required number of nodes (and in case of a single node, the number of cores on that node)
        nodes, ppn = parallelization.get_requirements(self.cores_per_node) # give the number of cores per node

        # In threads mode, show a warning message if the number of threads > the number of cores per node * hyperthreads per core
        # (we can't use multiple nodes in threads mode)
        if self.config.mode == "threads" and parallelization.threads > self.cores_per_node * self.threads_per_core:

            # Show a warning and return immediately
            log.warning("The number of threads " + str(parallelization.threads) + " exceeds the number of logical cores on a node: skipping")
            return

        # Debugging
        # Inform the user about the number of processors, processes, threads per process, nodes and processors per node
        log.debug(" - total number of processing cores = " + str(parallelization.cores))
        log.debug(" - number of parallel processes = " + str(parallelization.processes))
        log.debug(" - number of parallel threads per process = " + str(parallelization.threads))
        log.debug(" - number of nodes = " + str(nodes))
        log.debug(" - number of requested processors per node = " + str(ppn))

        # Add an entry to the simulations info table
        self.info_table.add_row([definition.name, parallelization.cores, parallelization.processes, parallelization.threads, nodes, ppn])

        # Add the simulation to the queue
        self.launcher.add_to_queue(definition, definition.name, parallelization)

        # Set the scheduling options if necessary
        if self.scheduler:

            # Calculate the expected walltime for this number of processors if a scheduling system is used
            # Estimate the walltime
            runtime = self.estimate_runtime(parallelization)
            log.info(" - expected runtime: " + str(runtime) + " seconds")

            # Create the job script. The name of the script indicates the mode in which we run this scaling test and
            # the current number of processors used. We enable the SKIRT verbose logging mode to be able to compare
            # the progress of the different parallel processes afterwards. Because for scaling tests, we don't want
            # processes to end up on different nodes or the SKIRT processes sensing interference from other programs,
            # we set the 'fullnode' flag to True, which makes sure we always request at least one full node, even when
            # the current number of processors is less than the number of cores per node
            jobscript_path = fs.join(self.run_paths.temp, "job_" + str(parallelization.cores) + ".sh")

            # Create a SchedulingOptions instance
            scheduling_options = SchedulingOptions()

            # Adjust the scheduling options
            scheduling_options.nodes = nodes
            scheduling_options.ppn = ppn
            scheduling_options.walltime = runtime
            scheduling_options.local_jobscript_path = jobscript_path
            scheduling_options.mail = False
            scheduling_options.full_node = True

            # Set scheduling options (for the different remote hosts with a scheduling system)
            self.launcher.set_scheduling_options(self.config.remote, definition.name, scheduling_options)

    # -----------------------------------------------------------------

    def estimate_runtime(self, parallelization):

        """
        This function ...
        :param parallelization:
        :return:
        """

        # Inform the user
        log.info("Estimating the walltime for the simulation ...")

        # Get, for this scaling run,
        #  - the number of photon packages
        #  - the number of wavelengths
        #  - the number of dust cells (?)

        # Estimate the runtime for the configured remote host and the simulation properties
        runtime = self.estimator.runtime_for(ski, parallelization, self.launcher.single_remote.host_id, self.launcher.single_remote.cluster_name)

        # Return the runtime
        return runtime

# -----------------------------------------------------------------
