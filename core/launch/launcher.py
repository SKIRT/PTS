#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.launcher Contains the SKIRTLauncher class, which can be used to launch SKIRT simulations
#  locally or remotely.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from ..simulation.execute import SkirtExec
from ..basics.configurable import Configurable
from ..tools import monitoring, introspection
from ..simulation.definition import SingleSimulationDefinition
from .options import LoggingOptions
from .analyser import SimulationAnalyser
from ..simulation.remote import SkirtRemote
from ..tools.logging import log
from .options import SchedulingOptions
from ..advanced.parallelizationtool import ParallelizationTool
from ..advanced.memoryestimator import MemoryEstimator
from ..simulation.parallelization import Parallelization
from .options import AnalysisOptions
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class SKIRTLauncher(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SKIRTLauncher, self).__init__(config, interactive)

        # -- Attributes --

        # Create the local SKIRT execution context
        self.skirt = SkirtExec()

        # Create the SKIRT remote execution context
        self.remote = None

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # The simulation definition
        self.definition = None

        # The logging options
        self.logging_options = None

        # The analysis options
        self.analysis_options = None

        # The parallelization scheme
        self.parallelization = None

        # The simulation object
        self.simulation = None

        # Initialize a list to contain the retrieved finished simulations
        self.simulations = []

        # Estimates of the memory requirement
        self.memory = None

    # -----------------------------------------------------------------

    @property
    def has_parallelization(self):

        """
        This function ...
        :return:
        """

        # Check whether the number of processes and the number of threads are both defined
        #return self.config.arguments.parallel.processes is not None and self.config.arguments.parallel.threads is not None
        return False

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the simulation definition
        self.create_definition()

        # 2. Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()
        else: self.check_parallelization()

        # 3. Launch the simulation
        self.launch()

        # 4. Retrieve the simulations that are finished
        if self.config.remote: self.retrieve()
        else: self.simulations.append(self.simulation) # add the locally run simulation to the list of simulations to be analysed

        # 5. Analyse the output of the retrieved simulations
        self.analyse()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SKIRTLauncher, self).setup(**kwargs)

        # Setup the remote execution context
        if self.config.remote is not None:
            self.remote = SkirtRemote()
            self.remote.setup(self.config.remote, self.config.cluster)

        # Create output directory
        if self.config.create_output and not fs.is_directory(self.config.output): fs.create_directory(self.config.output)

        # Create the logging options
        self.logging_options = LoggingOptions()
        self.logging_options.set_options(self.config.logging)

        # Create the analysis options
        if "analysis_options" in kwargs: self.set_analysis_options(kwargs.pop("analysis_options"))
        else: self.create_analysis_options()

        # Get the memory information passed to this instance
        self.memory = kwargs.pop("memory", None)

    # -----------------------------------------------------------------

    def create_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the analysis options ...")

        # Create the analysis options object
        self.analysis_options = AnalysisOptions()
        self.analysis_options.set_options(self.config.analysis)

        # Check the options
        self.analysis_options.check(self.logging_options, self.config.output)

    # -----------------------------------------------------------------

    def set_analysis_options(self, options):

        """
        This function ...
        :param options:
        :return:
        """

        # Inform the user
        log.info("Setting the analysis options ...")

        # Set
        self.analysis_options = options

        # Check the options
        self.analysis_options.check(self.logging_options, self.config.output)

    # -----------------------------------------------------------------

    def create_definition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the simulation definition ...")

        # Create the simulation definition
        self.definition = SingleSimulationDefinition(self.config.ski, self.config.output, self.config.input)

        #print(self.definition)

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        # Set parallelization
        if self.config.remote: self.set_parallelization_remote()
        else: self.set_parallelization_local()

    # -----------------------------------------------------------------

    def set_parallelization_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the optimal parallelization scheme ...")

        # Check whether MPI is available on this system
        if introspection.has_mpi():

            # If memory requirement is not set
            if self.memory is None:

                # The memory estimator
                estimator = MemoryEstimator()

                # Configure the memory estimator
                estimator.config.ski = self.definition.ski_path
                estimator.config.input = self.config.input
                estimator.config.show = False

                # Estimate the memory
                estimator.run()

                # Get the memory requirement
                self.memory = estimator.memory

            # Get the serial and parallel parts of the simulation's memory
            serial_memory = self.memory.serial
            parallel_memory = self.memory.parallel

            # Calculate the total memory of one process without data parallelization
            total_memory = serial_memory + parallel_memory

            # Calculate the maximum number of processes based on the memory requirements
            free_memory = monitoring.free_memory()
            processes = int(free_memory / total_memory)

            # If there is too little free memory for the simulation, the number of processes will be smaller than one
            if processes < 1:

                # Exit with an error
                log.error("Not enough memory available to run this simulation locally: free memory = " + str(free_memory) + ", required memory = " + str(total_memory))
                exit()

        # No MPI available
        else: processes = 1

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        threads = int(monitoring.free_cpus() / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:
            processes = max(int(monitoring.free_cpus()), 1)
            threads = 1

        # Set the parallelization options
        #self.config.arguments.parallel.processes = processes
        #self.config.arguments.parallel.threads = threads

        cores = processes * threads
        threads_per_core = 2

        # Debugging
        log.debug("The number of cores is " + str(cores))
        log.debug("The number of thread per core is " + str(threads_per_core))
        log.debug("The number of processes is " + str(processes))

        # Set the parallelization scheme
        self.parallelization = Parallelization(cores, threads_per_core, processes, data_parallel=False)

        # Debugging
        log.debug("The parallelization scheme is " + str(self.parallelization))

    # -----------------------------------------------------------------

    def set_parallelization_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for remote execution ...")

        # If the remote uses a scheduling system
        if self.remote.scheduler:

            # Set host properties
            nnodes = self.config.nnodes
            nsockets = self.remote.host.cluster.sockets_per_node
            ncores = self.remote.host.cluster.cores_per_socket
            memory = self.remote.host.cluster.memory

            mpi = True
            hyperthreading = self.remote.host.use_hyperthreading
            threads_per_core = self.remote.host.cluster.threads_per_core

        # Remote does not use a scheduling system
        else:

            # Get host properties
            nnodes = 1
            nsockets = int(math.floor(self.remote.free_sockets))
            ncores = self.remote.cores_per_socket
            memory = self.remote.free_memory

            mpi = True
            hyperthreading = self.remote.host.use_hyperthreading
            threads_per_core = self.remote.threads_per_core

        # Create the parallelization tool
        tool = ParallelizationTool()

        # Set configuration options
        tool.config.ski = self.definition.ski_path
        tool.config.input = self.definition.input_path

        # Set host properties
        tool.config.nnodes = nnodes
        tool.config.nsockets = nsockets
        tool.config.ncores = ncores
        tool.config.memory = memory

        # MPI available and used
        tool.config.mpi = mpi
        tool.config.hyperthreading = hyperthreading
        tool.config.threads_per_core = threads_per_core

        # Number of dust cells
        tool.config.ncells = None  # number of dust cells (relevant if ski file uses a tree dust grid)

        # Don't show the parallelization
        tool.config.show = False

        # Run the parallelization tool (passing the memory requirement of the simulation as an argument)
        tool.run(memory=self.memory)

        # Get the parallelization scheme
        self.parallelization = tool.parallelization

        # Debugging
        log.debug("The parallelization scheme is " + str(self.parallelization))

    # -----------------------------------------------------------------

    def check_parallelization(self):

        """
        This function checks whether the parallelization scheme that is asked by the user is possible given the
        number of cores and hyperthreads per core on the remote host.
        Returns:
        """

        # If the remote host uses a scheduling system, check whether the parallelization options are possible
        # based on the cluster properties defined in the configuration
        if self.remote.scheduler:

            # Determine the total number of hardware threads that can be used on the remote cluster
            hardware_threads_per_node = self.remote.cores_per_node
            if self.remote.use_hyperthreading: hardware_threads_per_node *= self.remote.threads_per_core

            # Raise an error if the number of requested threads per process exceeds the number of hardware threads
            # per node
            if self.config.arguments.parallel.threads > hardware_threads_per_node:
                raise RuntimeError("The number of requested threads per process exceeds the number of allowed threads per node")

            # Determine the number of processes per node (this same calculation is also done in JobScript)
            # self.remote.cores = cores per node
            processes_per_node = self.remote.cores_per_node // self.config.arguments.parallel.threads

            # Determine the amount of requested nodes based on the total number of processes and the number of processes per node
            requested_nodes = math.ceil(self.config.arguments.parallel.processes / processes_per_node)

            # Raise an error if the number of requested nodes exceeds the number of nodes of the system
            if requested_nodes > self.remote.nodes: raise RuntimeError("The required number of computing nodes for"
                                                                       "the requested number of processes and threads "
                                                                       "exceeds the existing number of nodes")

        # No scheduling system
        else:

            # Determine the total number of requested threads
            requested_threads = self.config.arguments.parallel.processes * self.config.arguments.parallel.threads

            # Determine the total number of hardware threads that can be used on the remote host
            hardware_threads = self.remote.cores_per_node
            if self.remote.use_hyperthreading: hardware_threads *= self.remote.threads_per_core

            # If the number of requested threads is greater than the allowed number of hardware threads, raise
            # an error
            if requested_threads > hardware_threads: raise RuntimeError("The requested number of processes and threads "
                                                                        "exceeds the total number of hardware threads")

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation ...")

        # Launch remotely or locally
        if self.config.remote is not None: self.launch_remote()
        else: self.launch_local()

    # -----------------------------------------------------------------

    def launch_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation locally...")

        # Run the simulation
        self.simulation = self.skirt.run(self.definition, logging_options=self.logging_options, silent=False, wait=True, progress_bar=self.config.progress_bar)

        # Set the simulation name
        self.simulation.name = self.definition.prefix

        # Set the analysis options for the simulation
        self.simulation.analysis = self.analysis_options

    # -----------------------------------------------------------------

    def launch_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation remotely ...")

        # Add the walltime to the scheduling options
        if self.config.walltime is not None:
            scheduling_options = SchedulingOptions()
            scheduling_options.walltime = self.config.walltime
        else: scheduling_options = None

        # Run the simulation
        self.simulation = self.remote.run(self.definition, self.logging_options, self.parallelization,
                                          scheduling_options=scheduling_options, attached=self.config.attached,
                                          analysis_options=self.analysis_options, progress_bar=self.config.progress_bar)

        # Set the analysis options for the simulation
        self.set_remote_simulation_options()

        # Save the simulation object
        self.simulation.save()

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Retrieving finished simulations...")

        # Get a list of the simulations that have been succesfully retrieved
        self.simulations = self.remote.retrieve()

    # -----------------------------------------------------------------

    def analyse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the output of retrieved simulations...")

        # Loop over the list of simulations and analyse them
        for simulation in self.simulations:

            # Run the analyser on the simulation
            self.analyser.run(simulation=simulation)

            # Clear the analyser
            self.analyser.clear()

    # -----------------------------------------------------------------

    def set_remote_simulation_options(self):

        """
        This function ...
        :return:
        """

        # Remove remote files
        self.simulation.remove_remote_input = not self.config.keep
        self.simulation.remove_remote_output = not self.config.keep
        self.simulation.remove_remote_simulation_directory = not self.config.keep

        # Retrieval
        self.simulation.retrieve_types = self.config.retrieve_types

# -----------------------------------------------------------------

class SimpleSKIRTLauncher(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # The SKIRT execution context
        self.skirt = SkirtExec()

    # -----------------------------------------------------------------

    def run(self, ski_path, out_path, wcs, total_flux, kernel, instrument_name=None):

        """
        This function ...
        :param ski_path:
        :param out_path:
        :param wcs:
        :param total_flux:
        :param kernel:
        :param instrument_name:
        :return:
        """

        from ...magic.core.frame import Frame
        from ..simulation.arguments import SkirtArguments

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True  # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running a SKIRT simulation with " + str(fs.name(ski_path)) + " ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Get the simulation prefix
        prefix = simulation.prefix()

        # Get the (frame)instrument name
        if instrument_name is None:

            # Get the name of the unique instrument (give an error if there are more instruments)
            instrument_names = simulation.parameters().get_instrument_names()
            assert len(instrument_names) == 1
            instrument_name = instrument_names[0]

        # Determine the name of the SKIRT output FITS file
        fits_name = prefix + "_" + instrument_name + "_total.fits"

        # Determine the path to the output FITS file
        fits_path = fs.join(out_path, fits_name)

        # Check if the output contains the "disk_earth_total.fits" file
        if not fs.is_file(fits_path): raise RuntimeError("Something went wrong with the " + prefix + " simulation: output FITS file missing")

        # Open the simulated frame
        simulated_frame = Frame.from_file(fits_path)

        # Set the coordinate system of the disk image
        simulated_frame.wcs = wcs

        # Debugging
        log.debug("Rescaling the " + prefix + " image to a flux density of " + str(total_flux) + " ...")

        # Rescale to the 3.6um flux density
        simulated_frame *= total_flux.value / simulated_frame.sum()
        simulated_frame.unit = total_flux.unit

        # Debugging
        log.debug("Convolving the " + prefix + " image ...")

        # Convolve the frame to the PACS 160 resolution
        simulated_frame.convolve(kernel)

        # Return the frame
        return simulated_frame

# -----------------------------------------------------------------
