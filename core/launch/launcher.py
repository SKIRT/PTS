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
from ..simulation.remote import SKIRTRemote
from ..basics.log import log
from .options import SchedulingOptions
from ..advanced.parallelizationtool import determine_parallelization
from ..advanced.memoryestimator import estimate_memory
from ..simulation.parallelization import Parallelization, get_possible_nprocesses_in_memory
from .options import AnalysisOptions
from ..tools import filesystem as fs
from ..tools import parallelization
from ..simulation.skifile import SkiFile
from ..tools import formatting as fmt
from ..remote.host import load_host
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

class SKIRTLauncher(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SKIRTLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Create the local SKIRT execution context
        self.skirt = SkirtExec()

        # Create the SKIRT remote execution context
        self.remote = None

        # Create a SimulationAnalyser instance
        self.analyser = SimulationAnalyser()

        # The simulation definition
        self.definition = None
        self.remote_input_path = None
        self.has_remote_input_files = False

        # The logging options
        self.logging_options = None

        # The analysis options
        self.analysis_options = None

        # Scheduling options
        self.scheduling_options = None

        # The parallelization scheme
        self.parallelization = None

        # The number of processes
        self.nprocesses = None
        self.nprocesses_per_node = None

        # ADVANCED
        self.local_script_path = None
        self.screen_output_path = None

        # The simulation object
        self.simulation = None

        # Initialize a list to contain the retrieved finished simulations
        self.simulations = []

        # The specified number of cells
        self.ncells = None

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
        #return False
        return self.parallelization is not None

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        return self.config.remote

    # -----------------------------------------------------------------

    @property
    def host(self):

        """
        This function returns the Host object
        :return:
        """

        # If the setup has not been called yet, load the host
        if self.remote is None:
            if self.host_id is None: return None # local
            else: return load_host(self.host_id) # remote

        # If the setup has already been called
        else: return self.remote.host

    # -----------------------------------------------------------------

    @property
    def cluster_name(self):

        """
        This function ...
        :return:
        """

        # Local execution
        if self.host_id is None: return None

        # Remote, but setup has not been called yet
        elif self.remote is None: # setup has not been called

            # Check cluster_name configuration setting
            if self.config.cluster_name is not None: return self.config.cluster_name

            # Get default cluster for host
            else: return self.host.clusters.default

        # Remote, and setup has been called (remote has been setup)
        else: return self.remote.cluster_name

    # -----------------------------------------------------------------

    @property
    def uses_remote(self):

        """
        This function ...
        :return:
        """

        return self.host is not None

    # -----------------------------------------------------------------

    @property
    def uses_scheduler(self):

        """
        This function ...
        :return:
        """

        if self.host is None: return False
        else: return self.host.scheduler

    # -----------------------------------------------------------------

    @property
    def do_retrieve(self):

        """
        This function ...
        :return:
        """

        return self.config.retrieve and self.config.remote

    # -----------------------------------------------------------------

    @property
    def do_show(self):

        """
        This function ...
        :return:
        """

        return self.config.show

    # -----------------------------------------------------------------

    @property
    def do_analyse(self):

        """
        This function ...
        :return:
        """

        return self.config.analyse and self.has_simulations

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Set or check parallelization scheme
        self.set_or_check_parallelization()

        # 3. Launch the simulation
        self.launch()

        # 4. Retrieve the simulations that are finished
        if self.do_retrieve: self.retrieve()

        # 5. Show
        if self.do_show: self.show()

        # 6. Analyse the output of the retrieved simulations
        if self.do_analyse: self.analyse()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SKIRTLauncher, self).setup(**kwargs)

        # Check 'remote' input
        if "remote" in kwargs:
            remote = kwargs.pop("remote")
            if not remote.connected:
                if self.config.remote is None: raise ValueError("Unconnected remote is passed but host ID is not specified in configuration")
                remote.setup(self.config.remote, self.config.cluster_name)
            elif self.config.remote is not None and remote.host_id != self.config.remote:
                raise ValueError("Remote is passed for host '" + remote.host_id + "' but configured host ID is '" + self.config.remote + "'")
            self.remote = SKIRTRemote.from_remote(remote)

        # Setup the remote execution context
        if self.config.remote is not None:
            self.remote = SKIRTRemote()
            self.remote.setup(self.config.remote, self.config.cluster_name)

        # Create output directory
        if self.config.create_output and not fs.is_directory(self.config.simulation_output): fs.create_directory(self.config.simulation_output)

        # Create the logging options
        if "logging_options" in kwargs and kwargs["logging_options"] is not None: self.logging_options = kwargs.pop("logging_options")
        else:
            self.logging_options = LoggingOptions()
            self.logging_options.set_options(self.config.logging)

        # Create the analysis options
        if "analysis_options" in kwargs and kwargs["analysis_options"] is not None: self.set_analysis_options(kwargs.pop("analysis_options"))
        else: self.create_analysis_options()

        # Get scheduling options
        if self.uses_scheduler:
            if "scheduling_options" in kwargs and kwargs["scheduling_options"] is not None: self.scheduling_options = kwargs.pop("scheduling_options")
            # Add the walltime to the scheduling options
            if self.config.walltime is not None:
                if self.scheduling_options is None: self.scheduling_options = SchedulingOptions()
                self.scheduling_options.walltime = self.config.walltime

        # Get the memory information passed to this instance
        self.memory = kwargs.pop("memory", None)

        # Get the definition
        if "definition" in kwargs: self.definition = kwargs.pop("definition")

        # Has remote input?
        if "has_remote_input_files" in kwargs: self.has_remote_input_files = kwargs.pop("has_remote_input_files")
        if self.has_remote_input_files and self.remote is None: raise ValueError("Cannot have remote input files when launching simulation locally")

        # Has remote input path?
        if "remote_input_path" in kwargs: self.remote_input_path = kwargs.pop("remote_input_path")
        if self.remote_input_path is not None and self.remote is None: raise ValueError("Cannot have remote input path when launching simulation locally")
        if self.remote_input_path is not None and self.has_remote_input_files: raise ValueError("Cannot have remote input path and have seperate remote input files simultaneously")

        # Get the parallelization
        if "parallelization" in kwargs: self.parallelization = kwargs.pop("parallelization")

        # Get the number of processes
        if "nprocesses" in kwargs: self.nprocesses = kwargs.pop("nprocesses")
        if "nprocesses_per_node" in kwargs: self.nprocesses_per_node = kwargs.pop("nprocesses_per_node")

        # ADVANCED
        if "local_script_path" in kwargs: self.local_script_path = kwargs.pop("local_script_path")
        if "screen_output_path" in kwargs: self.screen_output_path = kwargs.pop("screen_output_path")

        # Get the number of dust cells if given
        if "ncells" in kwargs: self.ncells = kwargs.pop("ncells")

        ##

        # Create the simulation definition (if necessary)
        if self.definition is None: self.create_definition()

    # -----------------------------------------------------------------

    @property
    def has_nprocesses(self):

        """
        This function ...
        :return:
        """

        return self.nprocesses is not None

    # -----------------------------------------------------------------

    @property
    def has_nprocesses_per_node(self):

        """
        This function ...
        :return:
        """

        return self.nprocesses_per_node is not None

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
        self.analysis_options.check(logging_options=self.logging_options, output_path=self.config.simulation_output, retrieve_types=self.config.retrieve_types)

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
        self.analysis_options.check(logging_options=self.logging_options, output_path=self.config.simulation_output, retrieve_types=self.config.retrieve_types)

    # -----------------------------------------------------------------

    def create_definition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the simulation definition ...")

        # Create the simulation definition
        self.definition = SingleSimulationDefinition(self.config.ski, self.config.simulation_output, self.config.simulation_input)

    # -----------------------------------------------------------------

    def set_or_check_parallelization(self):

        """
        This function ...
        :return:
        """

        # Set the parallelization scheme
        if not self.has_parallelization: self.set_parallelization()

        # Check
        elif self.config.check_parallelization: self.check_parallelization()

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        # Remote
        if self.config.remote: self.set_parallelization_remote()

        # Local
        else: self.set_parallelization_local()

    # -----------------------------------------------------------------

    def set_parallelization_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the optimal parallelization scheme for local execution ...")

        # Determine the number of processes
        processes = self.get_nprocesses_local()

        # Calculate the maximum number of threads per process based on the current cpu load of the system
        free_cpus = monitoring.free_cpus()
        threads = int(free_cpus / processes)

        # If there are too little free cpus for the amount of processes, the number of threads will be smaller than one
        if threads < 1:
            log.warning("The number of processes was " + str(processes) + " but the number of free CPU's is only " + str(free_cpus))
            processes = max(int(free_cpus), 1)
            log.warning("Adjusting the number of processes to " + str(processes) + " ...")
            threads = 1

        # Determine number of cores
        cores = processes * threads
        threads_per_core = 2

        # Debugging
        log.debug("The number of cores is " + str(cores))
        log.debug("The number of thread per core is " + str(threads_per_core))
        log.debug("The number of processes is " + str(processes))

        # Set the parallelization scheme
        self.parallelization = Parallelization(cores, threads_per_core, processes, data_parallel=self.config.data_parallel_local)

        # Debugging
        log.debug("The parallelization scheme is " + str(self.parallelization))

    # -----------------------------------------------------------------

    def get_nprocesses_local(self):

        """
        This function ...
        :return:
        """

        # Check whether MPI is available
        if not introspection.has_mpi():

            # Check nprocesses
            if self.has_nprocesses and self.nprocesses > 1: raise ValueError("The number of processes that is specified is not possible: MPI installation not present")

            # Set number of processes to 1
            processes = 1

        # MPI present and number of processes is defined
        elif self.has_nprocesses: processes = self.nprocesses

        # MPI present and number of processes not defined
        else:

            # If memory requirement is not set
            if self.memory is None: self.memory = estimate_memory(self.definition.ski_path, input_path=self.config.simulation_input, ncells=self.ncells)

            # Determine the number of possible nprocesses
            processes = get_possible_nprocesses_in_memory(monitoring.free_memory(), self.memory.serial,
                                                          self.memory.parallel, data_parallel=self.config.data_parallel_local)

        # Return
        return processes

    # -----------------------------------------------------------------

    @lazyproperty
    def ski(self):

        """
        This function ...
        :return:
        """

        return SkiFile(self.definition.ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        # No file wavelength grid
        if not self.ski.wavelengthsfile(): return self.ski.nwavelengths()

        # File wavelength grid
        # Remote
        elif self.uses_remote:

            # Some input files may be remote
            if self.has_remote_input_files:

                from ..simulation.input import find_input_filepath
                filename = self.ski.wavelengthsfilename()
                filepath = find_input_filepath(filename, self.definition.input_path)

                if fs.is_file(filepath): return self.ski.nwavelengthsfile(self.definition.input_path)
                elif self.remote.is_file(filepath):
                    nwavelengths = int(self.remote.read_first_line(filepath))
                    return nwavelengths
                else: raise ValueError("We shouldn't get here")

            # Remote input directory is specified
            elif self.remote_input_path is not None:

                filename = self.ski.wavelengthsfilename()
                filepath = fs.join(self.remote_input_path, filename)

                # Check
                if not self.remote.is_file(filepath): raise IOError("The remote input file '" + filename + "' does not exist in '" + self.remote_input_path + "'")

                # Get the number of wavelengths and return
                nwavelengths = int(self.remote.read_first_line(filepath))
                return nwavelengths

            # Nothing is remote
            else: return self.ski.nwavelengthsfile(self.definition.input_path)

        # No remote
        else: return self.ski.nwavelengthsfile(self.definition.input_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dustlib_dimension(self):

        """
        This function ...
        :return:
        """

        return self.ski.dustlib_dimension()

    # -----------------------------------------------------------------

    def set_parallelization_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the optimal parallelization scheme for remote execution ...")

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

        # The number of processes is defined
        if self.has_nprocesses:

            # Determine cores per node and total number of cores
            cores_per_node = nsockets * ncores
            total_ncores = nnodes * cores_per_node

            # Check number of processes
            if self.nprocesses > cores_per_node: raise ValueError("The number of processes cannot be larger than the number of cores per node (" + str(cores_per_node) + ")")

            # Determine other parameters
            ppn = nsockets * ncores
            nprocesses_per_node = int(self.nprocesses / nnodes)
            nprocesses = nprocesses_per_node * nnodes
            ncores_per_process = ppn / nprocesses_per_node
            threads_per_core = threads_per_core if hyperthreading else 1
            threads_per_process = threads_per_core * ncores_per_process

            # Determine data-parallel flag
            if self.config.data_parallel_remote is None:

                if self.nwavelengths >= 10 * nprocesses and self.dustlib_dimension == 3: data_parallel = True
                else: data_parallel = False

            else: data_parallel = self.config.data_parallel_remote

            # Create the parallelization object
            self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
                                                             threads_per_process=threads_per_process,
                                                             data_parallel=data_parallel)

        # The number of processes per node is defined
        elif self.has_nprocesses_per_node:

            # Determine other parameters
            ppn = nsockets * ncores
            nprocesses = self.nprocesses_per_node * self.config.nnodes
            ncores_per_process = ppn / self.nprocesses_per_node
            threads_per_core = threads_per_core if hyperthreading else 1
            threads_per_process = threads_per_core * ncores_per_process
            total_ncores = nnodes * nsockets * ncores

            # Determine data-parallel flag
            if self.config.data_parallel_remote is None:

                if self.nwavelengths >= 10 * nprocesses and self.dustlib_dimension == 3: data_parallel = True
                else: data_parallel = False

            else: data_parallel = self.config.data_parallel_remote

            # Create the parallelization object
            self.parallelization = Parallelization.from_mode("hybrid", total_ncores, threads_per_core,
                                                             threads_per_process=threads_per_process,
                                                             data_parallel=data_parallel)

        # Determine the parallelization scheme with the parallelization tool
        # ski_path, input_path, memory, nnodes, nsockets, ncores, host_memory, mpi, hyperthreading, threads_per_core, ncells=None
        else: self.parallelization = determine_parallelization(self.definition.ski_path, self.definition.input_path, self.memory, nnodes, nsockets, ncores, memory, mpi, hyperthreading, threads_per_core, ncells=self.ncells, nwavelengths=self.nwavelengths)

        # Debugging
        log.debug("The parallelization scheme is " + str(self.parallelization))

    # -----------------------------------------------------------------

    def check_parallelization(self):

        """
        This function checks whether the parallelization scheme that is asked by the user is possible given the
        number of cores and hyperthreads per core on the remote host.
        Returns:
        """

        # Check locally or remotely
        if self.config.remote: self.check_parallelization_remote()
        else: self.check_parallelization_local()

    # -----------------------------------------------------------------

    def check_parallelization_remote(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Checking the parallelization scheme ...")

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
            processes_per_node = self.remote.cores_per_node // self.parallelization.threads

            # Determine the amount of requested nodes based on the total number of processes and the number of processes per node
            requested_nodes = math.ceil(self.config.arguments.parallel.processes / processes_per_node)

            # Raise an error if the number of requested nodes exceeds the number of nodes of the system
            if requested_nodes > self.remote.nodes: raise RuntimeError("The required number of computing nodes for"
                                                                       "the requested number of processes and threads "
                                                                       "exceeds the existing number of nodes")

        # No scheduling system
        else:

            # Determine the total number of requested threads
            requested_threads = self.parallelization.processes * self.parallelization.threads

            # Determine the total number of hardware threads that can be used on the remote host
            hardware_threads = self.remote.cores_per_node
            if self.remote.use_hyperthreading: hardware_threads *= self.remote.threads_per_core

            # If the number of requested threads is greater than the allowed number of hardware threads, raise
            # an error
            if requested_threads > hardware_threads: raise RuntimeError("The requested number of processes and threads "
                                                                        "exceeds the total number of hardware threads")

    # -----------------------------------------------------------------

    def check_parallelization_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the parallelization scheme ...")

        # Determine the total number of requested threads
        requested_threads = self.parallelization.nthreads

        # Determine the total number of hardware threads that can be used on the remote host
        hardware_threads = parallelization.ncores()
        if parallelization.has_hyperthreading(): hardware_threads *= parallelization.nthreads_per_core()

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

        # Launch remotely
        if self.config.remote is not None: self.launch_remote()

        # Launch locally
        else: self.launch_local()

    # -----------------------------------------------------------------

    def launch_local(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation locally ...")

        # Run the simulation
        self.simulation = self.skirt.run(self.definition, logging_options=self.logging_options, silent=False, wait=True,
                                         show_progress=self.config.show_progress, parallelization=self.parallelization,
                                         finish_after=self.config.finish_after, finish_at=self.config.finish_at,
                                         debug_output=self.config.debug_output)

        # Set the simulation name
        self.simulation.name = self.definition.prefix

        # Set the analysis options for the simulation
        self.simulation.analysis = self.analysis_options

        # Add the locally run simulation to the list of simulations to be analysed
        self.simulations.append(self.simulation)

    # -----------------------------------------------------------------

    @property
    def remove_remote_input(self):

        """
        This function ...
        :return:
        """

        if self.remote_input_path is not None or self.has_remote_input_files: return False
        elif self.config.keep_input: return False
        else: return not self.config.keep

    # -----------------------------------------------------------------

    @property
    def remove_remote_output(self):

        """
        Thisn function ...
        :return:
        """

        return not self.config.keep

    # -----------------------------------------------------------------

    @property
    def remove_remote_simulation_directory(self):

        """
        Thisn function ...
        :return:
        """

        return not self.config.keep

    # -----------------------------------------------------------------

    def launch_remote(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation remotely ...")

        # Resolve the remote screen output directory path
        screen_output_path = self.screen_output_path.replace("$HOME", self.remote.home_directory).replace("$SKIRT", self.remote.skirt_root_path) if self.screen_output_path is not None else None

        # Create the necessary directories for the screen output file
        #screen_output_dirpath = fs.directory_of(screen_output_path)
        #if not self.remote.is_directory(screen_output_dirpath): self.remote.create_directory(screen_output_dirpath, recursive=True)
        if screen_output_path is not None and not self.remote.is_directory(screen_output_path): self.remote.create_directory(screen_output_path, recursive=True)

        # Run the simulation
        self.simulation = self.remote.run(self.definition, self.logging_options, self.parallelization,
                                          scheduling_options=self.scheduling_options, attached=self.config.attached,
                                          analysis_options=self.analysis_options, show_progress=self.config.show_progress,
                                          local_script_path=self.local_script_path, screen_output_path=screen_output_path,
                                          remote_input_path=self.remote_input_path, has_remote_input=self.has_remote_input_files,
                                          debug_output=self.config.debug_output, retrieve_types=self.config.retrieve_types,
                                          remove_remote_input=self.remove_remote_input, remove_remote_output=self.remove_remote_output,
                                          remove_remote_simulation_directory=self.remove_remote_simulation_directory)

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

    @property
    def nsimulations(self):

        """
        This function ...
        :return:
        """

        return len(self.simulations)

    # -----------------------------------------------------------------

    @property
    def has_simulations(self):

        """
        This function ...
        :return:
        """

        return self.nsimulations > 0

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show finished simulations
        if self.has_simulations and self.config.show_finished: self.show_finished()

    # -----------------------------------------------------------------

    def show_finished(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Showing the output of finished simulations ...")

        # Loop over the simulations
        print("")
        for simulation in self.simulations:

            # Print the simulation name
            print(fmt.blue + simulation.prefix() + fmt.reset + ":")

            # Show the output
            simulation.output.show(line_prefix="  ")

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

class SingleImageSKIRTLauncher(object):

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

    def run(self, ski_path, out_path, wcs, total_flux, kernel, instrument_name=None, show_progress=False):

        """
        This function ...
        :param ski_path:
        :param out_path:
        :param wcs:
        :param total_flux:
        :param kernel:
        :param instrument_name:
        :param show_progress:
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
        simulation = self.skirt.run(arguments, silent=False if log.is_debug else True, show_progress=show_progress)

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

        # Rescale to the flux density
        simulated_frame.normalize(to=total_flux)

        # Debugging
        log.debug("Convolving the " + prefix + " image ...")

        # Convolve the frame
        simulated_frame.convolve(kernel)

        # Return the frame
        return simulated_frame

# -----------------------------------------------------------------

def launch_single_image():

    """
    This function ...
    :return:
    """

    pass

# -----------------------------------------------------------------

def generate_frame(ski, projection):

    """
    This function ...
    :param ski:
    :param projection:
    :return:
    """

    # Create instrument
    from ...modeling.basics.instruments import FrameInstrument
    instrument = FrameInstrument.from_projection(projection)

# -----------------------------------------------------------------
