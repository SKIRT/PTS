#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.arguments Contains the SkirtArguments class, used for representing the set of
#  command-line arguments that can be passed to SKIRT.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import fnmatch

# Import the relevant PTS classes and modules
from .simulation import SkirtSimulation
from ..basics.map import Map
from ..tools import filesystem as fs
from .definition import SingleSimulationDefinition, MultiSimulationDefinition
from ..basics.log import log
from .input import SimulationInput
from ..tools import types
from ..tools import strings
from ..tools.utils import lazyproperty
from ..tools import numbers
from ..tools import introspection

# -----------------------------------------------------------------

default_skirt_path = "skirt"
default_mpirun_path = "mpirun"

# -----------------------------------------------------------------

class SkirtArguments(object):

    """
    This class ...
    """

    def __init__(self, definition=None, logging_options=None, parallelization=None, emulate=False, skirt_path=None,
                 mpirun_path=None, simulation_name=None):

        """
        The constructor ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param skirt_path:
        :param mpirun_path:
        :param simulation_name:
        :return:
        """

        # TODO: discriminate between different types of 'SimulationDefinition' (multiple or single)

        # SKIRT path and MPIRUN path
        self.skirt_path = skirt_path if skirt_path is not None else default_skirt_path
        self.mpirun_path = mpirun_path if mpirun_path is not None else default_mpirun_path

        # Options for the ski file pattern
        self.ski_pattern = definition.ski_path if definition is not None else None
        self.recursive = None
        self.relative = None

        # The input and output paths
        self.input_path = definition.input_path if definition is not None else None
        self.output_path = definition.output_path if definition is not None else None

        # Other options
        self.emulate = emulate    # Run in emulation mode
        self.single = False       # True if only a single simulation is expected

        # Options for logging
        self.logging = Map()
        self.logging.brief = logging_options.brief if logging_options is not None else False     # Brief console logging
        self.logging.verbose = logging_options.verbose if logging_options is not None else False  # Verbose logging
        self.logging.memory = logging_options.memory if logging_options is not None else False  # State the amount of used memory with each log message
        self.logging.allocation = logging_options.allocation if logging_options is not None else False # Write log messages with the amount of (de)allocated memory
        self.logging.allocation_limit = logging_options.allocation_limit if logging_options is not None else 1e-5  # The lower limit for the amount of (de)allocated memory to be logged

        # Options for parallelization
        self.parallel = Map()
        self.parallel.simulations = None  # The number of parallel simulations
        self.parallel.threads = parallelization.threads if parallelization is not None else None # The number of parallel threads per simulation
        self.parallel.processes = parallelization.processes if parallelization is not None else None # The number of parallel processes per simulation
        self.parallel.dataparallel = parallelization.data_parallel if parallelization is not None else False  # Run in data parallelization mode
        self.parallel.threads_per_core = parallelization.threads_per_core if parallelization is not None else 1 # default is no hyperthreading

        # Simulation name (not required)
        self.simulation_name = simulation_name

    # -----------------------------------------------------------------

    @property
    def parallelization(self):

        """
        Thisf unction ...
        :return:
        """

        from .parallelization import Parallelization
        return Parallelization.from_processes_and_threads(self.parallel.processes, self.parallel.threads, threads_per_core=self.parallel.threads_per_core, data_parallel=self.parallel.dataparallel)

    # -----------------------------------------------------------------

    @parallelization.setter
    def parallelization(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.parallel.threads = value.threads
        self.parallel.processes = value.processes
        self.parallel.dataparallel = value.data_parallel
        self.parallel.threads_per_core = value.threads_per_core

    # -----------------------------------------------------------------

    @property
    def prefix(self):

        """
        This function ...
        :return:
        """

        if not fs.is_file(self.ski_pattern): raise RuntimeError("Cannot determine the prefix for the ski pattern '" + self.ski_pattern + "'. Does it define multiple files?")
        return fs.strip_extension(fs.name(self.ski_pattern))

    # -----------------------------------------------------------------

    @property
    def definition(self):

        """
        This function ...
        :return:
        """

        return SingleSimulationDefinition(self.ski_pattern, self.output_path, input_path=self.input_path, name=self.simulation_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_command(cls, command, nnodes=None):

        """
        This function ...
        :param command:
        :param nnodes:
        :return:
        """

        from ..tools import sequences

        # Get and remove comment from the line
        if "# " in command: command, comment = strings.split_at_last(command, "# ")
        else: comment = None

        # Split the command and get the first part
        parts = command.split()
        first = parts[0]

        # SKIRT is called directly, no MPI
        if first.endswith("skirt"):

            # Set SKIRT and mpirun path
            mpirun_path = None
            skirt_path = first
            nprocesses = 1
            cores_per_process = None

        # Probably MPI call
        else:

            # Set SKIRT and mpirun path
            mpirun_path = first
            skirt_path = sequences.find_unique_endswith(parts, "skirt")

            # Get number of processes
            if "-np" in command: nprocesses = int(command.split("-np")[1].split()[0])
            elif "-n" in command: nprocesses = int(command.split("-n")[1].split()[0])
            else: #raise ValueError("Unknown MPI command")
                if "--hybrid" in command:
                    #and nnodes is not None:
                    if nnodes is not None:
                        nprocesses_per_node = int(command.split("--hybrid ")[1].split(" ")[0])
                        nprocesses = nprocesses_per_node * nnodes
                    else:
                        log.warning("Cannot determine the number of processes: number of nodes is not known")
                        nprocesses = None
                else:
                    log.warning("Unknown MPI command: '" + command.split(skirt_path)[0].strip() + "': cannot determine the number of processes")
                    nprocesses = None

            # Get cores per process
            if "--map-by socket" in command: cores_per_process = int(command.split("--map-by socket:pe=")[1].split()[0])
            elif "--cpus-per-proc" in command: cores_per_process = int(command.split("--cpus-per-proc")[1].split()[0])
            else: cores_per_process = None

        # Find ski file path
        ski_path = sequences.find_unique_endswith(parts, ".ski")

        # Input and output
        output_path = strings.unquote(command.split(" -o ")[1].split()[0].strip()) if "-o" in parts else None
        input_path = strings.unquote(command.split(" -i ")[1].split()[0].strip()) if "-i" in parts else None

        # Other parallelization options
        nthreads = int(command.split(" -t ")[1].split()[0]) if "-t" in parts else None
        data_parallel = "-d" in parts

        # Find number of threads per core
        if cores_per_process is None: threads_per_core = 1
        else:
            threads_per_core = nthreads / cores_per_process
            if not numbers.is_integer(threads_per_core): raise RuntimeError("Something went wrong")
            threads_per_core = int(threads_per_core)

        # Logging options
        verbose = "-v" in parts
        brief = "-b" in parts
        memory = "-m" in parts
        allocation = "-l" in parts
        allocation_limit = float(command.split(" -l ")[1].split()[0]) if "-l" in parts else 1e-5

        # Emulate?
        emulate = "-e" in parts

        # Determine simulation name
        simulation_name = comment.strip() if comment is not None else None

        # Create and return the skirt arguments object
        return cls.single(ski_path, input_path, output_path, processes=nprocesses, threads=nthreads, brief=brief,
                          verbose=verbose, memory=memory, allocation=allocation, allocation_limit=allocation_limit,
                          data_parallel=data_parallel, threads_per_core=threads_per_core, emulate=emulate,
                          skirt_path=skirt_path, mpirun_path=mpirun_path, simulation_name=simulation_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_definition(cls, definition, logging_options=None, parallelization=None, emulate=False):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param emulate:
        :return:
        """

        # If the first argument defines a single simulation
        if isinstance(definition, SingleSimulationDefinition):

            # Create the SkirtArguments object
            arguments = SkirtArguments(logging_options=logging_options, parallelization=parallelization, simulation_name=definition.name)

            # Set the base simulation options such as ski path, input path and output path (remote)
            arguments.ski_pattern = definition.ski_path
            arguments.input_path = definition.input_path
            arguments.output_path = definition.output_path

            # Set other
            arguments.single = True
            arguments.emulate = emulate

            # Return the arguments object
            return arguments

        # If the first argument defines multiple simulations
        elif isinstance(definition, MultiSimulationDefinition): raise NotImplementedError("Not implemented yet")

        # Invalid first argument
        else: raise ValueError("Invalid argument for 'definition'")

    # -----------------------------------------------------------------

    @classmethod
    def single(cls, ski_path, input_path, output_path, processes=None, threads=None, brief=False, verbose=False,
               memory=False, allocation=False, allocation_limit=1e-5, data_parallel=False, threads_per_core=1,
               emulate=False, skirt_path=None, mpirun_path=None, simulation_name=None):

        """
        This function ...
        :param ski_path:
        :param input_path:
        :param output_path:
        :param processes:
        :param threads:
        :param verbose:
        :param memory:
        :param allocation:
        :param allocation_limit:
        :param data_parallel:
        :param threads_per_core:
        :param emulate:
        :param skirt_path:
        :param mpirun_path:
        :param simulation_name:
        :return:
        """

        # Create a SkirtArguments instance
        arguments = cls(skirt_path=skirt_path, mpirun_path=mpirun_path, simulation_name=simulation_name)

        # Set the options
        arguments.ski_pattern = ski_path
        arguments.single = True
        arguments.emulate = emulate
        arguments.recursive = False
        arguments.relative = False

        # Input and output
        arguments.input_path = input_path
        arguments.output_path = output_path

        # Logging
        arguments.logging.brief = brief
        arguments.logging.verbose = verbose
        arguments.logging.memory = memory
        arguments.logging.allocation = allocation
        arguments.logging.allocation_limit = allocation_limit

        # Parallelization
        arguments.parallel.processes = processes
        arguments.parallel.threads = threads
        arguments.parallel.dataparallel = data_parallel
        arguments.parallel.threads_per_core = threads_per_core

        # Return the new instance
        return arguments

    # -----------------------------------------------------------------

    def simulations(self, simulation_names=None, simulation_name=None):

        """
        This function ...
        :param simulation_names:
        :param simulation_name:
        :return:
        """

        # Initialize a list to contain the simulation objects
        simulations = []

        # Loop over the seperate ski files defined in the ski pattern
        pattern = [self.ski_pattern] if types.is_string_type(self.ski_pattern) else self.ski_pattern
        for skifile in pattern:

            # Determine the directory path and the actual file descriptor
            root, name = fs.directory_and_name(skifile)

            # Construct the 'dirlist' variable; this is a list of 3-tuples (dirpath, dirnames, filenames)
            if self.recursive: dirlist = os.walk(root)
            else: dirlist = [(root, [], filter(lambda fn: os.path.isfile(os.path.join(root,fn)), os.listdir(root)))]

            # Search for ski files matching the pattern and construct SkirtSimulation objects
            for dirpath, dirnames, filenames in dirlist:
                for filename in fnmatch.filter(filenames, name):

                    # Determine input and output path
                    inp = os.path.join(dirpath, self.input_path) if (self.relative and self.input_path is not None) else self.input_path
                    out = os.path.join(dirpath, self.output_path) if (self.relative and self.output_path is not None) else self.output_path

                    # Create the simulation and add it to the list
                    filepath = fs.join(dirpath, filename)
                    sim_name = simulation_names[filepath] if simulation_names is not None and filepath in simulation_names else None
                    simulations.append(SkirtSimulation(filename, inpath=inp, outpath=out, ski_path=filepath, name=sim_name))

        # Check whether the ski pattern is ought to represent only one particular simulation
        if self.single:

            # If multiple matching ski files are found, raise an error
            if len(simulations) > 1: raise ValueError("The specified ski pattern defines multiple simulations")
            else: simulation = simulations[0]

            # Set name and return simulation
            if simulation_name is not None: simulation.name = simulation_name
            elif self.simulation_name is not None: simulation.name = self.simulation_name
            return simulation

        # Else, just return the list of simulations (even when containing only one item)
        else:
            if simulation_name is not None: raise ValueError("'simulation_name' cannot be specified if single=False")
            return simulations

    # -----------------------------------------------------------------

    @lazyproperty
    def input_directory_path(self):

        """
        This function ...
        :return:
        """

        # Check input directory
        if isinstance(self.input_path, SimulationInput):
            input_dir_path = self.input_path.to_single_directory()

        # If the input consists of a list of paths, check whether they represent files in the same directory
        elif types.is_sequence(self.input_path):
            # print(1)
            input_dir_path = SimulationInput(*self.input_path).to_single_directory()
        elif types.is_string_type(self.input_path):
            # print(2)
            input_dir_path = SimulationInput(self.input_path).to_single_directory()
        elif types.is_dictionary(self.input_path):
            # print(3)
            # print(self.input_path)
            input_dir_path = SimulationInput(**self.input_path).to_single_directory()
        elif self.input_path is None:
            # print(4)
            input_dir_path = None
        else: raise ValueError("Type of simulation input not recognized")

        # Return
        return input_dir_path

    # -----------------------------------------------------------------

    def to_command(self, scheduler, skirt_path=None, mpirun_path=None, bind_to_cores=True, to_string=False, remote=None,
                   report_bindings=None):

        """
        This function ...
        :param scheduler:
        :param skirt_path:
        :param mpirun_path:
        :param bind_to_cores:
        :param to_string:
        :param remote:
        :param report_bindings:
        :return:
        """

        # Set SKIRT path and mpirun path
        skirt_path = skirt_path if skirt_path is not None else self.skirt_path
        mpirun_path = mpirun_path if mpirun_path is not None else self.mpirun_path

        # Create the argument list
        arguments = skirt_command(skirt_path, mpirun_path, bind_to_cores, self.parallel.processes, self.parallel.threads, self.parallel.threads_per_core, scheduler, remote, report_bindings=report_bindings)

        # Parallelization options
        if self.parallel.threads > 0: arguments += ["-t", str(self.parallel.threads)]
        if self.parallel.simulations > 1 and self.parallel.processes <= 1: arguments += ["-s", str(self.parallel.simulations)]
        if self.parallel.dataparallel and self.parallel.processes > 1: arguments += ["-d"]

        # Logging options
        if self.logging.brief: arguments += ["-b"]
        if self.logging.verbose: arguments += ["-v"]
        if self.logging.memory: arguments += ["-m"]
        if self.logging.allocation: arguments += ["-l", str(self.logging.allocation_limit)]

        # Options for input and output
        if self.input_directory_path is not None: arguments += ["-i", strings.add_quotes_if_spaces(self.input_directory_path)]
        if self.output_path is not None: arguments += ["-o", strings.add_quotes_if_spaces(self.output_path)]

        # Other options
        if self.emulate: arguments += ["-e"]

        # Ski file pattern
        if self.relative: arguments += ["-k"]
        if self.recursive: arguments += ["-r"]
        if types.is_string_type(self.ski_pattern): arguments += [strings.add_quotes_if_spaces(self.ski_pattern)]
        elif isinstance(self.ski_pattern, list): arguments += strings.add_quotes_if_spaces(self.ski_pattern)
        else: raise ValueError("The ski pattern must consist of either a string or a list of strings")

        # If requested, convert the argument list into a string
        # Create the final command string for this simulation
        if to_string: return " ".join(arguments)

        # Otherwise, return the list of argument values
        else: return arguments

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function creates a copy of this SkirtArguments object
        :return:
        """

        # Create a new SkirtArguments object
        arguments = SkirtArguments()

        ## Set options identical to this instance

        # SKIRT and MPI paths
        arguments.skirt_path = self.skirt_path
        arguments.mpirun_path = self.mpirun_path

        # Options for the ski file pattern
        arguments.ski_pattern = self.ski_pattern
        arguments.recursive = self.recursive
        arguments.relative = self.relative

        # The input and output paths
        arguments.input_path = self.input_path
        arguments.output_path = self.output_path

        # Other options
        arguments.emulate = self.emulate    # Run in emulation mode
        arguments.single = self.single     # True if only a single simulation is expected

        # Options for logging
        arguments.logging.brief = self.logging.brief            # Brief console logging
        arguments.logging.verbose = self.logging.verbose        # Verbose logging
        arguments.logging.memory = self.logging.memory          # State the amount of used memory with each log message
        arguments.logging.allocation = self.logging.allocation  # Write log messages with the amount of (de)allocated memory
        arguments.logging.allocation_limit = self.logging.allocation_limit  # The lower limit for the amount of (de)allocated memory to be logged

        # Options for parallelization
        arguments.parallel.simulations = self.parallel.simulations   # The number of parallel simulations
        arguments.parallel.threads = self.parallel.threads           # The number of parallel threads per simulation
        arguments.parallel.processes = self.parallel.processes       # The number of parallel processes per simulation
        arguments.parallel.dataparallel = self.parallel.dataparallel # Run in data parallelization mode
        arguments.parallel.threads_per_core = self.parallel.threads_per_core # Number of threads per core

        # Simulation name
        arguments.simulation_name = self.simulation_name

        # Return the new object
        return arguments

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        """

        properties = []
        properties.append("SKIRT command: " + self.skirt_path)
        properties.append("MPI command: " + self.mpirun_path)
        properties.append("ski path: " + self.ski_pattern)
        properties.append("recursive: " + str(self.recursive))
        properties.append("relative: " + str(self.relative))
        properties.append("input path: " + str(self.input_path))
        properties.append("output path: " + str(self.output_path))
        properties.append("emulate: " + str(self.emulate))
        properties.append("single: " + str(self.single))
        properties.append("brief: " + str(self.logging.brief))
        properties.append("verbose: " + str(self.logging.verbose))
        properties.append("memory: " + str(self.logging.memory))
        properties.append("allocation: " + str(self.logging.allocation))
        properties.append("allocation_limit: " + str(self.logging.allocation_limit))
        properties.append("simulations: " + str(self.parallel.simulations))
        properties.append("threads: " + str(self.parallel.threads))
        properties.append("processes: " + str(self.parallel.processes))
        properties.append("data-parallelization: " + str(self.parallel.dataparallel))
        properties.append("number of threads per core: " + str(self.parallel.threads_per_core))
        if self.simulation_name is not None: properties.append("simulation name: " + str(self.simulation_name))

        return_str = self.__class__.__name__ + ":\n"
        for property in properties: return_str += " -" + property + "\n"
        return return_str

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        """

        return '<' + self.__class__.__name__ + " ski path: '" + self.ski_pattern + "'>"

# -----------------------------------------------------------------

def skirt_command(skirt_path, mpi_command, bind_to_cores, processes, threads, threads_per_core, scheduler, remote=None,
                  report_bindings=None):

    """
    This function ...
    :param skirt_path:
    :param mpi_command:
    :param bind_to_cores:
    :param processes:
    :param threads:
    :param threads_per_core:
    :param scheduler:
    :param remote:
    :param report_bindings:
    :return:
    """

    from ..remote.remote import Remote
    from ..remote.host import Host

    # Set remote
    if remote is not None:
        # Load remote host if only the host ID is passed
        if types.is_string_type(remote): remote = Remote(host_id=remote)
        elif isinstance(remote, Host): remote = Remote()
        elif not isinstance(remote, Remote): raise ValueError("Invalid value for 'remote': should be remote instance, host instance or host ID string")

    # Multiprocessing mode
    if processes > 1:

        # Determine the command based on whether or not a scheduling system is used
        if scheduler: command = mpi_command.split()
        else: command = mpi_command.split() + ["-np", str(processes)]

        # If 'process to core' binding must be enabled, add the 'cpus-per-proc' option
        # (see https://www.open-mpi.org/faq/?category=tuning)
        if bind_to_cores:

            # TODO: SHOULD WE ALWAYS USE MAY BY SOCKET, OR DO WE HAVE TO CHECK HOW THE NUMBER OF PROCESSES COMPARES TO THE NUMBER OF USED SOCKETS?

            # Hyperthreading: threads_per_core will be > 1
            # No hyperthreading: threads_per_core will be 1
            # cores / process = (cores / thread) * (threads / process)
            cores_per_process = threads / threads_per_core
            if not numbers.is_integer(cores_per_process): raise RuntimeError("Something went wrong")
            cores_per_process = int(cores_per_process)

            # Remote is not specified: assume most recent version
            if remote is None:

                # Only add the map by and bind to options if we are not on the MacOS platform
                if not introspection.is_macos(): command += ["--map-by", "socket:pe=" + str(cores_per_process), "--bind-to", "core"]

            # Check if --map-by and --bind-to options are available
            elif remote.mpi_has_bind_to_option and remote.mpi_has_map_by_option: command += ["--map-by", "socket:pe=" + str(cores_per_process), "--bind-to", "core"]

            # Check if cpus-per-proc option is possible
            elif remote is None or remote.mpi_has_cpus_per_proc_option: command += ["--cpus-per-proc", str(cores_per_process)] # "CPU'S per process" means "core per process" in our definitions

            # Give warning
            else: log.warning("The MPI version on the remote host does not know the 'cpus-per-proc' command. Processes cannot be bound to cores")

        # Add report bindings option
        if remote is None:
            if report_bindings is None: report_bindings = not introspection.is_macos() # not on MacOS, else assume most recent version
            elif introspection.is_macos(): raise ValueError("Reporting the process binding is not available on MacOS")
        elif report_bindings is None: report_bindings = remote.mpi_has_report_bindings_option
        if report_bindings: command += ["--report-bindings"]

        # Add the SKIRT path and return the final command list
        command += [skirt_path]
        return command

    # Singleprocessing mode, no MPI command or options
    else: return [skirt_path]

# -----------------------------------------------------------------
