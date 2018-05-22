#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.execute Executing the SKIRT command line application
#
# An instance of the SkirtExec class in this module represents a particular SKIRT executable,
# and allows invoking it with given command line arguments.

# -----------------------------------------------------------------

# Import standard modules
import os
import subprocess
import tempfile

# Import the relevant PTS classes and modules
from .arguments import SkirtArguments
from .definition import SimulationDefinition
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.log import log
from .definition import SingleSimulationDefinition
from .status import LogSimulationStatus, SpawnSimulationStatus
from ..tools import strings

# -----------------------------------------------------------------
#  SkirtExec class
# -----------------------------------------------------------------

# Convenience function
def run_simulation(definition, skirt_path=None, logging_options=None, parallelization=None, wait=True, silent=False,
                   show_progress=False, debug_output=False):

    """
    This function ...
    :param definition:
    :param skirt_path:
    :param logging_options:
    :param parallelization:
    :param wait:
    :param silent:
    :param show_progress:
    :param debug_output:
    :return:
    """

    # Create the SKIRT execution object
    skirt = SkirtExec(path=skirt_path)

    # Run the simulation
    simulation = skirt.run(definition, logging_options=logging_options, parallelization=parallelization, wait=wait,
                           silent=silent, show_progress=show_progress, debug_output=debug_output)

    # Return the simulation
    return simulation

# -----------------------------------------------------------------

## An instance of the SkirtExec class represents a particular SKIRT executable, and allows invoking it with
# given command line arguments. The current implementation only supports local execution; remote execution
# on another host will be added in the future.
class SkirtExec:

    ## The constructor accepts a single argument specifying the path to the <tt>SKIRT</tt> executable to be used.
    # - The path may or may not include the "skirt" filename and if not, it may or may not include the final "/" of a
    #   directory path. The path may be absolute, relative to a user's home folder, or relative to the current
    #   working directory.
    #   The default value is the empty path, which means <tt>SKIRT</tt> is looked for in the standard system path \c $PATH.
    # - mpistyle: the method or style to invoke the mpirun command; currently supported values are 'generic' (the
    #   default), which uses the standard -np switch useful for launching a number of MPI processes on a single
    #   computing node; 'lsf', which uses the -lsf switch supported by platform MPI under the LSF queueing system;
    #   and 'srun', which uses the -srun switch supported by platform MPI under the SLURM queueing system.
    def __init__(self, path="", mpi_style="generic"):

        # Set the SKIRT path
        self._path = path if path is not None else ""
        if not self._path.endswith("skirt"): self._path = os.path.join(self._path, "skirt")
        if self._path != "skirt": self._path = os.path.realpath(os.path.expanduser(self._path))
        if self._path == "skirt":
            if introspection.skirt_is_present(): self._path = introspection.skirt_path
            else: raise EnvironmentError("SKIRT is not installed or not in the PATH environment variable")

        # Indicate no simulations are running yet
        self._process = None

        # Set the MPI style
        self.mpi_style = mpi_style.lower()

        # Temporary files to fetch standard and error output
        self._output_file = None
        self._error_file = None

    @property
    def path(self):
        return self._path

    @property
    def has_pexpect(self):
        return introspection.is_present_package("pexpect")

    ## This function invokes the <tt>SKIRT</tt> executable with the simulations and command line options corresponding to the
    #  values of the function arguments:
    #
    # - skipattern: the paths/names of the ski file(s) to be executed
    # - recursive: if one or more ski filename patterns passed to the execute function contain wildcards,
    #   the recursive argument can be used to specify whether all directories recursively nested within the base
    #   path are searched as well, using the same filename pattern.
    # - inpath: a string specifying the absolute or relative path for simulation input files.
    # - outpath: a string specifying the absolute or relative path for simulation output files.
    # - skirel: if \c True, the simulation input/output paths are relative to the \em ski file being processed;
    #   if \c False or missing, they are relative to the current directory.
    # - threads: a positive integer specifying the number of parallel threads for each simulation; if zero or missing
    #   the number of logical cores on the computer is used.
    # - parallel: a positive integer specifying the number of simulations to be executed in parallel (by means of
    #   multithreading); if missing the default value is one.
    # - processes: the number of parallel MPI processes to be launched. The default value is one, which means MPI
    #   is not used. If the specified number of processes is larger than one, the value of \em parallel argument is
    #   ignored (i.e. you can't run multiple simulations in parallel when using MPI).
    # - dataparallel: if \c True and multiple processes are specified, the simulation will run in data-parallelized
    #   mode; if missing of False the simulation will run in task-parallelized mode.
    # - mpistyle: a string specifying the MPI implementation of the host system; if missing defaults to 'generic'.
    # - brief: if \c True, causes brief console logging; if missing defaults to False.
    # - verbose: This option has effect only if the number of processes is larger than one. If set to \c True, each
    #   process creates its own complete log file. If missing or set to \em False, only the root process creates a
    #   full log file, and the other processes only create a log file when there are errors or warnings.
    # - memory: This option enables memory logging for SKIRT.
    # - allocation: This option enables memory allocation logging for SKIRT.
    # - emulate: This option can be switched on to run SKIRT in 'emulation' mode to estimate the memory consumption.
    # - single: if \c True, only a single simulation is expected from the passed ski pattern, and the function returns
    #   a single SkirtSimulation object instead of a list.
    # - wait: if \c True or missing, the function waits until <tt>SKIRT</tt> execution completes and sends the brief
    #   SKIRT log to the standard console; if \c False the function returns immediately without waiting for SKIRT,
    #   and SKIRT's log messages are sent to the null device.
    # - silent: if \c True or missing, the SKIRT output will not be visible on the console.
    #
    # The function returns a list of SkirtSimulation instances corresponding to the simulations to be performed
    # (after processing any wildcards in the ski filenames), in arbitrary order.
    #
    def execute(self, skipattern, recursive=False, inpath=None, outpath=None, skirel=False,
                threads=0, parallel=1, processes=1, dataparallel=False, mpistyle='generic',
                brief=False, verbose=False, memory=False, allocation=False,
                emulate=False, single=False, wait=True, silent=False):

        self.mpi_style = mpistyle

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # The ski file pattern
        arguments.ski_pattern = skipattern
        arguments.recursive = recursive
        arguments.relative = skirel

        # Input and output
        arguments.input_path = inpath
        arguments.output_path = outpath

        # Parallelization settings
        arguments.parallel.threads = threads
        arguments.parallel.processes = processes
        arguments.parallel.simulations = parallel
        arguments.parallel.dataparallel = dataparallel

        # Logging settings
        arguments.logging.brief = brief
        arguments.logging.verbose = verbose
        arguments.logging.memory = memory
        arguments.logging.allocation = allocation

        # Other settings
        arguments.emulate = emulate
        arguments.single = single

        # Run SKIRT with the specified parameters
        return self.run(arguments, wait=wait, silent=silent)

    ## This function returns the MPI executable command
    @property
    def mpi_command(self):

        # Determine the MPI command based on MPI style
        if self.mpi_style == "lsf": mpi_command = "mpirun -lsf"
        elif self.mpi_style == "srun": mpi_command = "mpirun -srun"
        elif self.mpi_style == "generic": mpi_command = "mpirun"
        else: raise ValueError("Invalid MPI style")

    ## This property returns whether the target system uses a scheduling system
    @property
    def scheduler(self):

        # Based on MPI style
        if self.mpi_style == "lsf": return True
        elif self.mpi_style == "srun": return True
        elif self.mpi_style == "generic": return False
        else: raise ValueError("Invalid MPI style")

    ## This function does the same as the execute function, but obtains its arguments from a SkirtArguments object
    def run(self, definition_or_arguments, logging_options=None, parallelization=None, emulate=False, wait=True,
            silent=False, show_progress=False, finish_at=None, finish_after=None, debug_output=False):

        # Enable progress bar to true when finish_at or finish_after is defined because we'll have to follow up the simulation status
        if finish_at is not None or finish_after is not None: show_progress = True

        # Set the use_pexpect flag
        # Don't use pexpect when not waiting because I don't know a way to let pexpect run without doing child.expect(something)
        use_pexpect = self.has_pexpect and wait

        # The simulation names for different ski paths
        simulation_names = dict()

        # Simulation definition
        if isinstance(definition_or_arguments, SimulationDefinition):

            # The logging options cannot be None
            if logging_options is None: #raise ValueError("Logging options must be specified")
                log.warning("Logging options are not given: using default options ...")
                from ..launch.options import LoggingOptions
                logging_options = LoggingOptions()

            # Create the arguments
            arguments = SkirtArguments.from_definition(definition_or_arguments, logging_options, parallelization, emulate=emulate)

            # Set simulation name
            if isinstance(definition_or_arguments, SingleSimulationDefinition): simulation_names[definition_or_arguments.ski_path] = definition_or_arguments.name

        # Arguments are passed to the function
        elif isinstance(definition_or_arguments, SkirtArguments): arguments = definition_or_arguments
        else: raise ValueError("Invalid argument: should be simulation definition or SKIRT arguments instance")

        # Check whether MPI is present on this system if multiple processe are requested
        if arguments.parallel.processes > 1 and not introspection.has_mpi():
            log.warning("No mpirun executable: skipping simulations")
            return []

        # Get the command string
        command = arguments.to_command(self.scheduler, skirt_path=self._path, mpirun_path=self.mpi_command)

        # Not waiting and progress_bar don't go together!
        if show_progress and not wait: raise ValueError("Cannot show progress when 'wait' is False")
        if show_progress: wait = False

        # Debugging
        command_string = " ".join(command)
        log.debug("The command to launch SKIRT is: " + strings.add_other_quotes(command_string))

        # Using pexpect
        if use_pexpect: self.launch_pexpect(command_string, wait=wait, silent=silent)

        # Using subprocess
        else: self.launch_subprocess(command, wait=wait, silent=silent)

        # Show progress
        if show_progress: self.show_progress(arguments, finish_at=finish_at, finish_after=finish_after, debug_output=debug_output)
        elif self.using_pexpect: self.run_pexpect()

        # Return the list of simulations so that their results can be followed up
        simulations = arguments.simulations(simulation_names=simulation_names)

        # Check whether SKIRT has started
        if isinstance(self._process, subprocess.Popen):
            returncode = self._process.poll()
            has_finished = returncode is not None
        elif introspection.lazy_isinstance(self._process, "spawn", "pexpect", return_false_if_fail=True):
            has_finished = not self._process.isalive()
        elif self._process is None: has_finished = True
        else: raise RuntimeError("Unknown process handle")

        # when wait=True, or returncode is not None, SKIRT executable should have finished
        # so check whether indeed the log file(s) are present
        # otherwise something went wrong
        if wait or has_finished:

            # Check presence of log files
            if arguments.single:
                if not fs.is_file(simulations.logfilepath()): raise RuntimeError("SKIRT executable has stopped but log file is not present")
            else:
                for simulation in simulations:
                    if not fs.is_file(simulation.logfilepath()): raise RuntimeError("SKIRT executable has stopped but log file for simulation " + simulation.name + " is not present")

        # Return the list of simulations
        return simulations

    @property
    def using_pexpect(self):
        return introspection.lazy_isinstance(self._process, "spawn", "pexpect", return_false_if_fail=True)

    ## This function launches the SKIRT command using the Pexpect library
    def launch_pexpect(self, command_string, wait=True, silent=False):

        from ..tools.terminal import execute, launch_return

        if wait: execute(command_string, show_output=(not silent))
        else: self._process = launch_return(command_string)

        #from ..tools.terminal import launch_fetch_lines
        #for line in launch_fetch_lines(command_string): print(line)

    ## This function launches the SKIRT command using Python's subprocess
    def launch_subprocess(self, command, wait=True, silent=False):

        # Create a temporary file
        self._output_file = tempfile.TemporaryFile()
        self._error_file = tempfile.TemporaryFile()

        # Launch the SKIRT command
        if wait:

            self._process = None
            if silent: subprocess.call(command, stdout=self._output_file, stderr=self._error_file)
            else: subprocess.call(command)

        #else: self._process = subprocess.Popen(command, stdout=open(os.path.devnull, 'w'), stderr=subprocess.STDOUT)

        # CAUSES HANGING:
        # https://thraxil.org/users/anders/posts/2008/03/13/Subprocess-Hanging-PIPE-is-your-enemy/
        #else: self._process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1)

        # PROPER:
        else: self._process = subprocess.Popen(command, stdout=self._output_file, stderr=self._error_file)

    ## This function follows the progress of the simulation and
    def show_progress(self, arguments, finish_at=None, finish_after=None, debug_output=False):

        out_path = arguments.output_path if arguments.output_path is not None else fs.cwd()
        prefix = arguments.prefix
        log_path = fs.join(out_path, prefix + "_log.txt")

        ignore_output = ["Adding dust population", "Grain sizes range from", "Grain composition grid",
                         "Reading heat capacity data", "Reading grain composition", "closed.", "Reading SED data",
                         "Reading FITS file", "Frame dimensions:", "Writing grain size information", "created.",
                         "Writing optical dust population", "Writing dust population masses",
                         "Writing combined dust mix properties", "Reading wavelength grid data"]

        # Create the simulation status object
        if self.using_pexpect: status = SpawnSimulationStatus(self._process, debug_output=debug_output, ignore_output=ignore_output)
        else: status = LogSimulationStatus(log_path, debug_output=debug_output, ignore_output=ignore_output)

        # Show the simulation progress
        if self.using_pexpect: success = status.show_progress(finish_at=finish_at, finish_after=finish_after)
        else:
            with log.no_debugging(): success = status.show_progress(self._process, finish_at=finish_at, finish_after=finish_after)

        # Check whether not crashed
        if not success:

            # Show SKIRT error messages
            log.error("SKIRT output:")
            log.error("--------------------------")

            # Output was streamed to file
            if self._output_file is not None:
                for line in self._output_file: log.info(line)
                for line in self._error_file: log.error(line)
            else: self.show_errors(status)

            log.error("--------------------------")

            # Raise an error since the simulation crashed
            raise RuntimeError("The simulation crashed")

    ## This function shows the errors that appeared with the simulation
    def show_errors(self, status):

        # Subprocess
        if isinstance(self._process, subprocess.Popen):

            # Try to get standard and error output from the process
            out, err = self._process.communicate()

            for line in out:
                if "*** Error" in line:
                    line = line.split("*** Error: ")[1].split("\n")[0]
                    log.error(line)

            for line in err:
                if "*** Error" in line:
                    line = line.split("*** Error: ")[1].split("\n")[0]
                    log.error(line)

        # Pexpect spawn object
        elif introspection.lazy_isinstance(self._process, "spawn", "pexpect", return_false_if_fail=True):

            #print(self._process.logfile)
            #print(self._process.stdout)
            #for line in self._process.stdout: print(line)

            # Loop over the log lines
            for line in status.log_lines:
                if "*** Error" in line:
                    line = line.split("*** Error: ")[1].split("\n")[0]
                    log.error(line)

        # Invalid
        else: raise RuntimeError("Invalid state of the process object")

    ## This functions runs the simulation with Pexpect
    def run_pexpect(self):
        import pexpect
        self._process.expect(pexpect.EOF)

    ## This function returns True if the previously started SKIRT process is still running, False otherwise
    def isrunning(self):
        return (self._process != None and self._process.poll() is None)

    ## This function waits for the previously started SKIRT process to complete, if needed
    def wait(self):
        if self.isrunning(): self._process.wait()

    ## This function returns a string with version information on the SKIRT executable represented by this
    # object. The function invokes SKIRT with an incorrect command line argument to obtain this information.
    def version(self):
        # execute skirt with incorrect argument list and get its output
        process = subprocess.Popen([self._path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = process.communicate()[0]

        # return the relevant portion of the output
        return "SKIRT" + output.splitlines()[0].partition("SKIRT")[2]

    ## This function returns a string with the path of the SKIRT directory that is used, for example: HOME/SKIRT
    @property
    def root_directory(self):
        return os.path.dirname(os.path.dirname(os.path.dirname(self._path)))

    ## This function returns the path to the SKIRT repository directory (SKIRT/git)
    @property
    def repo_directory(self):
        return os.path.join(self.root_directory, "git")

    ## This function returns the path to the SKIRT run directory (SKIRT/run)
    @property
    def run_directory(self):
        return os.path.join(self.root_directory, "run")

# -----------------------------------------------------------------
