#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.remote Contains the SKIRTRemote class, used for launching, checking and retrieving
#  remote SKIRT simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import tempfile

# Import the relevant PTS classes and modules
from ..remote.remote import Remote
from .jobscript import JobScript, MultiJobScript, SKIRTJobScript
from ..tools import time, introspection
from ..tools import filesystem as fs
from .simulation import RemoteSimulation
from ..basics.log import log, no_debugging
from ..launch.options import SchedulingOptions
from ..simulation.parallelization import Parallelization
from ..simulation.arguments import SkirtArguments
from ..basics.handle import ExecutionHandle
from .status import LogSimulationStatus
from .input import SimulationInput
from ..tools import types
from .output import get_output_type, get_parent_type

# -----------------------------------------------------------------

def needs_retrieval(filename, types):

    """
    This function ...
    :param filename:
    :param types: retrieve types
    :return:
    """

    # Get the output type
    output_type = get_output_type(filename)

    # Check whether in specified retrieval types
    if output_type in types: return True

    # Get parent types
    parent_type = get_parent_type(output_type)

    # Check if (different from type itself) parent type is in the specified retrieval types
    if parent_type != output_type and parent_type in types: return True

    # Otherwise, return False
    return False

# -----------------------------------------------------------------

class SKIRTRemote(Remote):

    """
    This class ...
    """

    def __init__(self, host_id=None):

        """
        The constructor ...
        :param host_id:
        :return:
        """

        # Call the constructor of the base class
        super(SKIRTRemote, self).__init__()

        # -- Attributes --

        # Variables storing paths to the remote SKIRT installation location
        self.skirt_path = None
        self.skirt_dir = None
        self.skirt_repo_dir = None
        self.skirt_run_dir = None
        self.local_skirt_host_run_dir = None

        # Initialize an empty list for the simulation queue
        self.queue = []

        # Initialize a dictionary for the scheduling options
        self.scheduling_options = dict()

        # If host ID is given, setup
        if host_id is not None:
            if not self.setup(host_id): log.warning("The connection could not be made. Run setup().")

    # -----------------------------------------------------------------

    @classmethod
    def from_remote(cls, remote):

        """
        This function ...
        :param remote:
        :return:
        """

        # Create
        skirt = cls()

        # Give warning
        log.warning("When creating a SKIRTRemote instance from a regular Remote instance, the original Remote instance should not be used anymore")

        # Set attributes
        skirt.ssh = remote.ssh
        skirt.host = remote.host
        skirt.vpn = remote.vpn
        skirt.connected = remote.connected
        skirt.commands = remote.commands

        # Locate SKRIT
        success = skirt.locate_skirt()
        if not success: raise RuntimeError("Could not locate SKIRT on the remote host")

        # Return the instance
        return skirt

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster_name=None):

        """
        This function ...
        :param host_id:
        :param cluster_name:
        :return:
        """

        # Call the setup function of the base class
        success = super(SKIRTRemote, self).setup(host_id, cluster_name)
        if not success: return False

        # Locate SKIRT
        return self.locate_skirt()

    # -----------------------------------------------------------------

    def locate_skirt(self):

        """
        This function ...
        :return:
        """

        # Initialize
        success = True

        # Obtain some information about the SKIRT installation on the remote machine
        self.skirt_path = self.find_executable("skirt")

        # Check whether a SKIRT installation is found on the remote host
        if self.skirt_path is None: raise RuntimeError("SKIRT is not installed on the remote host")

        # We want absolute paths
        if self.skirt_path.startswith("~"):

            # Change the SKIRT path to the full, absolute path
            self.skirt_path = fs.join(self.home_directory, self.skirt_path[2:])

        # Determine the path to the SKIRT directory
        self.skirt_dir = self.skirt_path.split("/release")[0]

        # Determine the path to the SKIRT repository directory ('git')
        self.skirt_repo_dir = fs.join(self.skirt_dir, "git")

        # Determine the path to the SKIRT run directory
        self.skirt_run_dir = fs.join(self.skirt_dir, "run")

        # Determine the path to the local SKIRT run directory
        self.local_skirt_host_run_dir = fs.join(introspection.skirt_run_dir, self.host.id)

        # Create the local SKIRT run directory for this host if it doesn't already exist
        if not fs.is_directory(self.local_skirt_host_run_dir): fs.create_directory(self.local_skirt_host_run_dir, recursive=True)

        # Give a warning if the remote SKIRT version is different from the local SKIRT version
        local_version = introspection.skirt_version().split("built on")[0]
        remote_version = self.skirt_version.split("built on")[0]
        if remote_version != local_version:
            log.warning("Remote SKIRT version (" + remote_version + ") is different from local SKIRT version (" + local_version + ")")

        # Return succes
        return success

    # -----------------------------------------------------------------

    def add_to_queue(self, definition, logging_options, parallelization, name=None, scheduling_options=None,
                     remote_input_path=None, analysis_options=None, emulate=False, has_remote_input=False):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param name: a name given to the simulation
        :param scheduling_options:
        :param remote_input_path:
        :param analysis_options:
        :param emulate:
        :param has_remote_input:
        :return:
        """

        # Inform the user
        log.info("Adding simulation to the queue ...")

        # If a name is given for the simulation, check whether it doesn't contain spaces
        if name is not None and " " in name: raise ValueError("The simulation name cannot contain spaces")

        # Get local input and output path
        local_input_path = definition.input_path
        local_output_path = definition.output_path

        # Create the remote simulation directory
        remote_simulation_path = self.create_simulation_directory(definition)
        remote_simulation_name = fs.name(remote_simulation_path)

        # Set the name if none is given
        if name is None: name = remote_simulation_name

        # Make preparations for this simulation, create the SkirtArguments object
        arguments = self.prepare(definition, logging_options, parallelization, remote_simulation_path, remote_input_path, emulate=emulate, has_remote_input=has_remote_input)

        # Add the SkirtArguments object to the queue
        self.queue.append((arguments, name))

        # If scheduling options are defined, add them to the dictionary
        if scheduling_options is not None: self.scheduling_options[name] = scheduling_options

        # Generate a new simulation ID based on the ID's currently in use
        simulation_id = self._new_simulation_id()

        # Create a simulation object
        simulation = self.create_simulation_object(arguments, name, simulation_id, remote_simulation_path, definition.ski_path, local_input_path, local_output_path)

        # Set the parallelization properties to the simulation object
        processes = arguments.parallel.processes
        threads = arguments.parallel.threads
        threads_per_core = self.threads_per_core if self.use_hyperthreading else 1
        simulation.parallelization = Parallelization.from_processes_and_threads(processes, threads, threads_per_core)

        # If analysis options are defined, adjust the correponding settings of the simulation object
        if analysis_options is not None: simulation.set_analysis_options(analysis_options)

        # Save the simulation object
        simulation.save()

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def start_queue(self, queue_name=None, local_script_path=None, screen_output_path=None, group_simulations=False,
                    group_walltime=None, use_pts=False, jobscripts_path=None, attached=False, dry=False):

        """
        This function ...
        :param queue_name:
        :param local_script_path:
        :param screen_output_path:
        :param group_simulations:
        :param group_walltime:
        :param use_pts:
        :param jobscripts_path:
        :param attached:
        :param dry:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Inform the user
        log.info("Starting the queued simulations remotely ...")

        # Check whether attached mode is not requested for a scheduling remote
        if self.scheduler and attached: raise ValueError("Attached mode is not possible for a remote with scheduling system")

        # If the remote host uses a scheduling system, schedule all the simulations in the queue
        if self.scheduler: handles = self.start_queue_jobs(group_simulations=group_simulations,
                                                           group_walltime=group_walltime,
                                                           use_pts=use_pts,
                                                           jobscripts_path=jobscripts_path, queue_name=queue_name, dry=dry)

        # Else, initiate a screen session in which the simulations are executed
        else: handles = self.start_queue_screen(queue_name, local_script_path, screen_output_path, attached=attached, dry=dry)

        # Clear the queue
        self.clear_queue()

        # Return the execution handles
        return handles

    # -----------------------------------------------------------------

    def start_queue_jobs(self, group_simulations=False, group_walltime=None, use_pts=False, jobscripts_path=None,
                         queue_name=None, dry=False):

        """
        This function ...
        :param group_simulations:
        :param group_walltime:
        :param use_pts:
        :param jobscripts_path:
        :param queue_name:
        :param dry:
        :return:
        """

        # Inform the user
        log.info("Starting the queue by scheduling the simulations in one or multiple jobs ...")

        # Determine queue name
        if queue_name is None: queue_name = time.unique_name("simulations")

        # Group simulations in one job script
        if group_simulations:

            # Determine the preferred walltime per job
            preferred_walltime = group_walltime if group_walltime is not None else self.host.preferred_walltime

            # Use PTS or jobs
            if use_pts: handles = self._start_queue_jobs_pts(preferred_walltime, queue_name)
            else: handles = self._start_queue_jobs_groups(preferred_walltime, jobscripts_path)

        # Don't group simulations
        else:

            # Initialize a list to contain the execution handles
            handles = dict()

            # Loop over the items in the queue
            for arguments, name in self.queue:

                # Check whether scheduling options are defined for this simulation
                scheduling_options = self.scheduling_options[name] if name in self.scheduling_options else None

                # Submit the simulation to the remote scheduling system
                job_id = self.schedule(arguments, name, scheduling_options, local_ski_path=None, jobscript_dir_path=jobscripts_path, dry=dry)

                # Set the execution handle
                if job_id is not None: handle = ExecutionHandle.job(job_id, self.host_id)
                else: handle = ExecutionHandle.postponed(self.host_id)
                handles[name] = handle

        # Return execution handles
        return handles

    # -----------------------------------------------------------------

    def _start_queue_jobs_pts(self, preferred_walltime, queue_name):

        """
        This function ...
        :param preferred_walltime:
        :return:
        """

        # Initialize a list to contain the execution handles
        #handles = dict()

        # Start a remote python session
        session = self.start_python_session(assume_pts=True)

        # Initializing the simulation queue
        session.import_package("SimulationQueue", "pts.core.remote.queue")
        session.send_line("queue = SimulationQueue('" + queue_name + "')")
        session.send_line("queue.createtable()")
        #session.send_line("queue.close()")

        # Add the simulations
        with_statement = "with queue.transaction():"
        for_statement = "for record in records:"
        line = "db.insert(label, eaglesim, 0, record['galaxyid'], skitemplate)"
        lines = [line]

        # Execute the statements
        session.with_statement_and_loop(with_statement, for_statement, lines)

        # Close the queue
        session.send_line("queue.close()")

        # For each simulation:

        # update the records to indicate that a run has been scheduled and submit the job;
        # do this within a single transaction context to ensure that the scheduled job sees the updated records
        #print("Submitting job to queue", config.queue, "for run-ids", runids)
        #with db.transaction():
        #    db.updatestatus(runids, 'scheduled')
        #    db.updatefield(runids, 'queue', config.queue)
        #    subprocess.call(("bsub",), stdin=open(jobscriptname))

        #return handles

        # The execution handle is the same for all simulations (= the simulation queue database)
        handle = ExecutionHandle.sql(queue_name, self.host_id)

        # Return the execution handle
        return handle

    # -----------------------------------------------------------------

    def _start_queue_jobs_groups(self, preferred_walltime, jobscripts_path):

        """
        This function ...
        :param preferred_walltime:
        :param jobscripts_path:
        :return:
        """

        # Initialize a list to contain the execution handles
        handles = dict()

        # Keep track of the walltime
        current_walltime = 0.

        threads_for_job = None
        processes_for_job = None
        scheduling_options_for_job = SchedulingOptions()
        scheduling_options_for_job.walltime = preferred_walltime * 3600.  # in seconds
        simulations_for_job = []

        # Loop over the items in the queue
        for arguments, name in self.queue:

            # Get the estimated walltime
            estimated_walltime = self.scheduling_options[name].walltime

            # If this simulation doesn't fit in the current job anymore
            if current_walltime + estimated_walltime > preferred_walltime:

                # Schedule
                job_id = self.schedule_multisim(simulations_for_job, scheduling_options_for_job, jobscripts_path)
                handles[name] = ExecutionHandle.group_job(job_id, self.host_id)

                # Reset the current walltime
                current_walltime = 0.0

                # Reset threads and processes for job
                threads_for_job = None
                processes_for_job = None

                # Reset the scheduling options
                scheduling_options_for_job = SchedulingOptions()
                scheduling_options_for_job.walltime = preferred_walltime * 3600.  # in seconds

                # Reset the list of simulations for job
                simulations_for_job = []

            # If this simulation still fits in the current job
            else:

                # Check the scheduling options and parallelization
                threads = arguments.parallel.threads
                processes = arguments.parallel.processes

                if threads_for_job is None: threads_for_job = threads
                elif threads_for_job != threads: raise ValueError("Number of threads must be equal for all simulations in job")

                if processes_for_job is None: processes_for_job = processes
                elif processes_for_job != processes: raise ValueError("Number of processes must be equal for all simulations in job")

                # Add this simulation to the list for the current job
                simulations_for_job.append(arguments)

                # Update the current walltime
                current_walltime += estimated_walltime

        # Return the execution handles
        return handles

    # -----------------------------------------------------------------

    def start_queue_screen(self, screen_name, local_script_path, screen_output_path=None, attached=False, dry=False):

        """
        This function ...
        :param screen_name:
        :param local_script_path:
        :param screen_output_path:
        :param attached:
        :param dry:
        :return:
        """

        # Inform the user
        log.info("Starting the queue by initiating a screen session in which all simulations are executed ...")

        # Create a unique screen name indicating we are running SKIRT simulations if none is given
        if screen_name is None: screen_name = time.unique_name("SKIRT")

        # If the path for the shell script is not given, create a named temporary file
        if local_script_path is None:
            script_file = tempfile.NamedTemporaryFile()
            local_script_path = script_file.name

        # If a path is given, create a script file at the specified location
        else: script_file = open(local_script_path, 'w')

        # If no screen output path is set, create a directory
        if screen_output_path is None: screen_output_path = self.create_directory_in(self.pts_temp_path, screen_name)

        # Write a general header to the batch script
        remote_script_file_name = screen_name + ".sh"
        script_file.write("#!/bin/sh\n")
        script_file.write("# Batch script for running SKIRT on remote host " + self.host_id + "\n")
        script_file.write("# To execute manualy, upload this file to the remote filesystem in the following directory:\n")
        script_file.write("# " + screen_output_path + "\n")
        script_file.write("# under the name '" + remote_script_file_name + "' and enter the following commmands:\n")
        script_file.write("# cd '" + screen_output_path + "' # navigate to the screen output directory\n")
        script_file.write("# screen -S " + screen_name + " -L -d -m " + remote_script_file_name + "'\n")
        script_file.write("\n")

        # Loop over the items in the queue, add a line for each simulation
        for arguments, name in self.queue:

            # Write the command string to the job script
            threads_per_core = self.threads_per_core if self.use_hyperthreading else 1
            command = arguments.to_command(self.skirt_path, self.host.mpi_command, scheduler=False,
                                           bind_to_cores=self.host.force_process_binding,
                                           threads_per_core=threads_per_core, to_string=True, remote=self)
            script_file.write(command + "\n")

        # Write to disk
        script_file.flush()

        # Start a screen session, UNLESS DRY MODE IS ENABLED, IN WHICH CASE THE USER HAS TO UPLOAD AND RUN THE SCRIPT HIM/HERSELF
        if not dry: self.start_screen(screen_name, local_script_path, self.skirt_run_dir, screen_output_path, attached=attached)

        # Close the script file (if it is temporary it will automatically be removed)
        script_file.close()

        # Create the execution handle
        if attached: handle = ExecutionHandle.tty(self.tty, self.host_id)
        else: handle = ExecutionHandle.screen(screen_name, self.host_id, screen_output_path)

        # Return the execution handle(s)
        return handle

    # -----------------------------------------------------------------

    def clear_queue(self):

        """
        This function ...
        :return:
        """

        # Empty the queue and clear the scheduling options dictionary
        self.queue = []
        self.scheduling_options = dict()

    # -----------------------------------------------------------------

    def run(self, definition, logging_options, parallelization, name=None, scheduling_options=None,
            analysis_options=None, local_script_path=None, screen_output_path=None, attached=False,
            show_progress=False, has_remote_input=False):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param name:
        :param scheduling_options:
        :param analysis_options:
        :param local_script_path:
        :param screen_output_path:
        :param attached:
        :param show_progress:
        :param has_remote_input:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Raise an error if there are other simulations currently waiting in the queue
        if len(self.queue) > 0: raise RuntimeError("The simulation queue is not empty")

        # Add the simulation arguments to the queue
        simulation = self.add_to_queue(definition, logging_options, parallelization, name, scheduling_options, analysis_options=analysis_options, has_remote_input=has_remote_input)

        # Check whether attached mode is not requested for a scheduling remote
        if self.scheduler and attached: raise ValueError("Attached mode is not possible for a remote with scheduling system")

        # Progress bar: ask attached but set detached so that the start_queue function returns immediately
        if show_progress and not attached: raise ValueError("Cannot show progress when 'attached' is False")
        if show_progress: attached = False

        # Start the queue, get execution handle(s)
        handles = self.start_queue(name, local_script_path, screen_output_path, attached=attached)

        # Show progress bar with progress
        if show_progress:

            #out_path = arguments.output_path if arguments.output_path is not None else fs.cwd()
            #prefix = arguments.prefix
            #log_path = fs.join(out_path, prefix + "_log.txt")
            status = LogSimulationStatus(simulation.remote_log_file_path, remote=self)

            # Get the execution handle for the simulation
            handle = handles

            # Show the simulation progress
            with no_debugging(): success = status.show_progress(handle)

            # Check whether not crashed
            if not success: raise RuntimeError("The simulation crashed")

        # Set the execution handle for the simulation
        simulation.handle = handles if isinstance(handles, ExecutionHandle) else handles[0]
        simulation.save()

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def create_simulation_directory(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Create a unique name for the simulation directory
        skifile_name = fs.name(definition.ski_path).split(".ski")[0]
        remote_simulation_name = time.unique_name(skifile_name, separator="__")

        # Determine the full path of the simulation directory on the remote system
        remote_simulation_path = fs.join(self.skirt_run_dir, remote_simulation_name)

        # Create the remote simulation directory
        self.create_directory(remote_simulation_path)

        # Return the path to the remote simulation directory
        return remote_simulation_path

    # -----------------------------------------------------------------

    def prepare(self, definition, logging_options, parallelization, remote_simulation_path, remote_input_path=None,
                emulate=False, has_remote_input=False):

        """
        This function ...
        :param definition:
        :param logging_options:
        :param parallelization:
        :param remote_simulation_path:
        :param remote_input_path:
        :param emulate:
        :param has_remote_input:
        :return:
        """

        # Create the SkirtArguments object
        arguments = SkirtArguments(logging_options=logging_options, parallelization=parallelization, emulate=emulate)

        # If an output path is defined in the remote host configuration file, use it for the simulation output
        if self.host.output_path is not None:

            # Create the output directory
            output_path = self.absolute_path(self.host.output_path)
            if not self.is_directory(output_path): self.create_directory(output_path)

            # Get the name of the remote simulation directory and use use that name for the output directory
            remote_simulation_name = fs.name(remote_simulation_path)
            remote_output_path = fs.join(output_path, remote_simulation_name)

            # Expand the alias to the user's home directory
            #remote_output_path = self.absolute_path(remote_output_path)

            # If the remote output path is the same as the remote simulation path, use a folder called 'out' inside
            # the simulation directory instead for the output
            if remote_output_path == remote_simulation_path: remote_output_path = fs.join(remote_output_path, "out")

        # If an output path is not specified by the user, place a directory called 'out' next to the simulation's 'in' directory
        else: remote_output_path = fs.join(remote_simulation_path, "out")

        # The simulation does not require input
        if definition.input_path is None: remote_input_path = None

        # The simulation input is defined in terms of a single local directory
        elif isinstance(definition.input_path, basestring):

            # Check has_remote_input flag
            if has_remote_input: raise ValueError("Cannot enable 'has_remote_input' flag when input is a directory, only when input is defined in terms of a list or dictionary of filepaths (or SimulationInput object)")

            # A remote input path is not specified, this means that we have yet to copy the input
            if remote_input_path is None:

                # Determine the full path to the input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Copy the input directory to the remote host
                self.upload(definition.input_path, remote_input_path, show_output=True)

            # The specified remote input directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # If the simulation input is defined as a list of seperate file paths
        elif isinstance(definition.input_path, list):

            local_input_file_paths = definition.input_path # the list of file paths

            # A remote input path is not specified, this means that we have yet to copy the input files to a new
            # remote directory
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                self.create_directory(remote_input_path)

                # Check which files are actually local and which are already remote
                if has_remote_input:
                    nfiles = len(local_input_file_paths)
                    remote_input_file_paths = [filepath for filepath in local_input_file_paths if self.is_file(filepath)]
                    local_input_file_paths = [filepath for filepath in local_input_file_paths if fs.is_file(filepath)]
                    if len(remote_input_file_paths) + len(local_input_file_paths) != nfiles: raise ValueError("Some input files were not found")
                else: remote_input_file_paths = None

                # Upload the local input files to the new remote directory
                self.upload(local_input_file_paths, remote_input_path)

                # Copy the already remote files to the new remote directory
                if remote_input_file_paths is not None: self.copy_files(remote_input_file_paths, remote_input_path)

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # If we have a SimulationInput instance
        elif isinstance(definition.input_path, SimulationInput):

            # Check
            if has_remote_input: raise ValueError("Currently, has_remote_input is not possible with SimulationInput objects (due to internal checking in the latter class)")

            # If a remote input path is not specified
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                self.create_directory(remote_input_path)

                # Upload the local input files to the new remote directory
                for name, path in definition.input_path:

                    # Debugging
                    log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                    # Upload the file, giving it the desired name
                    self.upload(path, remote_input_path, new_name=name)

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # We have a dictionary
        elif types.is_dictionary(definition.input_path):

            # If a remote input path is not specified
            if remote_input_path is None:

                # Determine the full path to the remote input directory on the remote system
                remote_input_path = fs.join(remote_simulation_path, "in")

                # Create the remote directory
                self.create_directory(remote_input_path)

                # Upload the local input files to the new remote directory
                for name, path in definition.input_path.items():

                    # Check
                    if has_remote_input:

                        # Is local
                        if fs.is_file(path):

                            # Debugging
                            log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                            # Upload
                            self.upload(path, remote_input_path, new_name=name)

                        # Is remote
                        elif self.is_file(path):

                            # Debugging
                            log.debug("Copying the '" + path + "' file from '" + fs.directory_of(path) + "' to '" + remote_input_path + "' under the name '" + name + "'")

                            # Copy
                            self.copy_file(path, remote_input_path, new_name=name)

                        # Not found
                        else: raise ValueError("The input file '" + name + "' is not found remotely or locally")

                    # No checking
                    else:

                        # Debugging
                        log.debug("Uploading the '" + path + "' file to '" + remote_input_path + "' under the name '" + name + "'")

                        # Upload the file, giving it the desired name
                        self.upload(path, remote_input_path, new_name=name)

            # The specified remote directory (for re-usage of already uploaded input) does not exist
            else:
                if has_remote_input: raise ValueError("Cannot specify 'has_remote_input' and 'remote_input_path' simultaneously")
                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist: '" + remote_input_path + "'")

        # Invalid format for arguments.input_path
        else: raise ValueError("Invalid value for 'input_path': must be None, local directory path, list of file paths or SimulationInput object")

        # Create the remote output directory
        self.create_directory(remote_output_path)

        # Set the remote ski file path
        local_ski_path = definition.ski_path
        ski_name = fs.name(local_ski_path)
        remote_ski_path = fs.join(remote_simulation_path, ski_name)

        # Set the base simulation options such as ski path, input path and output path (remote)
        arguments.ski_pattern = remote_ski_path
        arguments.input_path = remote_input_path
        arguments.output_path = remote_output_path

        # Copy the input directory and the ski file to the remote host
        self.upload(local_ski_path, remote_simulation_path)

        # Return the SKIRT arguments instance
        return arguments

    # -----------------------------------------------------------------

    def create_simulation_object(self, arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path):

        """
        This function ...
        :param arguments:
        :param name:
        :param simulation_id:
        :param remote_simulation_path:
        :param local_ski_path:
        :param local_input_path:
        :param local_output_path:
        :return:
        """

        # Create a new remote simulation object
        simulation = RemoteSimulation(local_ski_path, local_input_path, local_output_path)

        # Determine and set the simulation file path
        simulation_file_path = fs.join(self.local_skirt_host_run_dir, str(simulation_id) + ".sim")

        # Set the host ID and cluster name (if applicable)
        simulation.host_id = self.host_id
        simulation.cluster_name = self.cluster_name

        # Set other attributes
        simulation.path = simulation_file_path
        simulation.id = simulation_id
        simulation.name = name
        simulation.remote_ski_path = arguments.ski_pattern
        simulation.remote_simulation_path = remote_simulation_path
        simulation.remote_input_path = arguments.input_path
        simulation.remote_output_path = arguments.output_path
        simulation.submitted_at = time.timestamp()

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def schedule(self, arguments, name, scheduling_options, local_ski_path, jobscript_dir_path=None, dry=False):

        """
        This function ...
        :param arguments:
        :param name:
        :param scheduling_options:
        :param local_ski_path:
        :param jobscript_dir_path:
        :param dry:
        :return:
        """

        # Inform the suer
        log.info("Scheduling simulation '" + name + "' on the remote host ...")

        if jobscript_dir_path is None:

            # Determine the jobscript path
            local_simulation_path = fs.directory_of(local_ski_path)
            jobscript_dir_path = local_simulation_path

        # Verify the scheduling options
        scheduling_options = self._verify_scheduling_options(scheduling_options, arguments, jobscript_dir_path)

        # Now get the options
        nodes = scheduling_options.nodes
        ppn = scheduling_options.ppn
        mail = scheduling_options.mail
        full_node = scheduling_options.full_node
        walltime = scheduling_options.walltime
        local_jobscript_path = scheduling_options.local_jobscript_path

        modules = []
        # module spider mympirun
        #modules.append("vsc-mympirun/3.4.3-intel-2016b-Python-2.7.12")
        modules.append("iimpi/2016b") # this loads GCC, icc, impi, python, iimpi, vsc-base, vsc-mympirun etc.

        # Create a job script next to the (local) simulation's ski file
        jobscript_name = fs.name(local_jobscript_path)
        #jobscript = JobScript(local_jobscript_path, arguments, self.host.cluster, self.skirt_path,
        #                      self.host.mpi_command, modules, walltime, nodes, ppn, name=name, mail=mail,
        #                      full_node=full_node, bind_to_cores=self.host.force_process_binding)

        # Determine remote jobscript location
        remote_simulation_path = fs.directory_of(arguments.ski_pattern)  # NEW, to avoid having to pass this as an argument
        remote_jobscript_path = fs.join(remote_simulation_path, jobscript_name)

        # Set header lines
        header_lines = []
        header_lines.append("To submit manualy, upload this file to the remote filesystem in the directory:")
        header_lines.append("'" + remote_simulation_path + "'")
        header_lines.append("with the name '" + jobscript_name + "' and enter the following commmands:")
        header_lines.append("cd '" + remote_simulation_path + "' # navigate to the simulation directory")
        header_lines.append("module swap cluster/[cluster_name] # swap to cluster [cluster_name], if desired")
        header_lines.append("qsub " + jobscript_name)
        header_lines.append("when the job has been submitted, copy the job ID that is returned [XXXX] and run the following command (locally):")
        header_lines.append("pts set_postponed_job_id " + self.host_id + " " + name + " XXXX")
        header_lines.append("")

        # Create the job script
        jobscript = SKIRTJobScript(name, arguments, self.host.cluster, self.skirt_path, self.host.mpi_command, walltime,
                                   modules, mail=mail, bind_to_cores=self.host.force_process_binding, extra_header_lines=header_lines)

        # Save the job script locally
        jobscript.saveto(local_jobscript_path)

        # Upload and submit
        if dry:
            log.warning("Dry mode is enabled, run job '" + jobscript_name + "' by locating the job script file at '" + local_jobscript_path + "' and following the instructions therein")
            job_id = None
        else: job_id = self.upload_and_submit_job(local_jobscript_path, remote_jobscript_path, remote_simulation_path)

        # Return the job ID
        return job_id

    # -----------------------------------------------------------------

    def upload_and_submit_job(self, local_jobscript_path, remote_jobscript_path, remote_simulation_path):

        """
        This function ...
        :return: 
        """

        # Debugging
        log.debug("Uploading and submitting job from script '" + remote_jobscript_path + "'")

        # Copy the job script to the remote simulation directory
        self.upload(local_jobscript_path, remote_simulation_path)

        ## Swap clusters
        # Then, swap to the desired cluster and launch the job script
        # output = subprocess.check_output("module swap cluster/" + self._clustername + "; qsub " + self._path, shell=True, stderr=subprocess.STDOUT)

        # Submit the job script to the remote scheduling system
        # output = self.execute("qsub " + remote_jobscript_path, contains_extra_eof=True)
        output = self.execute("qsub " + remote_jobscript_path)

        # The queue number of the submitted job is used to identify this simulation
        job_id = int(output[0].split(".")[0])
        return job_id

    # -----------------------------------------------------------------

    def schedule_multisim(self, arguments, scheduling_options, jobscripts_path):

        """
        This function ...
        :param arguments:
        :param scheduling_options:
        :return:
        """

        # Inform the suer
        log.info("Scheduling a job of " + str(len(arguments)) +  " simulations on the remote host ...")

        arguments_of_first_simulation = arguments[0] # for verifying the scheduling options (ALL SIMULATION HAVE THE SAME PARALLELIZATION MODE IN THIS MULTISIM FUNCTION)

        # Verify the scheduling options
        scheduling_options = self._verify_scheduling_options(scheduling_options, arguments_of_first_simulation, jobscripts_path)

        # Now get the options
        nodes = scheduling_options.nodes
        ppn = scheduling_options.ppn
        mail = scheduling_options.mail
        full_node = scheduling_options.full_node
        walltime = scheduling_options.walltime
        local_jobscript_path = scheduling_options.local_jobscript_path

        # Jobscript name
        jobscript_name = fs.name(local_jobscript_path)
        jobscript = MultiJobScript(local_jobscript_path, )

        # Copy the job script to the remote simulation directory
        remote_simulation_path = fs.directory_of(arguments.ski_pattern)  # NEW, to avoid having to pass this as an argument
        remote_jobscript_path = fs.join(remote_simulation_path, jobscript_name)
        self.upload(local_jobscript_path, remote_simulation_path)

        # SUBMIT THE JOB
        output = self.execute("qsub " + remote_jobscript_path)

        # The queue number of the submitted job is used to identify this simulation
        job_id = int(output[0].split(".")[0])

        # Return the job ID
        return job_id

    # -----------------------------------------------------------------

    def start(self, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Inform the user
        log.info("Starting simulation on the remote host")

        # Send the command to the remote machine using a screen session so that we can safely detach from the
        # remote shell
        threads_per_core = self.threads_per_core if self.use_hyperthreading else 1
        command = arguments.to_command(self.skirt_path, self.host.mpi_command, self.scheduler,
                                       self.host.force_process_binding, threads_per_core=threads_per_core,
                                       to_string=True, remote=self)
        self.execute("screen -d -m " + command, output=False)

        # Generate a new simulation ID based on the ID's currently in use
        simulation_id = self._new_simulation_id()

        # Return the simulation ID
        return simulation_id

    # -----------------------------------------------------------------

    def _simulation_ids_in_use(self):

        """
        This function ...
        :return:
        """

        # Check the contents of the local run directory to see which simulation id's are currently in use
        current_ids = []
        for name in fs.files_in_path(self.local_skirt_host_run_dir, extension="sim", returns="name"):

            # Get the simulation ID and add it to the list
            current_ids.append(int(name))

        # Return the list of currently used ID's
        return current_ids

    # -----------------------------------------------------------------

    def _new_simulation_id(self):

        """
        This function ...
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._simulation_ids_in_use()

        # Sort the current simulation ID's and find the lowest 'missing' integer number
        if len(current_ids) > 0:
            current_ids = sorted(current_ids)
            simulation_id = max(current_ids)+1
            for index in range(max(current_ids)):
                if current_ids[index] != index:
                    simulation_id = index
                    break

            # Return the simulation ID
            return simulation_id

        # If no simulation ID's are currently in use, return 0
        else: return 0

    # -----------------------------------------------------------------

    def _new_simulation_ids(self, count):

        """
        This function ...
        :param count:
        :return:
        """

        # Get a list of the ID's currently in use
        current_ids = self._simulation_ids_in_use()

        # Initialize a list to contain the new ID's
        new_ids = []

        # Sort the current simulation ID's and find the lowest 'missing' integer number
        if len(current_ids) > 0:
            current_ids = sorted(current_ids)
            for index in range(max(current_ids)):
                if current_ids[index] != index:
                    new_ids.append(index)
                    if len(new_ids) == count: return new_ids

            # Complement with new ID's
            max_id = max(new_ids)
            missing = count - len(new_ids)

            for index in range(max_id+1, max_id+1+missing):

                new_ids.append(index)

            return new_ids

        # If no simulation ID's are currently in use, return a list of the integers from 0 to count-1
        else: return range(count)

    # -----------------------------------------------------------------

    def retrieve(self):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Initialize a list to contain the simulations that have been retrieved
        simulations = []

        # Loop over the different entries of the status list
        for path, simulation_status in self.get_status():

            # Skip already retrieved and analysed simulations
            if simulation_status == "analysed": continue

            # If a simulation has been retrieved earlier, but is not yet analysed, also add it to the list of retrieved
            # simulations (again) so that its results can be analysed
            elif simulation_status == "retrieved":

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Add the retrieved simulation to the list
                simulations.append(simulation)

            # Finished simulations
            elif simulation_status == "finished":

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Debug info
                log.debug("Retrieving simulation " + str(simulation.name) + " with id " + str(simulation.id) + " ...")

                # If retrieve file types are not defined, download the complete output directory
                if simulation.retrieve_types is None or simulation.retrieve_types == "None":

                    # Debug info
                    log.debug("Retrieve file types are not defined, retrieving complete remote output directory ...")
                    log.debug("Local output directory: " + simulation.output_path)

                    # Check whether the output directory exists; if not, create it
                    if not fs.is_directory(simulation.output_path): fs.create_directory(simulation.output_path)

                    # Download the simulation output
                    self.download(simulation.remote_output_path, simulation.output_path)

                # If retrieve file types are defined, download these files seperately to the local fs
                else:

                    # Create a list for the paths of the files that have to be copied to the local fs
                    copy_paths = []

                    # Loop over the files that are present in the remoute output directory
                    for filename in self.files_in_path(simulation.remote_output_path):

                        # Determine the full path to the output file
                        filepath = fs.join(simulation.remote_output_path, filename)

                        # Check whether the file has to be retrieved
                        if needs_retrieval(filename, simulation.retrieve_types): copy_paths.append(filepath)

                    # Debugging
                    log.debug("Retrieving files: " + str(copy_paths))
                    log.debug("Local output directory: " + simulation.output_path)

                    # Check whether the output directory exists; if not, create it
                    if not fs.is_directory(simulation.output_path): fs.create_directory(simulation.output_path)

                    # Download the list of files to the local output directory
                    self.download(copy_paths, simulation.output_path)

                # If retrieval was succesful, add this information to the simulation file
                simulation.retrieved = True
                simulation.save()

                # Debug info
                log.debug("Successfully retrieved the necessary simulation output")

                # Remove the simulation from the remote
                simulation.remove_from_remote(self)

                # Add the simulation to the list of retrieved simulations
                simulations.append(simulation)

        # Return the list of retrieved simulations
        return simulations

    # -----------------------------------------------------------------

    def _get_simulation_status_not_scheduler(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # The name of the ski file (the simulation prefix)
        ski_name = simulation.prefix()

        # Determine the path to the remote log file
        remote_log_file_path = simulation.remote_log_file_path

        # Check whether the simulation has already been analysed
        if simulation.analysed: simulation_status = "analysed"

        # Check whether the simulation has already been retrieved
        elif simulation.retrieved: simulation_status = "retrieved"

        # Get the simulation status from the remote log file if not yet retrieved
        else: simulation_status = self.status_from_log_file(remote_log_file_path, simulation.handle, ski_name)

        # Return the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def _get_simulation_status_scheduler(self, simulation, jobs_status, python_session):

        """
        This function ...
        :param simulation:
        :param python_session:
        :return:
        """

        # The name of the ski file (the simulation prefix)
        ski_name = simulation.prefix()

        # The path to the simulation log file
        remote_log_file_path = simulation.remote_log_file_path

        # Check if the simulation has already been analysed
        if simulation.analysed: simulation_status = "analysed"

        # Check if the simulation has already been retrieved
        elif simulation.retrieved: simulation_status = "retrieved"

        # Simulation is a job
        elif simulation.handle.type == "job":

            # Get the job ID
            job_id = simulation.handle.value

            # Check if the job ID is in the list of queued or running jobs
            if job_id in jobs_status:

                # Check the status of this simulation
                job_status = jobs_status[job_id]

                # This simulation is still queued
                if job_status == 'Q': simulation_status = "queued"

                # This simulation is currently running
                elif job_status == 'R': simulation_status = self.running_status_from_log_file(remote_log_file_path)

                # If the job has been cancelled, check whether some part of the log file was already present
                # (the simulation was running but was aborted) or the log file is not present (the simulation is cancelled)
                elif job_status == "C":

                    if self.is_file(remote_log_file_path): simulation_status = "aborted"
                    else: simulation_status = "cancelled"

                # This simulation has an unknown status, check the log file
                else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

            # Job not present in the queue anymore: finished, crashed or aborted
            else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

        # Simulation is part of a group of simulations in a job
        elif simulation.handle.type == "group-job":

            # Get the job ID
            job_id = simulation.handle.value

            # Check if the job ID is in the list of queued or running jobs
            if job_id in jobs_status:

                # Check the status of the job
                job_status = jobs_status[job_id]

                # If the job is still queued, the simulation is also still queued
                if job_status == 'Q': simulation_status = "queued"

                # If the job is running
                elif job_status == 'R':

                    # Check if the log file exists
                    if self.is_file(remote_log_file_path):

                        # Get the last two lines of the remote log file
                        output = self.read_last_lines(remote_log_file_path, 2)

                        # Get the last line of the actual simulation
                        if len(output) == 0: return "invalid: cannot read log file" #simulation_status = "invalid: cannot read log file"
                        elif len(output) == 1: last = output[0]
                        elif " Available memory: " in output[1]: last = output[0]
                        else: last = output[1]

                        # Interpret the content of the last line
                        if " Finished simulation " + ski_name in last: simulation_status = "finished"
                        elif " *** Error: " in last: simulation_status = "crashed"
                        else: simulation_status = self.running_status_from_log_file(remote_log_file_path)

                    # The job is running but this simulation does not have a log file yet
                    else: simulation_status = "queued"

                # If the job has been cancelled
                elif job_status == 'C':

                    # Check if the log file exists
                    if self.is_file(remote_log_file_path): simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)
                    else: simulation_status = "cancelled"

                # This simulation has an unknown status, check the log file
                else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

        # Simulation is managed with an SQL database
        elif simulation.handle.type == "sql":

            # Get the simulation queue name
            queue_name = simulation.handle.value

            # Get the queue variable name in the remote python session
            queue_variable_name = "queue__" + queue_name

            # Get the status
            simulation_status = python_session.get_simple_property(queue_variable_name, "select('name=?', (simulation_name,))[0]['status']")

        # Invalid simulation handle
        else: raise RuntimeError("Unrecognized simulation handle")

        # Return the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def get_jobs_status(self):

        """
        This function ...
        :return:
        """

        # Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
        queue_status = dict()

        # Obtain job status information through the 'qstat' command
        output = self.execute("qstat")

        # Check every line in the output
        for line in output:

            # If this line mentions a job
            if "master15" in line:

                # Get the job ID
                jobid = int(line.split(".")[0])

                # Split the line
                splitted_line = line.split(" ")

                # Get the status (Q=queued, R=running)
                if "short" in splitted_line: position = splitted_line.index("short")
                elif "long" in splitted_line: position = splitted_line.index("long")
                else: continue
                jobstatus = splitted_line[position - 1]

                # Add the status of this job to the dictionary
                queue_status[jobid] = jobstatus

        # Return the queue status
        return queue_status

    # -----------------------------------------------------------------

    def get_status(self):

        """
        This function ..
        :return:
        """

        # Initialize a list to contain the statuses
        entries = []

        # If the remote host does not use a scheduling system
        if not self.scheduler:

            # Search for simulation files in the local SKIRT run/host_id directory
            for path in fs.files_in_path(self.local_skirt_host_run_dir, extension="sim", sort=int):

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Get the status
                simulation_status = self._get_simulation_status_not_scheduler(simulation)

                # Add the simulation properties to the list
                entries.append((path, simulation_status))

        # If the remote has a scheduling system for launching jobs
        else:

            # Get the status of the jobs
            jobs_status = self.get_jobs_status()

            # Start a remote python session
            if self.has_pts:

                # Start python session
                session = self.start_python_session(assume_pts=True)

                # Import statements
                session.import_package("SimulationQueue", from_name="pts.core.remote.queue")

                # Load all simulation queues
                for path, name in session.files_in_path(self.pts_run_path, extension="queue", returns=["path", "name"]):

                    # Load the queue
                    session.define_simple_variable("queue__" + name, "SimulationQueue.from_file('" + path + "'")

            # Do not open a session
            else: session = None

            # Search for simulation files in the SKIRT run directory
            for path, name in fs.files_in_path(self.local_skirt_host_run_dir, extension="sim", returns=["path", "name"]):

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # Get the status
                simulation_status = self._get_simulation_status_scheduler(simulation, jobs_status, session)

                # Add the simulation properties to the list
                entries.append((path, simulation_status))

        # Return the list of simulation properties
        return entries

    # -----------------------------------------------------------------

    def status_from_log_file(self, file_path, handle, simulation_prefix):

        """
        This function ...
        :param file_path:
        :param handle:
        :param simulation_prefix:
        :return:
        """

        # If the log file exists
        if self.is_file(file_path):

            # Get the last two lines of the remote log file
            output = self.read_last_lines(file_path, 2)

            # Get the last line of the actual simulation
            if len(output) == 0: return "invalid: cannot read log file"
            elif len(output) == 1: last = output[0]
            elif " Available memory: " in output[1]: last = output[0]
            else: last = output[1]

            # Interpret the content of the last line
            if " Finished simulation " + simulation_prefix in last: simulation_status = "finished"
            elif " *** Error: " in last: simulation_status = "crashed"
            else:

                # Screen session
                if handle.type == "screen":

                    screen_name = handle.value
                    if self.is_active_screen(screen_name): simulation_status = self.running_status_from_log_file(file_path)
                    else: simulation_status = "aborted"

                # Attached terminal session
                elif handle.type == "tty":

                    session_rank = handle.value
                    if session_rank in self.ttys: simulation_status = self.running_status_from_log_file(file_path)
                    else: simulation_status = "aborted"

                # Invalid execution handle
                else: raise ValueError("Invalid execution handle")

        # If the log file does not exist, the simulation has not started yet or has been cancelled
        else:

            # Screen session
            if handle.type == "screen":

                # The simulation has not started or it's screen session has been cancelled
                screen_name = handle.value
                if self.is_active_screen(screen_name): simulation_status = "queued"
                else: simulation_status = "cancelled"

            # Attached terminal session
            elif handle.type == "tty":

                session_rank = handle.value
                if session_rank in self.ttys: simulation_status = "queued"
                else: simulation_status = "cancelled"

            # Invalid execution handle
            else: raise ValueError("Invalid execution handle")

        # Return the string that indicates the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def status_from_log_file_job(self, file_path, simulation_prefix):

        """
        This function ...
        :param file_path:
        :param simulation_prefix:
        :return:
        """

        # Check whether the log file exists
        if self.is_file(file_path):

            # Get the last two lines of the remote log file
            output = self.read_last_lines(file_path, 2)

            # Get the last line of the actual simulation
            if len(output) == 0: return "invalid: cannot read log file"
            if len(output) == 1: last = output[0]
            elif " Available memory: " in output[1]: last = output[0]
            else: last = output[1]

            # Interpret the content of the last line
            if " Finished simulation " + simulation_prefix in last: simulation_status = "finished"
            elif " *** Error: " in last: simulation_status = "crashed"

            # The simulation cannot be running because we would have seen it in the qstat output
            # So with a partial log file, it must have been aborted
            else: simulation_status = "aborted"

        # If the log file does not exist, the simulation has been cancelled (if it would just be still scheduled
        # we would have encountered its job ID in the qstat output)
        else: simulation_status = "cancelled"

        # Return the string that indicates the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def running_status_from_log_file(self, file_path):

        """
        This function ...
        :return:
        """

        # Return string from simulation status
        return str(LogSimulationStatus(file_path, self))

    # -----------------------------------------------------------------

    def _verify_scheduling_options(self, options, arguments, jobscript_dir_path=None):

        """
        This function ...
        :param options:
        :param arguments:
        :param jobscript_dir_path:
        :return:
        """

        # If scheduling options is not defined, create a new SchedulingOptions object
        if options is None: options = SchedulingOptions()

        # Test the presence of the 'nodes' and 'ppn' options
        if options.nodes is None or options.ppn is None:

            # The number of threads per core that should be used
            threads_per_core = self.threads_per_core if self.use_hyperthreading else 1

            # Get the requirements in number of nodes and ppn
            processors = arguments.parallel.processes * arguments.parallel.threads / threads_per_core
            assert math.ceil(processors) == math.floor(processors) # make sure is integer
            processors = int(processors)
            nodes, ppn = self.get_requirements(processors)

            # Set the nodes and pppn
            options.nodes = nodes
            options.ppn = ppn

        # Check if 'mail' option is defined
        if options.mail is None: options.mail = False

        # Check if 'full_node' option is defined
        if options.full_node is None: options.full_node = True

        # We want to estimate the walltime here if it is not defined in the options
        if options.walltime is None: raise RuntimeError("The walltime is not defined in the scheduling options")

        # Check if job script path is defined
        if options.local_jobscript_path is None:

            # Check whether path to directory where jobscript should be saved locally is given
            if jobscript_dir_path is None: raise ValueError("Jobscript directory path must be specified in this function if the local jobscript path has not been set in the scheduling options")

            # Determine the job script path
            local_jobscript_path = fs.join(jobscript_dir_path, time.unique_name("job") + ".sh")

            # Set the local jobscript path
            options.local_jobscript_path = local_jobscript_path

        # Return the scheduling options
        return options

# -----------------------------------------------------------------
