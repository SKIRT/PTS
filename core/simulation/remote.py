#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.remote Contains the SkirtRemote class, used for launching, checking and retrieving
#  remote SKIRT simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import tempfile

# Import the relevant PTS classes and modules
from ..basics.remote import Remote
from .jobscript import JobScript
from ..tools import time, inspection, filesystem
from .simulation import RemoteSimulation
from ..tools.logging import log
from ..launch.options import SchedulingOptions
from ..launch.parallelization import Parallelization

# -----------------------------------------------------------------

class SkirtRemote(Remote):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(SkirtRemote, self).__init__()

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

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster=None):

        """
        This function ...
        :param host_id:
        :param cluster:
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemote, self).setup(host_id, cluster)

        # Obtain some information about the SKIRT installation on the remote machine
        self.skirt_path = self.find_executable("skirt")

        # Check whether a SKIRT installation is found on the remote host
        if self.skirt_path is None: raise RuntimeError("SKIRT is not installed on the remote host")

        # We want absolute paths
        if self.skirt_path.startswith("~"):

            # Change the SKIRT path to the full, absolute path
            self.skirt_path = filesystem.join(self.home_directory, self.skirt_path[2:])

        # Determine the path to the SKIRT directory
        self.skirt_dir = self.skirt_path.split("/release")[0]

        # Determine the path to the SKIRT repository directory ('git')
        self.skirt_repo_dir = filesystem.join(self.skirt_dir, "git")

        # Determine the path to the SKIRT run directory
        self.skirt_run_dir = filesystem.join(self.skirt_dir, "run")

        # Determine the path to the local SKIRT run directory
        self.local_skirt_host_run_dir = filesystem.join(inspection.skirt_run_dir, self.host.id)

        # Create the local SKIRT run directory for this host if it doesn't already exist
        if not filesystem.is_directory(self.local_skirt_host_run_dir): filesystem.create_directory(self.local_skirt_host_run_dir, recursive=True)

        # Give a warning if the remote SKIRT version is different from the local SKIRT version
        local_version = inspection.skirt_version().split("built on")[0]
        remote_version = self.skirt_version.split("built on")[0]
        if remote_version != local_version:
            log.warning("Remote SKIRT version (" + remote_version + ") is different from local SKIRT version (" + local_version + ")")

    # -----------------------------------------------------------------

    def add_to_queue(self, arguments, name=None, scheduling_options=None, remote_input_path=None, analysis_options=None):

        """
        This function ...
        :param arguments:
        :param name: a name given to the simulation
        :param scheduling_options:
        :param remote_input_path:
        :param analysis_options:
        :return:
        """

        # Inform the user
        log.info("Adding simulation to the queue ...")

        # First create a copy of the arguments
        arguments = arguments.copy()

        # Create the remote simulation directory
        remote_simulation_path = self.create_simulation_directory(arguments)
        remote_simulation_name = filesystem.name(remote_simulation_path)

        # Set the name if none is given
        if name is None: name = remote_simulation_name

        # Make preparations for this simulation
        local_ski_path, local_input_path, local_output_path = self.prepare(arguments, remote_simulation_path, remote_input_path)

        # Add the SkirtArguments object to the queue
        self.queue.append((arguments, name))

        # If scheduling options are defined, add them to the dictionary
        if scheduling_options is not None: self.scheduling_options[name] = scheduling_options

        # Generate a new simulation ID based on the ID's currently in use
        simulation_id = self._new_simulation_id()

        # Create a simulation object
        simulation = self.create_simulation_object(arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path)

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

    def start_queue(self, screen_name=None, local_script_path=None, screen_output_path=None, group_simulations=False):

        """
        This function ...
        :param screen_name:
        :param local_script_path:
        :param screen_output_path:
        :param group_simulations:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Inform the user
        log.info("Starting the queued simulations remotely ...")

        # If the remote host uses a scheduling system, schedule all the simulations in the queue
        if self.scheduler: screen_name = self.start_queue_jobs(group_simulations)

        # Else, initiate a screen session in which the simulations are executed
        else: screen_name = self.start_queue_screen(screen_name, local_script_path, screen_output_path)

        # Clear the queue
        self.clear_queue()

        # Return the screen name
        return screen_name

    # -----------------------------------------------------------------

    def start_queue_jobs(self, group_simulations=False):

        """
        This function ...
        :param group_simulations:
        :return:
        """

        # Inform the user
        log.info("Starting the queue by scheduling the simulation as seperate jobs ...")

        # Group simulations in one job script
        if group_simulations: raise NotImplementedError("Group simulations is not implemented yet")

        # Don't group simulations
        else:

            # Loop over the items in the queue
            for arguments, name in self.queue:

                # Check whether scheduling options are defined for this simulation
                scheduling_options = self.scheduling_options[name] if name in self.scheduling_options else None

                # Submit the simulation to the remote scheduling system
                job_id = self.schedule(arguments, name, scheduling_options, local_ski_path=None, remote_simulation_path=None)

        # Return the screen name = None
        return None

    # -----------------------------------------------------------------

    def start_queue_screen(self, screen_name, local_script_path, screen_output_path):

        """
        This function ...
        :param screen_name:
        :param local_script_path:
        :param screen_output_path:
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

        # Write a general header to the batch script
        script_file.write("#!/bin/sh\n")
        script_file.write("# Batch script for running SKIRT on a remote system\n")
        script_file.write("\n")

        # Loop over the items in the queue
        for arguments, name in self.queue:

            print(arguments)

            # Write the command string to the job script
            threads_per_core = self.threads_per_core if self.use_hyperthreading else 1
            print("threasd per core", threads_per_core, type(threads_per_core))
            command = arguments.to_command(self.skirt_path, self.host.mpi_command, scheduler=False, bind_to_cores=self.host.force_process_binding, threads_per_core=threads_per_core, to_string=True)
            print("command:", command)
            script_file.write(command + "\n")

        # Write to disk
        script_file.flush()

        # Start a screen session
        self.start_screen(screen_name, local_script_path, self.skirt_run_dir, screen_output_path)

        # Close the script file (if it is temporary it will automatically be removed)
        script_file.close()

        # Return the screen name
        return screen_name

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

    def run(self, arguments, name=None, scheduling_options=None, analysis_options=None, local_script_path=None, screen_output_path=None):

        """
        This function ...
        :param arguments:
        :param name:
        :param scheduling_options:
        :param analysis_options:
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Raise an error if there are other simulations currently waiting in the queue
        if len(self.queue) > 0: raise RuntimeError("The simulation queue is not empty")

        # Add the simulation arguments to the queue
        simulation = self.add_to_queue(arguments, name, scheduling_options, analysis_options=analysis_options)

        # Start the queue
        screen_name = self.start_queue(name, local_script_path, screen_output_path)
        simulation.screen_name = screen_name

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def create_simulation_directory(self, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a unique name for the simulation directory
        skifile_name = filesystem.name(arguments.ski_pattern).split(".ski")[0]
        remote_simulation_name = time.unique_name(skifile_name, separator="__")

        # Determine the full path of the simulation directory on the remote system
        remote_simulation_path = filesystem.join(self.skirt_run_dir, remote_simulation_name)

        # Create the remote simulation directory
        self.create_directory(remote_simulation_path)

        # Return the path to the remote simulation directory
        return remote_simulation_path

    # -----------------------------------------------------------------

    def prepare(self, arguments, remote_simulation_path, remote_input_path=None):

        """
        This function ...
        :param arguments:
        :param remote_simulation_path:
        :param remote_input_path:
        :return:
        """

        # If an output path is defined in the remote host configuration file, use it for the simulation output
        if self.host.output_path is not None:

            # Get the name of the remote simulation directory and use use that name for the output directory
            remote_simulation_name = filesystem.name(remote_simulation_path)
            remote_output_path = filesystem.join(self.host.output_path, remote_simulation_name)

            # Expand the alias to the user's home directory
            remote_output_path = self.expand_user_path(remote_output_path)

            # If the remote output path is the same as the remote simulation path, use a folder called 'out' inside
            # the simulation directory instead for the output
            if remote_output_path == remote_simulation_path: remote_output_path = filesystem.join(remote_output_path, "out")

        # If an output path is not specified by the user, place a directory called 'out' next to the simulation's 'in' directory
        else: remote_output_path = filesystem.join(remote_simulation_path, "out")

        # Change the parameters to accomodate for the fact that we are running remotely
        # but store the paths to the local output directory because we want to copy the
        # results later
        local_input_path = arguments.input_path
        local_output_path = arguments.output_path

        # The simulation does not require input
        if local_input_path is None: remote_input_path = None

        # The simulation does require input
        else:

            # A remote input path is not specified, this means that we have yet to copy the input
            if remote_input_path is None:

                # Determine the full path to the input directory on the remote system
                remote_input_path = filesystem.join(remote_simulation_path, "in")

                # Create the remote input directory
                #self.create_directory(remote_input_path)

                # Copy the input directory to the remote host
                self.upload(local_input_path, remote_input_path)

            else:

                if not self.is_directory(remote_input_path): raise RuntimeError("The remote input directory does not exist")

        # Set the remote input and output path
        arguments.input_path = remote_input_path
        arguments.output_path = remote_output_path

        # Create the remote output directory
        self.create_directory(remote_output_path)

        # Set the remote ski file path
        local_ski_path = arguments.ski_pattern
        ski_name = filesystem.name(local_ski_path)
        remote_ski_path = filesystem.join(remote_simulation_path, ski_name)
        arguments.ski_pattern = remote_ski_path

        # Copy the input directory and the ski file to the remote host
        self.upload(local_ski_path, remote_simulation_path)

        # Return the paths of the local ski file and the local input and output directories
        return local_ski_path, local_input_path, local_output_path

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
        simulation_file_path = filesystem.join(self.local_skirt_host_run_dir, str(simulation_id) + ".sim")

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

    def schedule(self, arguments, name, scheduling_options, local_ski_path, remote_simulation_path):

        """
        This function ...
        :param arguments:
        :param name:
        :param scheduling_options:
        :param local_ski_path:
        :param remote_simulation_path:
        :return:
        """

        # Inform the suer
        log.info("Scheduling simulation '" + name + "' on the remote host ...")

        # Verify the scheduling options
        scheduling_options = self._verify_scheduling_options(scheduling_options, arguments, local_ski_path)

        # Now get the options
        nodes = scheduling_options.nodes
        ppn = scheduling_options.ppn
        mail = scheduling_options.mail
        full_node = scheduling_options.full_node
        walltime = scheduling_options.walltime
        local_jobscript_path = scheduling_options.local_jobscript_path

        # Create a job script next to the (local) simulation's ski file
        jobscript_name = filesystem.name(local_jobscript_path)
        jobscript = JobScript(local_jobscript_path, arguments, self.host.clusters[self.host.cluster_name],
                              self.skirt_path, self.host.mpi_command, self.host.modules, walltime, nodes, ppn,
                              name=name, mail=mail, full_node=full_node, bind_to_cores=self.host.force_process_binding)

        # Copy the job script to the remote simulation directory
        remote_simulation_path = filesystem.directory_of(arguments.ski_pattern) # NEW, to avoid having to pass this as an argument
        remote_jobscript_path = filesystem.join(remote_simulation_path, jobscript_name)
        self.upload(local_jobscript_path, remote_simulation_path)

        ## Swap clusters
        # Then, swap to the desired cluster and launch the job script
        #output = subprocess.check_output("module swap cluster/" + self._clustername + "; qsub " + self._path, shell=True, stderr=subprocess.STDOUT)

        # Submit the job script to the remote scheduling system
        #output = self.execute("qsub " + remote_jobscript_path, contains_extra_eof=True)
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
        command = arguments.to_command(self.skirt_path, self.host.mpi_command, self.scheduler, self.host.force_process_binding, threads_per_core=threads_per_core, to_string=True)
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
        for name in filesystem.files_in_path(self.local_skirt_host_run_dir, extension="sim", returns="name"):
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

                    # Download the simulation output
                    self.download(simulation.remote_output_path, simulation.output_path)

                # If retrieve file types are defined, download these files seperately to the local filesystem
                else:

                    # Create a list for the paths of the files that have to be copied to the local filesystem
                    copy_paths = []

                    # Loop over the files that are present in the remoute output directory
                    for filename in self.files_in_path(simulation.remote_output_path):

                        # Determine the full path to the output file
                        filepath = filesystem.join(simulation.remote_output_path, filename)

                        # Loop over the different possible file types and add the filepath if the particular type is in the list of types to retrieve
                        if filename.endswith("_ds_isrf.dat"):
                            if "isrf" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif "_ds_temp" in filename and filename.endswith(".fits"):
                            if "temp" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif filename.endswith("_sed.dat"):
                            if "sed" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif filename.endswith("_total.fits"):
                            if "image" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif filename.endswith("_ds_celltemps.dat"):
                            if "celltemp" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif "_log" in filename and filename.endswith(".txt"):
                            if "log" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif filename.endswith("_wavelengths.dat"):
                            if "wavelengths" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif "_ds_grid" in filename and filename.endswith(".dat"):
                            if "grid" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif "_ds_grho" in filename and filename.endswith(".fits"):
                            if "grho" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif "_ds_trho" in filename and filename.endswith(".fits"):
                            if "trho" in simulation.retrieve_types: copy_paths.append(filepath)
                        elif filename.endswith("_ds_convergence.dat"):
                            if "convergence" in simulation.retrieve_types: copy_paths.append(filepath)

                    # Debugging
                    log.debug("Retrieving files: " + str(copy_paths))
                    log.debug("Local output directory: " + simulation.output_path)

                    # Download the list of files to the local output directory
                    self.download(copy_paths, simulation.output_path)

                # If retrieval was succesful, add this information to the simulation file
                simulation.retrieved = True
                simulation.save()

                # Debug info
                log.debug("Successfully retrieved the necessary simulation output")

                # Remove the remote input, if present, if requested
                if simulation.remove_remote_input and simulation.has_input: self.remove_directory(simulation.remote_input_path)

                # Remove the remote output, if requested
                if simulation.remove_remote_output: self.remove_directory(simulation.remote_output_path)

                # If both the input and output directories have to be removed, the remote simulation directory
                # can be removed too
                if simulation.remove_remote_simulation_directory: self.remove_directory(simulation.remote_simulation_path)

                # Add the simulation to the list of retrieved simulations
                simulations.append(simulation)

        # Return the list of retrieved simulations
        return simulations

    # -----------------------------------------------------------------

    @property
    def skirt_version(self):

        """
        This function ...
        :return:
        """

        # Execute SKIRT with incorrect argument list and get its output
        output = self.execute("skirt --version")

        # Return the relevant portion of the output
        return "SKIRT" + output[0].partition("SKIRT")[2]

    # -----------------------------------------------------------------

    def update(self):

        """
        This function ...
        :return:
        """

        # Navigate to the SKIRT repository directory
        self.execute("cd " + self.skirt_repo_dir, output=False)

        # Update SKIRT
        self.execute("git pull origin master", output=False)

        # Compile the SKIRT code
        self.execute("./makeSKIRT.sh", output=False)

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
            for path in filesystem.files_in_path(self.local_skirt_host_run_dir, extension="sim", sort=int):

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # The name of the ski file (the simulation prefix)
                ski_name = simulation.prefix()

                # The path to the simulation log file
                remote_log_file_path = filesystem.join(simulation.remote_output_path, ski_name + "_log.txt")

                # Check whether the simulation has already been analysed
                if simulation.analysed: simulation_status = "analysed"

                # Check whether the simulation has already been retrieved
                elif simulation.retrieved: simulation_status = "retrieved"

                # Get the simulation status from the remote log file if not yet retrieved
                else: simulation_status = self.status_from_log_file(remote_log_file_path, simulation.screen_name, ski_name)

                # Add the simulation properties to the list
                entries.append((path, simulation_status))

        # If the remote has a scheduling system for launching jobs
        else:

            # Obtain job status information through the 'qstat' command
            output = self.execute("qstat")

            # Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
            queue_status = dict()

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
                    jobstatus = splitted_line[position-1]

                    # Add the status of this job to the dictionary
                    queue_status[jobid] = jobstatus

            # Search for files in the SKIRT run directory
            for path, name in filesystem.files_in_path(self.local_skirt_host_run_dir, extension="sim", returns=["path", "name"]):

                # Open the simulation file
                simulation = RemoteSimulation.from_file(path)

                # The name of the ski file (the simulation prefix)
                ski_name = simulation.prefix()

                # The path to the simulation log file
                remote_log_file_path = filesystem.join(simulation.remote_output_path, ski_name + "_log.txt")

                # Get the job ID from the name of the simulation file
                job_id = int(name)

                # Check if the simulation has already been analysed
                if simulation.analysed: simulation_status = "analysed"

                # Check if the simulation has already been retrieved
                elif simulation.retrieved: simulation_status = "retrieved"

                # Check if the job ID is in the list of queued or running jobs
                elif job_id in queue_status:

                    # Check the status of this simulation
                    job_status = queue_status[job_id]

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

                # If the simulation is not in the list of jobs
                else: simulation_status = self.status_from_log_file_job(remote_log_file_path, ski_name)

                # Add the simulation properties to the list
                entries.append((path, simulation_status))

        # Return the list of simulation properties
        return entries

    # -----------------------------------------------------------------

    def status_from_log_file(self, file_path, screen_name, simulation_prefix):

        """
        This function ...
        :param file_path:
        :param screen_name:
        :param simulation_prefix:
        :return:
        """

        # If the log file exists
        if self.is_file(file_path):

            # Get the last two lines of the remote log file
            output = self.execute("tail -2 " + file_path)

            # Get the last line of the actual simulation
            if " Available memory: " in output[1]: last = output[0]
            else: last = output[1]

            # Interpret the content of the last line
            if " Finished simulation " + simulation_prefix in last: simulation_status = "finished"
            elif " *** Error: " in last: simulation_status = "crashed"
            else:
                # The simulation is either still running or has been aborted
                if self.is_active_screen(screen_name): simulation_status = self.running_status_from_log_file(file_path)
                else: simulation_status = "aborted"

        # If the log file does not exist, the simulation has not started yet or has been cancelled
        else:

            # The simulation has not started or it's screen session has been cancelled
            if self.is_active_screen(screen_name): simulation_status = "queued"
            else: simulation_status = "cancelled"

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
            output = self.execute("tail -2 " + file_path)

            # Get the last line of the actual simulation
            if " Available memory: " in output[1]: last = output[0]
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

        output = self.read_text_file(file_path)

        phase = None
        cycle = None
        progress = None

        for line in output:

            if "Starting setup" in line: phase = "setup"
            elif "Starting the stellar emission phase" in line: phase = "stellar emission"
            elif "Launched stellar emission photon packages" in line:

                progress = float(line.split("packages: ")[1].split("%")[0])

            elif "Starting the first-stage dust self-absorption cycle" in line: phase = "self-absorption [stage 1"
            elif "Launched first-stage dust self-absorption cycle" in line:

                cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                progress = float(line.split("packages: ")[1].split("%")[0])

            elif "Starting the second-stage dust self-absorption cycle" in line: phase = "self-absorption [stage 2"
            elif "Launched second-stage dust self-absorption cycle" in line:

                cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                progress = float(line.split("packages: ")[1].split("%")[0])

            elif "Starting the last-stage dust self-absorption cycle" in line: phase = "self-absorption [stage 3"
            elif "Launched last-stage dust self-absorption cycle" in line:

                cycle = int(line.split("cycle ")[1].split(" photon packages")[0])
                progress = float(line.split("packages: ")[1].split("%")[0])

            elif "Starting the dust emission phase" in line: phase = "dust emission"
            elif "Launched dust emission photon packages" in line: progress = float(line.split("packages: ")[1].split("%")[0])
            elif "Starting writing results" in line: phase = "writing"

        if phase is None: return "running"
        elif "self-absorption" in phase: return "running: " + str(phase) + ", cycle " + str(cycle) + "] " + str(progress) + "%"
        elif "stellar emission" in phase or "dust emission" in phase: return "running: " + str(phase) + " " + str(progress) + "%"
        else: return "running: " + str(phase)

    # -----------------------------------------------------------------

    def _verify_scheduling_options(self, options, arguments, local_ski_path):

        """
        This function ...
        :param options:
        :param arguments:
        :param local_ski_path:
        :return:
        """

        # If scheduling options is not defined, create a new SchedulingOptions object
        if options is None: options = SchedulingOptions()

        # Test the presence of the 'nodes' and 'ppn' options
        if options.nodes is None or options.ppn is None:

            # Get the requirements in number of nodes and ppn
            processors = arguments.parallel.processes * arguments.parallel.threads
            nodes, ppn = self.get_requirements(processors)

            # Set the nodes and pppn
            options.nodes = nodes
            options.ppn = ppn

        # Check if 'mail' option is defined
        if options.mail is None: options.mail = False

        # Check if 'full_node' option is defined
        if options.full_node is None: options.full_node = True

        # We want to estimate the walltime here if it is not defined in the options
        if options.walltime is None:

            #factor = 1.2

            # Create and run a ResourceEstimator instance
            #estimator = ResourceEstimator()
            ##estimator.run(local_ski_path, arguments.parallel.processes, arguments.parallel.threads)
            #estimator.run(local_ski_path, 1, 1)

            # Return the estimated walltime
            ##walltime = estimator.walltime * factor
            #options.walltime = estimator.walltime_for(arguments.parallel.processes, arguments.parallel.threads) * factor

            raise RuntimeError("The walltime is not defined in the scheduling options")

        # Check if job script path is defined
        if options.local_jobscript_path is None:

            # Determine the jobscript path
            #local_simulation_path = filesystem.directory_of(local_ski_path)
            #local_jobscript_path = filesystem.join(local_simulation_path, "job.sh")
            local_jobscript_path = filesystem.join(filesystem.home(), time.unique_name("job") + ".sh")

            # Set the local jobscript path
            options.local_jobscript_path = local_jobscript_path

        # Return the scheduling options
        return options

# -----------------------------------------------------------------
