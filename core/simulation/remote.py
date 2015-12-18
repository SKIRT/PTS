#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.remote Contains the SkirtRemote class, used for launching, checking and retreiving
#  remote SKIRT simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import re
import json
import pxssh
import pexpect
import tempfile

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from .jobscript import JobScript
from ..tools import time, inspection
from .simulation import SkirtSimulation
from ..tools import configuration
from ..basics.map import Map
from ..test.resources import ResourceEstimator

# -----------------------------------------------------------------

class SkirtRemote(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SkirtRemote, self).__init__(config, "core")

        ## Attributes

        # Create the SSH interface
        self.ssh = pxssh.pxssh()

        # Set the host configuration to None initially
        self.host = None

        # Set the connected flag to False initially
        self.connected = False

        # Set the SKIRT paths to None initially
        self.skirt_path = None
        self.skirt_dir = None
        self.skirt_run_dir = None
        self.local_skirt_run_dir = None

        # Generate a regular expression object to be used on the remote console output
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

        # Initialize an empty list for the simulation queue
        self.queue = []

    # -----------------------------------------------------------------

    def setup(self, host_id, cluster=None, pre_installation=False):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemote, self).setup()

        # Determine the path to the configuration file for the specified host and check if it is present
        host_file_path = os.path.join(inspection.pts_user_dir, "hosts", host_id + ".cfg")
        if not os.path.isfile(host_file_path): raise ValueError("The configuration settings for remote host " + host_id + " could not be found in the PTS/user/hosts directory")

        # Open the host configuration file
        self.host = configuration.open(host_file_path)

        # Set the host ID and cluster name (if a scheduling system is used)
        self.config.host_id = host_id
        if self.host.scheduler: # If no scheduling system is used, self.config.cluster_name stays at None
            if cluster is None:  # If no cluster name is given, use the default cluster (defined in host configuration file)
                self.config.cluster_name = self.host.clusters.default
            else: self.config.cluster_name = cluster

        # Make the connection
        self.login()

        # Load the necessary modules
        if self.host.scheduler:
            self.log.info("Loading necessary modules...")
            self.execute("module load " + " ".join(self.host.modules), output=False)

        # Check whether the output directory exists
        if not self.is_directory(self.host.output_path): raise ValueError("The specified output path does not exist")

        # Skip some steps in the setup when SKIRT has yet to be installed on the remote system
        if not pre_installation:

            # Obtain some information about the SKIRT installation on the remote machine
            output = self.execute("which skirt")

            # Determine the path to the SKIRT executable
            self.skirt_path = output[0]  # only one line is expected

            # We want absolute paths
            if self.skirt_path.startswith("~"):

                # Change the SKIRT path to the full, absolute path
                self.skirt_path = os.path.join(self.home_directory, self.skirt_path[2:])

            # Determine the path to the SKIRT directory
            self.skirt_dir = self.skirt_path.split("/release")[0]

            # Determine the path to the SKIRT run directory
            self.skirt_run_dir = os.path.join(self.skirt_dir, "run")

            # Determine the path to the local SKIRT run directory
            self.local_skirt_host_run_dir = os.path.join(inspection.skirt_run_dir, self.config.host_id)

            # Create the local SKIRT run directory for this host if it doesn't already exist
            if not os.path.isdir(self.local_skirt_host_run_dir): os.makedirs(self.local_skirt_host_run_dir)

            # Give a warning if the remote SKIRT version is different from the local SKIRT version
            local_version = inspection.skirt_version().split("built on")[0]
            remote_version = self.skirt_version.split("built on")[0]
            if remote_version != local_version: self.log.warning("Remote SKIRT version is different from local SKIRT version")

    # -----------------------------------------------------------------

    def login(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Logging in to the remote SKIRT environment on host '" + self.config.host_id + "'")

        # Connect to the remote host
        self.connected = self.ssh.login(self.host.name, self.host.user, self.host.password)

        # Check whether connection was succesful
        if not self.connected: raise RuntimeError("Connection failed")

    # -----------------------------------------------------------------

    def logout(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Logging out from the remote SKIRT environment")

        # Disconnect
        if self.connected: self.ssh.logout()

    # -----------------------------------------------------------------

    def add_to_queue(self, arguments, name=None, scheduling_options=None):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Inform the user
        self.log.info("Adding simulation to the queue...")

        # First create a copy of the arguments
        arguments = arguments.copy()

        # Create the remote simulation directory
        remote_simulation_path = self.create_simulation_directory(arguments)
        remote_simulation_name = os.path.basename(remote_simulation_path)

        # Set the name if none is given
        if name is None: name = remote_simulation_name

        # Make preparations for this simulation
        local_ski_path, local_input_path, local_output_path = self.prepare(arguments, remote_simulation_path)

        # If the remote host uses a scheduling system, submit the simulation right away
        if self.host.scheduler:

            # Submit the simulation to the remote scheduling system
            simulation_id = self.schedule(arguments, name, scheduling_options, local_ski_path, remote_simulation_path)

        # If no scheduling system is used, just store the SKIRT arguments in a list for now and execute the complete
        # list of simulations later on (when 'start_queue' is called)
        else:

            # Add the SkirtArguments object to the queue
            self.queue.append(arguments)

            # Generate a new simulation ID based on the ID's currently in use
            simulation_id = self._new_simulation_id()

        # Create a simulation file and return its path
        simulation_file_path = self.create_simulation_file(arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path)
        return simulation_file_path

    # -----------------------------------------------------------------

    def start_queue(self, screen_name=None, local_script_path=None, screen_output_path=None):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # If a scheduling system is used by the remote host, we don't need to do anything, simulations added to the queue
        # are already waiting to be executed (or are already being executed)
        if self.host.scheduler:
            self.log.warning("The remote host uses its own scheduling system so calling 'start_queue' will have no effect")
            return

        # Inform the user
        self.log.info("Starting the queued simulations remotely...")

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
        for arguments in self.queue:

            # Write the command string to the job script
            script_file.write(arguments.to_command(self.skirt_path, self.host.mpi_command, scheduler=False, to_string=True) + "\n")

        # Write to disk
        script_file.flush()

        # Copy the script to the remote host
        script_name = os.path.basename(local_script_path)
        remote_script_path = os.path.join(self.skirt_run_dir, script_name)
        self.upload(local_script_path, self.skirt_run_dir)

        # Close the script file (if it is temporary it will automatically be removed)
        script_file.close()

        # Make the shell script executable
        self.execute("chmod +x " + remote_script_path, output=False)

        # Create a unique screen name indicating we are running SKIRT simulations if none is given
        if screen_name is None: screen_name = time.unique_name("SKIRT")

        # Record the screen output: 'script' command
        if screen_output_path is not None: self.execute("script " + screen_output_path)

        # Create the screen session and execute the batch script
        self.execute("screen -S " + screen_name + " -d -m " + remote_script_path, output=False)

        # Remove the remote shell script
        self.execute("rm " + remote_script_path, output=False)

        # Clear the queue
        self.clear_queue()

        # Return the screen name
        return screen_name

    # -----------------------------------------------------------------

    def clear_queue(self):

        """
        This function ...
        :return:
        """

        self.queue = []

    # -----------------------------------------------------------------

    def kill_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Quit the specified screen session
        self.execute("screen -S " + name + " -X quit", output=False)

    # -----------------------------------------------------------------

    def kill_job(self, id):

        """
        This function ...
        :param id:
        :return:
        """

        # Stop the job with the specified ID
        self.execute("qdel " + str(id), output=False)

    # -----------------------------------------------------------------

    def screen_state(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Execute the 'screen -ls' command
        output = self.execute("screen -ls")

        # Loop over the different active screen sessions
        for entry in output:

            # Check if the specified screen name corresponds to the current entry
            if name in entry:

                # Check the state of this screen session
                if "Detached" in entry: return "detached"
                elif "Attached" in entry: return "attached"
                else: raise ValueError("Screen " + name + " in unkown state")

        # If the screen name was not encountered, return None (the screen session has finished or has been aborted)
        return None

    # -----------------------------------------------------------------

    def is_active_screen(self, name):

        """
        This function ...
        :return:
        """

        state = self.screen_state(name)
        return state == "detached" or state == "attached"

    # -----------------------------------------------------------------

    def is_attached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "attached"

    # -----------------------------------------------------------------

    def is_detached_screen(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.screen_state(name) == "detached"

    # -----------------------------------------------------------------

    def run(self, arguments, name=None, scheduling_options=None):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Raise an error if there are other simulations currently waiting in the queue
        if len(self.queue) > 0: raise RuntimeError("The simulation queue is not empty")

        # Add the simulation arguments to the queue
        simulation_file_path = self.add_to_queue(arguments, name, scheduling_options)

        # Start the queue if that is not left up to the remote's own scheduling system
        if not self.host.scheduler:
            screen_name = self.start_queue(name)
            with open(simulation_file_path, 'a') as simulation_file:
                simulation_file.write("launched within screen session " + screen_name + "\n")

        # Return the path to the simulation file
        return simulation_file_path

    # -----------------------------------------------------------------

    def create_simulation_directory(self, arguments):

        """
        This function ...
        :return:
        """

        # Create a unique name for the simulation directory
        skifile_name = os.path.basename(arguments.ski_pattern).split(".ski")[0]
        remote_simulation_name = time.unique_name(skifile_name, separator="__")

        # Determine the full path of the simulation directory on the remote system
        remote_simulation_path = os.path.join(self.skirt_run_dir, remote_simulation_name)

        # Create the remote simulation directory
        self.execute("mkdir " + remote_simulation_path, output=False)

        # Return the path to the remote simulation directory
        return remote_simulation_path

    # -----------------------------------------------------------------

    def prepare(self, arguments, remote_simulation_path):

        """
        This function ...
        :return:
        """

        # Determine the full paths to the input and output directories on the remote system
        remote_input_path = os.path.join(remote_simulation_path, "in")

        # If an output path is defined in the remote host configuration file, use it for the simulation output
        if self.host.output_path is not None:

            # Get the name of the remote simulation directory and use use that name for the output directory
            remote_simulation_name = os.path.basename(remote_simulation_path)
            remote_output_path = os.path.join(self.host.output_path, remote_simulation_name)

        # If an output path is not specified by the user, place a directory called 'out' next to the simulation's 'in' directory
        else: remote_output_path = os.path.join(remote_simulation_path, "out")

        # Change the parameters to accomodate for the fact that we are running remotely
        # but store the paths to the local output directory because we want to copy the
        # results later
        local_input_path = arguments.input_path
        local_output_path = arguments.output_path

        if local_input_path is None: remote_input_path = None

        arguments.input_path = remote_input_path
        arguments.output_path = remote_output_path

        # Create the remote input directory if necessary
        if remote_input_path is not None: self.execute("mkdir " + remote_output_path, output=False)

        # Create the remote output directory
        self.execute("mkdir " + remote_output_path, output=False)

        local_ski_path = arguments.ski_pattern
        ski_name = os.path.basename(local_ski_path)
        remote_ski_path = os.path.join(remote_simulation_path, ski_name)
        arguments.ski_pattern = remote_ski_path

        # Copy the input directory and the ski file to the remote host
        self.upload(local_ski_path, remote_simulation_path)
        if local_input_path is not None: self.upload(local_input_path, remote_input_path)

        # Return the paths of the local ski file and the local input and output directories
        return local_ski_path, local_input_path, local_output_path

    # -----------------------------------------------------------------

    def create_simulation_file(self, arguments, name, simulation_id, remote_simulation_path, local_ski_path, local_input_path, local_output_path):

        """
        This function ...
        :return:
        """

        # Create the simulation file
        simulation_file_path = os.path.join(self.local_skirt_host_run_dir, str(simulation_id) + ".sim")
        simulation_file = open(simulation_file_path, 'w')

        # Determine the simulation name from the ski file name if none is given
        #if name is None: name = os.path.basename(arguments.ski_pattern).split(".")[0]

        # Add the simulation information
        simulation_file.write("simulation name: " + name + "\n")
        simulation_file.write("local skifile path: " + local_ski_path + "\n")
        simulation_file.write("remote skifile path: " + arguments.ski_pattern + "\n")
        simulation_file.write("local input directory: " + str(local_input_path) + "\n") # can be None
        simulation_file.write("local output directory: " + local_output_path + "\n")
        simulation_file.write("remote input directory: " + str(arguments.input_path) + "\n") # can be None
        simulation_file.write("remote simulation directory: " + str(remote_simulation_path) + "\n")
        simulation_file.write("remote output directory: " + arguments.output_path + "\n")
        simulation_file.write("submitted at: " + time.timestamp() + "\n")

        # Close the file
        simulation_file.close()

        # Return the path to the simulation file
        return simulation_file_path

    # -----------------------------------------------------------------

    def schedule(self, arguments, name, scheduling_options, local_ski_path, remote_simulation_path):

        """
        This function ...
        :return:
        """

        # Inform the suer
        self.log.info("Scheduling simulation on the remote host")

        if scheduling_options is None:

            processors = arguments.parallel.processes * arguments.parallel.threads
            nodes, ppn = self.get_requirements(processors)

            scheduling_options = {}
            scheduling_options["nodes"] = nodes
            scheduling_options["ppn"] = ppn
            #scheduling_options["walltime"] = ...
            #scheduling_options["jobscript_path"] = ...
            scheduling_options["mail"] = False
            scheduling_options["full_node"] = True

        # We want to estimate the walltime here if it is not passed to this function
        if "walltime" not in scheduling_options:

            factor = 1.2

            # Create and run a ResourceEstimator instance
            estimator = ResourceEstimator()
            #estimator.run(local_ski_path, arguments.parallel.processes, arguments.parallel.threads)
            estimator.run(local_ski_path, 1, 1)

            # Return the estimated walltime
            #walltime = estimator.walltime * factor
            walltime = estimator.walltime_for(arguments.parallel.processes, arguments.parallel.threads) * factor

        else: walltime = scheduling_options["walltime"]

        if "nodes" not in scheduling_options: raise ValueError("The number of nodes is not defined in the scheduling options")
        if "ppn" not in scheduling_options: raise ValueError("The number of processors per node is not defined in the scheduling options")
        nodes = scheduling_options["nodes"]
        ppn = scheduling_options["ppn"]

        mail = scheduling_options["mail"] if "mail" in scheduling_options else False
        full_node = scheduling_options["full_node"] if "full_node" in scheduling_options else False

        # Create a job script next to the (local) simulation's ski file
        if "jobscript_path" not in scheduling_options:
            local_simulation_path = os.path.dirname(local_ski_path)
            local_jobscript_path = os.path.join(local_simulation_path, "job.sh")
        else: local_jobscript_path = scheduling_options["jobscript_path"]
        jobscript_name = os.path.basename(local_jobscript_path)
        jobscript = JobScript(local_jobscript_path, arguments, self.host.clusters[self.config.cluster_name], self.skirt_path, self.host.mpi_command, self.host.modules, walltime, nodes, ppn, name=name, mail=mail, full_node=full_node)

        # Copy the job script to the remote simulation directory
        remote_jobscript_path = os.path.join(remote_simulation_path, jobscript_name)
        self.upload(local_jobscript_path, remote_simulation_path)

        ## Swap clusters
        # Then, swap to the desired cluster and launch the job script
        #output = subprocess.check_output("module swap cluster/" + self._clustername + "; qsub " + self._path, shell=True, stderr=subprocess.STDOUT)

        # Submit the job script to the remote scheduling system
        #output = self.execute("qsub " + remote_jobscript_path, contains_extra_eof=True)
        output = self.execute("qsub " + remote_jobscript_path)

        # The queu number of the submitted job is used to identify this simulation
        simulation_id = int(output[0].split(".")[0])

        # Return the simulation ID
        return simulation_id

    # -----------------------------------------------------------------

    def start(self, arguments):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Starting simulation on the remote host")

        # Send the command to the remote machine using a screen session so that we can safely detach from the
        # remote shell
        command = arguments.to_command(self.skirt_path, self.host.mpi_command, self.host.scheduler, to_string=True)
        self.execute("screen -d -m " + command, output=False)

        # Generate a new simulation ID based on the ID's currently in use
        simulation_id = self._new_simulation_id()

        # Return the simulation ID
        return simulation_id

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

    def _simulation_ids_in_use(self):

        """
        This function ...
        :return:
        """

        # Check the contents of the local run directory to see which simulation id's are currently in use
        current_ids = []
        for item in os.listdir(self.local_skirt_host_run_dir):

            # Determine the full path to this item
            path = os.path.join(self.local_skirt_host_run_dir, item)

            # If this item is a directory or it is hidden, skip it
            if os.path.isdir(path) or item.startswith("."): continue

            # If the file has the 'sim' extension, get the simulation ID and add it to the list
            current_ids.append(int(item.split(".sim")[0]))

        # Return the list of currently used ID's
        return current_ids

    # -----------------------------------------------------------------

    def _new_simulation_id(self):

        """
        This function ...
        :param count:
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

    def execute(self, command, output=True, expect_eof=True, contains_extra_eof=False):

        """
        This function ...
        :param command:
        :return:
        """

        # Send the command
        self.ssh.sendline(command)

        # Retreive the output if requested
        eof = self.ssh.prompt()

        # If an extra EOF is used before the actual output line (don't ask me why but I encounter this on the HPC UGent infrastructure), do prompt() again
        if contains_extra_eof: eof = self.ssh.prompt()

        # If the command could not be sent, raise an error
        if not eof and expect_eof and not contains_extra_eof: raise RuntimeError("The command could not be sent")

        # Ignore the first and the last line (the first is the command itself, the last is always empty)
        if output:
            # Trial and error to get it right for HPC UGent login nodes; don't know what is happening
            if contains_extra_eof: return self.ssh.before.replace('\x1b[K', '').split("\r\n")[1:-1]
            else: return self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]

    # -----------------------------------------------------------------

    def remove_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.execute("rm -rf " + path, output=False)

    # -----------------------------------------------------------------

    def remove_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.execute("rm " + path, output=False)

    # -----------------------------------------------------------------

    def change_cwd(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        self.execute("cd " + path)

    # -----------------------------------------------------------------

    def directories_in_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the path to the current working directory
        working_directory = self.working_directory

        # Change the working directory to the provided path
        self.change_cwd(path)

        # List the directories in the provided path
        output = self.execute("for i in $(ls -d */); do echo ${i%%/}; done")

        # Change the working directory back to the original working directory
        self.change_cwd(working_directory)

        # Return the list of directories
        return output

    # -----------------------------------------------------------------

    def files_in_path(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the path to the current working directory
        working_directory = self.working_directory

        # Change the working directory to the provided path
        self.change_cwd(path)

        # List the files in the provided path
        output = self.execute("for f in *; do [[ -d $f ]] || echo $f; done")

        # Change the working directory back to the original working directory
        self.change_cwd(working_directory)

        # Return the list of directories
        return output

    # -----------------------------------------------------------------

    def download(self, origin, destination, timeout=30):

        """
        This function ...
        :return:
        """

        # Construct the command string
        copy_command = "scp "

        # Add the host address
        copy_command += self.host.user + "@" + self.host.name + ":"

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            if self.is_file(origin): copy_command += origin.replace(" ", "\ ") + " "

            # Check if it represents a directory
            #elif self.is_directory(origin): copy_command += origin.replace(" ", "\ ") + "/* " + "-r "
            elif self.is_directory(origin): copy_command += origin.replace(" ", "\ ") + "/* "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file on the remote host")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist remotely
            for file_path in origin:
                if not self.is_file(file_path): raise ValueError("The file " + file_path + " does not exist on the remote host")

            # Escape possible space characters
            origin = [path.replace(" ", "\ ") for path in origin]

            # Add a quotation mark character because the seperate file paths are going to be separated by spaces
            # (the command is going to be of the form scp username@ip.of.server.copyfrom:"file1.log file2.log" "~/yourpathtocopy")
            copy_command += '"'

            # Add the file paths to the command string
            copy_command += " ".join(origin)

            # Add another quotation mark to identify the end of the filepath list
            copy_command += '" '

        # Add the destination path to the command
        copy_command += destination.replace(" ", "\ ") + "/"

        self.log.debug("Copy command: " + copy_command)

        # Temporary file for output of the scp command
        temp_file_path = tempfile.mktemp()
        temp_file = open(temp_file_path, 'w')

        # Execute the command
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)
        child.logfile = temp_file
        child.expect(pexpect.EOF, timeout=None)
        child.close()

        # Close the temporary output file
        temp_file.close()

        # Open the output file again and read the contents
        temp_file = open(temp_file_path, 'r')
        stdout = temp_file.read()
        temp_file.close()

        # Raise an error if something went wrong
        if child.exitstatus != 0: raise RuntimeError(stdout)

        self.log.debug("Copy stdout: " + str(stdout))

    # -----------------------------------------------------------------

    def upload(self, origin, destination, timeout=30):

        """
        This function ...
        :return:
        """

        # Construct the command string
        copy_command = "scp "

        # If the origin is a string, we assume it represents a single file path or directory path
        if isinstance(origin, basestring):

            # Check if the origin represents a file
            if os.path.isfile(origin): copy_command += origin.replace(" ", "\ ") + " "

            # Check if it represents a directory
            elif os.path.isdir(origin): copy_command += "-r " + origin.replace(" ", "\ ") + "/ "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist locally
            for file_path in origin:
                if not os.path.isfile(file_path): raise ValueError("The file " + file_path + " does not exist")

            # Escape possible space characters
            origin = [path.replace(" ", "\ ") for path in origin]

            # Add the file paths to the command string
            copy_command += " ".join(origin) + " "

        # Invalid argument
        else: raise ValueError("The origin must be a string or a list of strings")

        # Add the host address and the destination directory
        copy_command += self.host.user + "@" + self.host.name + ":" + destination.replace(" ", "\ ") + "/"
        self.log.debug("Copy command: " + copy_command)

        # Temporary file for output of the scp command
        temp_file_path = tempfile.mktemp()
        temp_file = open(temp_file_path, 'w')

        # Execute the command
        child = pexpect.spawn(copy_command, timeout=timeout)
        if self.host.password is not None:
            child.expect(['password: '])
            child.sendline(self.host.password)
        child.logfile = temp_file
        try:
            child.expect(pexpect.EOF, timeout=None)
        except pexpect.EOF:
            pass
        child.close()

        # Close the temporary output file
        temp_file.close()

        # Open the output file again and read the contents
        temp_file = open(temp_file_path, 'r')
        stdout = temp_file.read()
        temp_file.close()

        # Raise an error if something went wrong
        if child.exitstatus != 0: raise RuntimeError(stdout)

        self.log.debug("Copy stdout: " + str(stdout))

    # -----------------------------------------------------------------

    def retreive(self):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Initialize a list to contain the simulations that have been retreived
        simulations = []

        # Loop over the different entries of the status list
        for entry in self.status:

            # The path to the simulation file
            path = entry.file_path

            # The ski file path
            remote_ski_path = entry.remote_ski_path

            # The remote output path
            remote_output_path = entry.remote_output_path

            # The simulation status
            simulation_status = entry.status

            # Skip already retreived simulations
            if simulation_status == "retreived": continue

            # Finished simulations
            elif simulation_status == "finished":

                prefix = os.path.basename(remote_ski_path).split(".")[0]

                # Properties obtained from the simulation file
                local_ski_path = None
                local_input_path = None
                local_output_path = None
                remote_input_path = None
                remote_simulation_path = None
                extract_progress = None
                extract_timeline = None
                extract_memory = None
                plot_seds = None
                plot_grids = None
                plot_progress = None
                plot_timeline = None
                plot_memory = None
                make_rgb = None
                make_wave = None
                remove_remote_input = None
                remove_remote_output = None
                retreive_types = None
                extraction_directory = None
                plotting_directory = None
                scaling_run_name = None
                scaling_file_path = None
                scaling_plot_path = None
                screen_session = None

                # Get simulation properties
                with open(path) as simulation_file:

                    # Loop over the lines in the file
                    for line in simulation_file:

                        if "local skifile path" in line: local_ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        if "local input directory" in line: local_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "local output directory" in line: local_output_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote input directory" in line: remote_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote simulation directory" in line: remote_simulation_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "extract progress" in line: extract_progress = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "extract timeline" in line: extract_timeline = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "extract memory" in line: extract_memory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "plot seds" in line: plot_seds = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "plot grids" in line: plot_grids = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "plot progress" in line: plot_progress = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "plot timeline" in line: plot_timeline = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "plot memory" in line: plot_memory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "make rgb images" in line: make_rgb = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "make wave movie" in line: make_wave = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "remove remote input" in line: remove_remote_input = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "remove remote output" in line: remove_remote_output = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip() == "True"
                        elif "retreive types" in line: retreive_types = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "extraction directory" in line: extraction_directory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "plotting directory" in line: plotting_directory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "part of scaling test run" in line: scaling_run_name = line.split("scaling test run ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "scaling data file" in line: scaling_file_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "scaling plot path" in line: scaling_plot_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "launched within screen session" in line: screen_session = line.split("screen session ")[1].replace('\n', ' ').replace('\r', '').strip()

                # If retreive file types are not defined, download the complete output directory
                if retreive_types is None or retreive_types == "None" or retreive_types == "null": # 'null' is what json makes of 'None'

                    # Download the simulation output
                    self.download(remote_output_path, local_output_path)

                # If retreive file types are defined, download these files seperately to the local filesystem
                else:

                    # Try to create a list from the string that represents the retreive types
                    try: retreive_types_list = json.loads(retreive_types)
                    except ValueError: raise ValueError("The format of the retreive types string is invalid")

                    # Create a list for the paths of the files that have to be copied to the local filesystem
                    copy_paths = []

                    # Loop over the files that are present in the remoute output directory
                    for filename in self.files_in_path(remote_output_path):

                        # Determine the full path to the output file
                        filepath = os.path.join(remote_output_path, filename)

                        # Loop over the different possible file types and add the filepath if the particular type is in the list of types to retreive
                        if filename.endswith("_ds_isrf.dat"):
                            if "isrf" in retreive_types_list: copy_paths.append(filepath)
                        elif "_ds_temp" in filename and filename.endswith(".fits"):
                            if "temp" in retreive_types_list: copy_paths.append(filepath)
                        elif filename.endswith("_sed.dat"):
                            if "sed" in retreive_types_list: copy_paths.append(filepath)
                        elif filename.endswith("_total.fits"):
                            if "image" in retreive_types_list: copy_paths.append(filepath)
                        elif filename.endswith("_ds_celltemps.dat"):
                            if "celltemp" in retreive_types_list: copy_paths.append(filepath)
                        elif "_log" in filename and filename.endswith(".txt"):
                            if "log" in retreive_types_list: copy_paths.append(filepath)
                        elif filename.endswith("_wavelengths.dat"):
                            if "wavelengths" in retreive_types_list: copy_paths.append(filepath)
                        elif "_ds_grid" in filename and filename.endswith(".dat"):
                            if "grid" in retreive_types_list: copy_paths.append(filepath)
                        elif "_ds_grho" in filename and filename.endswith(".fits"):
                            if "grho" in retreive_types_list: copy_paths.append(filepath)
                        elif "_ds_trho" in filename and filename.endswith(".fits"):
                            if "trho" in retreive_types_list: copy_paths.append(filepath)
                        elif filename.endswith("_ds_convergence.dat"):
                            if "convergence" in retreive_types_list: copy_paths.append(filepath)

                    # Download the list of files to the local output directory
                    self.download(copy_paths, local_output_path)

                # If retreival was succesful, add this information to the simulation file
                with open(path, "a") as simulation_file:
                    simulation_file.write("retreived at: " + time.timestamp() + "\n")

                # Remove the remote input, if requested
                if remove_remote_input: self.remove_directory(remote_input_path)

                # Remove the remote output, if requested
                if remove_remote_output: self.remove_directory(remote_output_path)

                # If both the input and output directories have to be removed, the remote simulation directory
                # can be removed too
                if remove_remote_input and remove_remote_output: self.remove_directory(remote_simulation_path)

                # Create a simulation object and add it to the list
                simulation = SkirtSimulation(prefix, inpath=local_input_path, outpath=local_output_path, ski_path=local_ski_path)
                simulation.extract_progress = extract_progress
                simulation.extract_timeline = extract_timeline
                simulation.extract_memory = extract_memory
                simulation.plot_seds = plot_seds
                simulation.plot_grids = plot_grids
                simulation.plot_progress = plot_progress
                simulation.plot_timeline = plot_timeline
                simulation.plot_memory = plot_memory
                simulation.make_rgb = make_rgb
                simulation.make_wave = make_wave
                simulation.extraction_path = extraction_directory
                simulation.plot_path = plotting_directory
                if scaling_run_name is not None: simulation.scaling_run_name = scaling_run_name
                if scaling_file_path is not None: simulation.scaling_file_path = scaling_file_path
                if scaling_plot_path is not None: simulation.scaling_plot_path = scaling_plot_path
                if screen_session is not None: simulation.screen_session = screen_session

                # Add the simulation to the list of retreived simulations
                simulations.append(simulation)

        # Return the list of retreived simulations
        return simulations

    # -----------------------------------------------------------------

    def install(self, private=False, key_password=None):

        """
        This function ...
        :return:
        """

        # Navigate to the home directory
        self.execute("cd ~", output=False)

        # Create the SKIRT directory
        self.execute("mkdir SKIRT", output=False)

        # In the SKIRT directory, create the necessary subdirectories
        self.execute("cd SKIRT", output=False)
        self.execute("mkdir git run release", output=False)

        # Clone the SKIRT repository
        if private:
            output = self.execute("git clone git@github.ugent.be:SKIRT/SKIRT.git git", expect_eof=False)
            self.ssh.expect(['id_rsa: '])
            self.ssh.sendline(key_password)

        else: self.execute("git clone https://github.com/SKIRT/SKIRT.git git", output=False)

        # Compile the SKIRT code
        self.execute("./makeSKIRT.sh", output=False)

        # Put SKIRT in the PATH environment variable


    # -----------------------------------------------------------------

    def update(self):

        """
        This function ...
        :return:
        """

        # Navigate to the SKIRT repository directory
        skirt_repo_dir = os.path.join(self.skirt_dir, "git")
        self.execute("cd " + skirt_repo_dir, output=False)

        # Update SKIRT
        self.execute("git pull origin master", output=False)

        # Compile the SKIRT code
        self.execute("./makeSKIRT.sh", output=False)

    # -----------------------------------------------------------------

    @property
    def system_name(self):

        """
        This function ...
        :return:
        """

        name = self.config.host_id
        if self.host.scheduler: name += "-" + self.config.cluster_name
        return name

    # -----------------------------------------------------------------

    @property
    def home_directory(self):

        """
        This function ...
        :return:
        """

        # Find out the path to the user's home directory and return it
        output = self.execute("echo $HOME")
        return output[0]

    # -----------------------------------------------------------------

    @property
    def working_directory(self):

        """
        This function ...
        :return:
        """

        # Find out the path to the current working directory and return it
        output = self.execute("echo $PWD")
        return output[0]

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

    @property
    def free_cores(self):

        """
        This function ...
        :return:
        """

        return self.cores * self.cpu_load

    # -----------------------------------------------------------------

    @property
    def free_memory(self):

        """
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        #output = self.execute("free -t | grep 'Total'")
        output = self.execute("free -t | grep buffers/cache")
        splitted = output[0].split(":")[1].split()

        # Calculate the free amount of memory in gigabytes
        free = float(splitted[1]) / 1e6

        # Return the free amount of virtual memory in gigabytes
        return free

    # -----------------------------------------------------------------

    @property
    def free_space(self):

        """
        This function ...
        :return:
        """

        # Use the 'df' command to obtain information about the free disk space
        output = self.execute("df -lh")

        total = 0.0
        used = 0.0
        free = 0.0

        for entry in output[1:]:

            if not entry.startswith("/dev/"): continue

            splitted = entry.split()

            total += float(splitted[1].split("G")[0]) if "G" in splitted[1] else float(splitted[1].split("T")[0]) * 1e3
            used += float(splitted[2].split("G")[0]) if "G" in splitted[2] else float(splitted[2].split("T")[0]) * 1e3
            free += float(splitted[3].split("G")[0]) if "G" in splitted[3] else float(splitted[3].split("T")[0]) * 1e3

        # Return the amount of free memory in gigabytes
        return free

    # -----------------------------------------------------------------

    @property
    def cores(self):

        """
        This function ...
        :return:
        """

        # If the remote host uses a scheduling system, the number of cores on the computing nodes is defined in the configuration
        if self.host.scheduler: return self.host.clusters[self.config.cluster_name].cores

        # If no scheduler is used, the computing node is the actual node we are logged in to
        else:

            # Use the 'lscpu' command to obtain the number of CPU's (hardware threads)
            output = self.execute("lscpu | grep '^CPU(s)'")
            cpus = int(float(output[0].split(":")[1]))

            # Use the 'lscpu' command to get the number of hardware threads per core
            output = self.execute("lscpu | grep '^Thread(s) per core'")
            threads_per_core = int(float(output[0].split(":")[1]))

            # Return the number of physical cores
            return cpus / threads_per_core

    # -----------------------------------------------------------------

    @property
    def cpu_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'top' command to get the current CPU load
        output = self.execute("top -b -n1 | grep 'Cpu(s)' | awk '{print $2 + $4}'")

        # Convert the output into the fraction of CPU that is used
        load = float(output[0]) / 100.0

        # Return the current CPU load
        return load

    # -----------------------------------------------------------------

    @property
    def memory_load(self):

        """
        This function ...
        :return:
        """

        # Use the 'free' command to get information about the virtual memory usage
        output = self.execute("free -t | grep 'Total'")
        splitted = output[0].split(":")[1].split()

        # Calculate the total and used amount of memory in gigabytes
        total = float(splitted[0]) / 1e6
        used = float(splitted[1]) / 1e6
        free = float(splitted[2]) / 1e6

        # Return the fraction of virtual memory that is currently in use
        return used / total

    # -----------------------------------------------------------------

    def file_or_directory(self, path):

        """
        This function ...
        :return:
        """

        print(path)

        # Launch a bash command to check whether the path exists as either a file or a directory on the remote file system
        output = self.execute("if [ -f " + path + " ]; then echo file; elif [ -d " + path + " ]; then echo directory; else echo None; fi")

        # Return the result
        if output[0] == "None": return None
        else: return output[0]

    # -----------------------------------------------------------------

    def is_directory(self, path):

        """
        This function ...
        :return:
        """

        # Launch a bash command to check whether the path exists as a directory on the remote file system
        output = self.execute("if [ -d " + path + " ]; then echo True; else echo False; fi")

        # Return the result
        return output[0] == "True"

    # -----------------------------------------------------------------

    def is_file(self, path):

        """
        This function ...
        :return:
        """

        # Launch a bash command to check whether the path exists as a regular file on the remote file system
        output = self.execute("if [ -f " + path + " ]; then echo True; else echo False; fi")

        # Return the result
        return output[0] == "True"

    # -----------------------------------------------------------------

    @property
    def host_name(self):

        """
        This function ...
        :return:
        """

        return self.config.host_id

    # -----------------------------------------------------------------

    @property
    def status(self):

        """
        This function ..
        :return:
        """

        # Initialize a list to contain the simulations
        simulations = []

        # If the remote host does not use a scheduling system
        if not self.host.scheduler:

            # Search for files in the local SKIRT run/host_id directory
            for item in os.listdir(self.local_skirt_host_run_dir):

                # If the item is not a simulation file or it is hidden, skip it
                if not item.endswith(".sim") or item.startswith("."): continue

                # Determine the full path to the simulation file
                path = os.path.join(self.local_skirt_host_run_dir, item)

                # Determine the simulation ID
                simulation_id = int(os.path.basename(path).split(".sim")[0])

                # Open the file
                name = None
                local_ski_path = None
                remote_ski_path = None
                remote_input_path = None
                remote_output_path = None
                remote_simulation_path = None
                screen_name = None

                retreived = False

                with open(path) as simulation_file:

                    # Loop over the lines in the file
                    for line in simulation_file:

                        if "simulation name" in line: name = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "local skifile path" in line: local_ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote skifile path" in line: remote_ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote input directory" in line: remote_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote output directory" in line: remote_output_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote simulation directory" in line: remote_simulation_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "launched within screen session" in line: screen_name = line.split("session ")[1].replace('\n', ' ').replace('\r', '').strip()

                        elif "retreived at" in line: retreived = True

                    if remote_ski_path is None: raise ValueError("Not a valid simulation file")
                    if remote_output_path is None: raise ValueError("Not a valid simulation file")

                # The name of the ski file (the simulation prefix)
                ski_name = os.path.basename(remote_ski_path).split(".")[0]

                # The path to the simulation log file
                remote_log_file_path = os.path.join(remote_output_path, ski_name + "_log.txt")

                if retreived: simulation_status = "retreived"
                else:
                    # Get the simulation status from the remote log file
                    simulation_status = self.status_from_log_file(remote_log_file_path, screen_name, ski_name)

                # Add the simulation properties to the list
                simulation = Map()
                simulation.id = simulation_id
                simulation.file_path = path
                simulation.name = name
                simulation.local_ski_path = local_ski_path
                simulation.remote_ski_path = remote_ski_path
                simulation.remote_input_path = remote_input_path
                simulation.remote_output_path = remote_output_path
                simulation.remote_simulation_path = remote_simulation_path
                simulation.status = simulation_status
                simulations.append(simulation)

            # Return the list of simulation properties
            return simulations

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
            for item in os.listdir(self.local_skirt_host_run_dir):

                # If the item is not a simulation file or it is hidden, skip it
                if not item.endswith(".sim") or item.startswith("."): continue

                # Determine the full path to the simulation file
                path = os.path.join(self.local_skirt_host_run_dir, item)

                # Open the file
                name = None
                local_ski_path = None
                remote_ski_path = None
                remote_input_path = None
                remote_output_path = None
                remote_simulation_path = None

                retreived = False

                with open(path) as simulation_file:

                    # Loop over the lines in the file
                    for line in simulation_file:

                        if "simulation name" in line: name = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "local skifile path" in line: local_ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote skifile path" in line: remote_ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote input directory" in line: remote_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote output directory" in line: remote_output_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote simulation directory" in line: remote_simulation_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()

                        elif "retreived at" in line: retreived = True

                    if remote_ski_path is None: raise ValueError("Not a valid simulation file")
                    if remote_output_path is None: raise ValueError("Not a valid simulation file")

                # The name of the ski file (the simulation prefix)
                ski_name = os.path.basename(remote_ski_path).split(".")[0]

                # The path to the simulation log file
                remote_log_file_path = os.path.join(remote_output_path, ski_name + "_log.txt")

                # Get the job ID from the name of the simulation file
                job_id = int(os.path.splitext(item)[0])

                # Check if the simulation has already been retreived
                if retreived: simulation_status = "retreived"

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
                simulation = Map()
                simulation.id = job_id
                simulation.file_path = path
                simulation.name = name
                simulation.local_ski_path = local_ski_path
                simulation.remote_ski_path = remote_ski_path
                simulation.remote_input_path = remote_input_path
                simulation.remote_output_path = remote_output_path
                simulation.remote_simulation_path = remote_simulation_path
                simulation.status = simulation_status
                simulations.append(simulation)

            # Return the list of simulation properties
            return simulations

    # -----------------------------------------------------------------

    def status_from_log_file(self, file_path, screen_name, simulation_prefix):

        """
        This function ...
        :param file_path:
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
        :param job_id:
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
        # we would have encountered it's job ID in the qstat output)
        else: simulation_status = "cancelled"

        # Return the string that indicates the simulation status
        return simulation_status

    # -----------------------------------------------------------------

    def read_text_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the text file into a variable
        self.execute("value='cat " + path + "'")

        # Print the variable to the console, and obtain the output
        return self.execute('echo "$($value)"')

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

        if "self-absorption" in phase: return "running " + str(phase) + ", cycle " + str(cycle) + "] " + str(progress) + "%)"
        elif "stellar emission" in phase or "dust emission" in phase: return "running (" + str(phase) + " " + str(progress) + "%)"
        else: return "running (" + str(phase) + ")"

# -----------------------------------------------------------------
