#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module can be used to launch SKIRT/FitSKIRT simulations remotely
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import re
import pxssh
import pexpect
import tempfile

from operator import itemgetter
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics import Configurable
from .jobscript import JobScript
from ..tools import time, inspection
from .simulation import SkirtSimulation

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
        super(SkirtRemote, self).__init__(config)

        ## Attributes

        # Create the SSH interface
        self.ssh = pxssh.pxssh()

        # Set the connected flag to False initially
        self.connected = False

        # Set the SKIRT paths to None initially
        self.skirt_path = None
        self.skirt_dir = None
        self.skirt_run_dir = None
        self.local_skirt_run_dir = None

        # Generate a regular expression object to be used on the remote console output
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkirtRemote, self).setup()

        # Make the connection
        self.login()

        # Obtain some information about the SKIRT installation on the remote machine
        output = self.execute("which skirt")

        # Determine the path to the SKIRT executable
        self.skirt_path = output[0]  # only one line is expected

        # Determine the path to the SKIRT directory
        self.skirt_dir = self.skirt_path.split("/release")[0]

        # Determine the path to the SKIRT run directory
        self.skirt_run_dir = os.path.join(self.skirt_dir, "run")

        # Determine the path to the local SKIRT run directory
        self.local_skirt_host_run_dir = os.path.join(inspection.skirt_run_dir, self.config.host_id)

        # Create the local SKIRT run directory for this host if it doesn't already exist
        if not os.path.isdir(self.local_skirt_host_run_dir): os.makedirs(self.local_skirt_host_run_dir)

    # -----------------------------------------------------------------

    def login(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Logging in to the remote SKIRT environment")

        # Connect to the remote host
        self.connected = self.ssh.login(self.config.host, self.config.user, self.config.password)

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

    def run(self, parameters):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected to the remote")

        # Create a unique name for the simulation directory
        skifile_name = os.path.basename(parameters.ski_pattern).split(".ski")[0]
        remote_simulation_name = time.unique_name(skifile_name)

        # Determine the full path of the simulation directory on the remote system
        remote_simulation_path = os.path.join(self.skirt_run_dir, remote_simulation_name)

        # Create the remote simulation directory
        self.execute("mkdir " + remote_simulation_path, output=False)

        # Determine the full paths to the input and output directories on the remote system
        remote_input_path = os.path.join(remote_simulation_path, "in")
        remote_output_path = os.path.join(remote_simulation_path, "out")

        # Change the parameters to accomodate for the fact that we are running remotely
        # but store the paths to the local output directory because we want to copy the
        # results later
        local_input_path = parameters.input_path
        local_output_path = parameters.output_path

        if local_input_path is None: remote_input_path = None

        parameters.input_path = remote_input_path
        parameters.output_path = remote_output_path

        # Create the remote input directory if necessary
        if remote_input_path is not None: self.execute("mkdir " + remote_output_path, output=False)

        # Create the remote output directory
        self.execute("mkdir " + remote_output_path, output=False)

        local_ski_path = parameters.ski_pattern
        ski_name = os.path.basename(local_ski_path)
        remote_ski_path = os.path.join(remote_simulation_path, ski_name)
        parameters.ski_pattern = remote_ski_path

        # Copy the input directory and the ski file to the remote host
        self.upload(local_ski_path, remote_simulation_path)
        if local_input_path is not None: self.upload(local_input_path, remote_input_path)

        # If the remote has a scheduling system for launching jobs
        simulation_id = None
        if self.config.scheduler:

            # Inform the suer
            self.log.info("Scheduling simulation on the remote host")

            # Create a jobscript
            #jobscript_path = None
            #walltime = None
            #jobscript = JobScript(path, parameters, walltime)

            # Copy the job script and ski file to the remote system
            #copy_command = ["scp"]
            #copy_command += [jobscript_path, skifile_path]
            #copy_command += destination
            #process = subprocess.Popen(copy_command)

            # Launch the job script on the remote machine
            #self.ssh.sendline("qsub " + path)

            # Read the job ID from the qsub output
            #job_id = int(output.split(".")[0])
            job_id = None

            simulation_id = job_id

        # No scheduling system
        else:

            # Inform the user
            self.log.info("Starting simulation on the remote host")

            # Send the command to the remote machine using a screen session so that we can safely detach from the
            # remote shell
            command = parameters.to_command(self.skirt_path, self.config.mpi_command, self.config.scheduler)
            self.execute("screen -d -m " + command, output=False)

            # Check the contents of the local run directory to see which simulation id's are currently in use
            current_ids = []
            for item in os.listdir(self.local_skirt_host_run_dir):

                # Determine the full path to this item
                path = os.path.join(self.local_skirt_host_run_dir, item)

                # If this item is a directory or it is hidden, skip it
                if os.path.isdir(path) or item.startswith("."): continue

                # If the file has the 'sim' extension, get the simulation ID and add it to the list
                current_ids.append(int(item.split(".sim")[0]))

            # Sort the current simulation ID's and find the lowest 'missing' integer number
            if len(current_ids) > 0:
                current_ids = sorted(current_ids)
                simulation_id = max(current_ids)
                for index in range(max(current_ids)):
                    if current_ids[index] != index:
                        simulation_id = index
                        break

            if simulation_id is None: simulation_id = 0

        # Create the simulation file
        simulation_file_path = os.path.join(self.local_skirt_host_run_dir, str(simulation_id) + ".sim")
        simulation_file = open(simulation_file_path, 'w')

        # Add the simulation information
        simulation_file.write("skifile path: " + parameters.ski_pattern + "\n")
        simulation_file.write("local input directory: " + str(local_input_path) + "\n") # can be None
        simulation_file.write("local output directory: " + local_output_path + "\n")
        simulation_file.write("remote input directory: " + str(remote_input_path) + "\n") # can be None
        simulation_file.write("remote output directory: " + remote_output_path + "\n")
        simulation_file.write("submitted at: " + time.timestamp() + "\n")

        # Close the file
        simulation_file.close()

        # Return the path to the simulation file
        return simulation_file_path

    # -----------------------------------------------------------------

    def execute(self, command, output=True):

        """
        This function ...
        :param command:
        :return:
        """

        # Send the command
        success = self.ssh.sendline(command)

        # If the command could not be sent, raise an error
        if not success: raise RuntimeError("The command could not be sent")

        # Retreive the output if requested
        self.ssh.prompt()

        # Ignore the first and the last line (the first is the command itself, the last is always empty)
        if output: return self.ansi_escape.sub('', self.ssh.before).replace('\x1b[K', '').split("\r\n")[1:-1]

    # -----------------------------------------------------------------

    def download(self, origin, destination, timeout=30):

        """
        This function ...
        :return:
        """

        # Download only certain files from remote to local
        #scp your_username@remote.edu:/some/remote/directory/\{a,b,c\} ./

        # Construct the command string
        copy_command = "scp " + self.config.user + "@" + self.config.host + ":" + origin
        copy_command += destination
        self.log.info("Copy command: " + copy_command)

        # Temporary file for output of the scp command
        temp_file_path = tempfile.mktemp()
        temp_file = open(temp_file_path, 'w')

        # Execute the command
        child = pexpect.spawn(copy_command, timeout=timeout)
        child.expect(['password: '])
        child.sendline(self.config.password)
        child.logfile = temp_file
        child.expect(pexpect.EOF)
        child.close()

        # Close the temporary output file
        temp_file.close()

        # Open the output file again and read the contents
        temp_file = open(temp_file_path, 'r')
        stdout = temp_file.read()
        temp_file.close()

        # Raise an error if something went wrong
        if child.exitstatus != 0: raise RuntimeError(stdout)

        self.log.info("Copy stdout: " + str(stdout))

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
            if os.path.isfile(origin): copy_command += origin + " "

            # Check if it represents a directory
            elif os.path.isdir(origin): copy_command += "-r " + origin + " "

            # The origin does not exist
            else: raise ValueError("The specified path " + origin + " does not represent an existing directory or file")

        # If the origin is a list, we assume it contains multiple file paths
        elif isinstance(origin, list):

            # Check whether the files exist locally
            for file_path in origin:
                if not os.path.isfile(origin): raise ValueError("The file " + file_path + " does not exist")

            # Add the file paths to the command string
            copy_command += " ".join(origin)

        # Invalid argument
        else: raise ValueError("The origin must be a string or a list of strings")

        # Add the host address and the destination directory
        copy_command += self.config.user + "@" + self.config.host + ":" + destination + "/"
        self.log.info("Copy command: " + copy_command)

        # Temporary file for output of the scp command
        temp_file_path = tempfile.mktemp()
        temp_file = open(temp_file_path, 'w')

        # Execute the command
        child = pexpect.spawn(copy_command, timeout=timeout)
        child.expect(['password: '])
        child.sendline(self.config.password)
        child.logfile = temp_file
        child.expect(pexpect.EOF)
        child.close()

        # Close the temporary output file
        temp_file.close()

        # Open the output file again and read the contents
        temp_file = open(temp_file_path, 'r')
        stdout = temp_file.read()
        temp_file.close()

        # Raise an error if something went wrong
        if child.exitstatus != 0: raise RuntimeError(stdout)

        self.log.info("Copy stdout: " + str(stdout))

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
            path = entry[0]

            # The ski file path
            ski_path = entry[1]

            # The remote output path
            remote_output_path = entry[2]

            # The simulation status
            simulation_status = entry[3]

            if simulation_status == "finished":

                prefix = os.path.basename(ski_path).split(".")[0]

                # Properties obtained from the simulation file
                local_input_path = None
                local_output_path = None
                remote_input_path = None
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
                extraction_directory = None
                plotting_directory = None

                # Get simulation properties
                with open(path) as simulation_file:

                    # Loop over the lines in the file
                    for line in simulation_file:

                        if "local input directory" in line: local_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "local output directory" in line: local_output_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "remote input directory" in line: remote_input_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        #elif "submitted at" in line: submit_time = time.parse(line.split(": ")[1])
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
                        elif "extraction directory" in line: extraction_directory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                        elif "plotting directory" in line: plotting_directory = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()

                # Download the simulation output
                self.download(remote_output_path, local_output_path)

                # Remove the remote input, if requested
                if remove_remote_input:
                    remove_command = "rm -rf " + remote_input_path
                    self.execute(remove_command)

                # Remove the remote output, if requested
                if remove_remote_output:
                    remove_command = "rm -rf " + remote_output_path
                    self.execute(remove_command)

                # Create a simulation object and add it to the list
                simulation = SkirtSimulation(prefix, inpath=local_input_path, outpath=local_output_path)
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
                simulations.append(simulation)

        # Return the list of retreived simulations
        return simulations

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

        # Use the 'lscpu' command to obtain the number of CPU's (hardware threads)
        output = self.execute("lscpu | grep '^CPU(s)'")
        cpus = float(output[0].split(":")[1])

        # Use the 'lscpu' command to get the number of hardware threads per core
        output = self.execute("lscpu | grep '^Thread(s) per core'")
        threads_per_core = float(output[0].split(":")[1])

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

    @property
    def status(self):

        """
        This function ..
        :return:
        """

        # Initialize a list to contain the
        simulations = []

        # Search for files in the local SKIRT run/host_id directory
        for item in os.listdir(self.local_skirt_host_run_dir):

            # If the item is not a simulation file or it is hidden, skip it
            if not item.endswith(".sim") or item.startswith("."): continue

            # Determine the full path to the simulation file
            path = os.path.join(self.local_skirt_host_run_dir, item)

            # Open the file
            ski_path = None
            remote_output_path = None
            with open(path) as simulation_file:

                # Loop over the lines in the file
                for line in simulation_file:

                    if "skifile path" in line: ski_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()
                    elif "remote output directory" in line: remote_output_path = line.split(": ")[1].replace('\n', ' ').replace('\r', '').strip()

                if ski_path is None: raise ValueError("Not a valid simulation file")
                if remote_output_path is None: raise ValueError("Not a valid simulation file")

            # The name of the ski file (the simulation prefix)
            ski_name = os.path.basename(ski_path).split(".")[0]

            # The path to the simulation log file
            remote_log_file_path = os.path.join(remote_output_path, ski_name + "_log.txt")

            # Check whether the log file exists remotely
            output = self.execute("if [ -e " + remote_log_file_path + " ]; then echo True; else echo False; fi")

            # If the log file does not exist, the simulation has not started yet
            simulation_status = None
            if output[0] == "False": simulation_status = "not started"
            else:

                # Get the last two lines of the remote log file
                output = self.execute("tail -2 " + remote_log_file_path)

                print(output)

                # Get the last line of the actual simulation
                if " Available memory: " in output[1]: last = output[0]
                else: last = output[1]

                # Interpret the content of the last line
                if " Finished simulation " + ski_name in last: simulation_status = "finished"
                elif " *** Error: " in last: simulation_status = "crashed"
                else: simulation_status = "running"

            # Add the simulation properties to the list
            simulation = (path, ski_path, remote_output_path, simulation_status)
            simulations.append(simulation)

        # Return the list of simulation properties
        return simulations

        # If the remote has a scheduling system for launching jobs
        if self.config.scheduler:

            # Obtain job status information through the 'qstat' command
            output = self.execute("qstat")

            # Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
            qstat = dict()

            # Check every line in the output
            for line in output:

                # If this line mentions a job
                if "master15" in line:

                    # Get the job ID
                    jobid = int(line.split(".")[0])

                    parsedline = line.split(" ")

                    # Get the status (Q=queued, R=running)
                    if "short" in parsedline: position = parsedline.index("short")
                    if "long" in parsedline: position = parsedline.index("long")
                    jobstatus = parsedline[position-1]

                    # Add the status of this job to the dictionary
                    qstat[jobid] = jobstatus

            # Create a list that contains the status of all the jobs (scheduled, running or finished)
            jobs = []

            # Search for files in the SKIRT run directory
            for itemname in os.listdir(self.local_skirt_host_run_dir):

                # Define the full path to this item
                itempath = os.path.join(rundir, itemname)

                # Check whether this item is a txt file and it is not hidden
                if os.path.isfile(itempath) and itemname.endswith(".txt") and not itemname.startswith("."):

                    # Open the file
                    with open(itempath) as jobfile:

                        # Get the lines
                        lines = jobfile.read().splitlines()

                        # Get the name of the ski file
                        skifilename = lines[0].split("ski file name: ")[1]

                        # Get the path of the ski file
                        skifiledir = lines[1].split("ski file directory: ")[1]

                        # Get the output directory of the simulation
                        simulationoutputdir = lines[3].split("output directory: ")[1]
                        simulationid = os.path.basename(simulationoutputdir)

                        # Get the name of the scaling test and the name of the specific run of that test
                        runoutputdir = os.path.dirname(simulationoutputdir)
                        run = os.path.basename(runoutputdir)
                        scalingoutdir = os.path.dirname(runoutputdir)
                        scalingtest = os.path.basename(scalingoutdir)

                        # Get the time of submitting
                        timestamp = lines[4].split("submitted at: ")[1]
                        time = _timestamp(timestamp)

                    # Get the job ID from the file name
                    jobid = int(os.path.splitext(itemname)[0])

                    # Check if the job ID is in the list of queued or running jobs
                    if jobid in qstat:

                        # Check the status of this job
                        stat = qstat[jobid]

                        # This simulation is still queued
                        if stat == 'Q': status = "queued"

                        # This simulation is currently running
                        elif stat == 'R': status = "running"

                        # This simulation has an unknown status
                        else: status = "unknown"

                    else:

                        # Check the status from the simulation's log file (finished or crashed)
                        simulation = SkirtSimulation(prefix=skifilename+".ski", outpath=simulationoutputdir)
                        stat = simulation.status()

                        # This simulation has sucessfullly finished
                        if stat == "Finished": status = "finished"

                        # This simulation has crashed with an error
                        elif stat == "Crashed": status = "crashed"

                        # This simulation has an unkown status
                        else: status = "unknown"

                    # Append the properties of this job to the list
                    jobs.append([skifilename, scalingtest, run, simulationid, time, status, itempath])

            # Sort the total list of jobs based on the time (the oldest has index 0)
            jobs.sort(key=itemgetter(4))

            # Loop over all the jobs and put them in a dictionary based on the scaling test
            scalingtests = defaultdict(list)

            # Iterate over the jobs and get the index (from the ordering)
            for rank, job in enumerate(jobs):

                test = job[1]
                run = job[2]
                simulation = job[3]
                status = job[5]
                jobfilepath = job[6]

                # Add the information to the dictionary
                scalingtests[test].append([run, simulation, status, rank, jobfilepath])

            # Finally print out the information
            for scalingtest, simulations in scalingtests.items():

                log.info(scalingtest + ":")

                for simulation in simulations:

                    run = simulation[0]
                    sim = simulation[1]
                    status = simulation[2]
                    rank = simulation[3]
                    jobfilepath = simulation[4]

                    if delete is not None and rank in delete:

                        tag = "    [X] "
                        deletethisfile = True

                    else:

                        tag = "    [" + str(rank) + "] "
                        deletethisfile = False

                    if status == "queued":

                        # TODO: figure out what to do when queued jobs should be deleted

                        log.info(tag + run + " > " + sim + ": " + status)

                    elif status == "running":

                        # TODO: figure out what to do when running jobs should be deleted

                        log.warning(tag + run + " > " + sim + ": " + status)

                    elif status == "finished":

                        # Delete the job file if requested
                        if deletethisfile: os.remove(jobfilepath)

                        log.success(tag + run + " > " + sim + ": " + status)

                    elif status == "crashed":

                        # Delete the job file if requested
                        if deletethisfile: os.remove(jobfilepath)

                        log.failure(tag + run + " > " + sim + ": " + status)

                    elif status == "unknown":

                        # TODO: figure out what to do when jobs with unknown status should be deleted

                        log.failure(tag + run + " > " + sim + ": " + status)

        # No scheduling system
        else: pass

# -----------------------------------------------------------------