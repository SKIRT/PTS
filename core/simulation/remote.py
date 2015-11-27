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
import pxssh
import subprocess

from operator import itemgetter
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics import Configurable
from .jobscript import JobScript

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
        self.ssh.sendline("which skirt")
        if self.ssh.prompt():

            # Determine the path to the SKIRT executable
            self.skirt_path = self.ssh.before.split("\n")[1].split("\r")[0]

            # Determine the path to the SKIRT directory
            self.skirt_dir = self.skirt_path.split("/release")[0]

            # Determine the path to the SKIRT run directory
            self.skirt_run_dir = os.path.join(self.skirt_dir, "run")

    # -----------------------------------------------------------------

    def login(self):

        """
        This function ...
        :return:
        """

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

        # Disconnect
        if self.connected: self.ssh.logout()

    # -----------------------------------------------------------------

    def run(self, parameters):

        """
        This function ...
        :return:
        """

        # Raise an error if a connection to the remote has not been made
        if not self.connected: raise RuntimeError("Not connected yet")

        # If the remote has a scheduling system for launching jobs
        if self.config.scheduler:

            # Create a jobscript
            jobscript_path = None
            walltime = None
            jobscript = JobScript(path, parameters, walltime)

            # Copy the job script and ski file to the remote system
            copy_command = ["scp"]
            copy_command += [jobscript_path, skifile_path]
            copy_command += destination
            process = subprocess.Popen(copy_command)

            # Launch the job script on the remote machine
            #self.ssh.sendline("qsub " + path)

        # No scheduling system
        else:

            # Send the command to the remote machine using a screen session so that we can safely detach from the
            # remote shell
            command = parameters.to_command(self.skirt_path, self.config.mpi_command, self.config.scheduler)
            self.ssh.sendline("screen -d -m " + command)

    # -----------------------------------------------------------------

    def retreive(self):

        """
        This function ...
        :return:
        """

        # Copy files from remote to local
        #scp your_username@remote.edu:/some/remote/directory/\{a,b,c\} ./

        # Copy directory from remote to local
        #scp -r user@your.server.example.com:/path/to/foo /home/user/Desktop/

    # -----------------------------------------------------------------

    def status(self):

        """
        This function ..
        :return:
        """

        # If the remote has a scheduling system for launching jobs
        if self.config.scheduler:

            # Obtain job status information through the 'qstat' command
            self.ssh.sendline("qstat")
            if self.ssh.prompt():

            # Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
            qstat = dict()

            # Check every line in the output
            for line in output.splitlines():

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
            for itemname in os.listdir(rundir):

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