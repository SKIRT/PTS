#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.status Check the status of running SKIRT simulations on the cluster

# -----------------------------------------------------------------

# Import standard modules
import os.path
from distutils.spawn import find_executable
import subprocess

# Import the relevant PTS modules
from pts.log import Log
from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Determine the full path to the SKIRT executable
skirtpath = find_executable("skirt")

# Determine the path to the SKIRT directory
skirtdir = os.path.dirname(os.path.dirname(os.path.dirname(skirtpath)))

# Determine the path to the SKIRT/run directory
rundir = os.path.join(skirtdir, "run")

# Execute the qstat command and get the output
output = subprocess.check_output("qstat -a", shell=True, stderr=subprocess.STDOUT)

# Create a dictionary that contains the status of the different jobs that are scheduled or running on the cluster
jobs = dict()

# Check every line in the output
for line in output.splitlines():

    # If this line mentions a job
    if "delca" in line:

        # Get the job ID
        jobid = int(line.split(".")[0])

        # Get the status (Q=queued, R=running)
        jobstatus = line.split(" ")[-3]

        # Add the status of this job to the dictionary
        jobs[jobid] = jobstatus

# Search for files in the SKIRT run directory
for itemname in os.listdir(rundir):

    # Define the full path to this item
    itempath = os.path.join(rundir, itemname)

    # Check whether this item is a txt file and it is not hidden
    if os.path.isfile(itempath) and itemname.endswith(".txt") and not itemname.startswith("."):

        # Open the file
        with open(itempath) as jobfile:

            # Get the lines
            lines = jobfile.readline()

            # Get the name of the ski file
            skifilename = lines[0].split("ski file name: ")[1]

            # Get the path of the ski file
            skifiledir = lines[1].split("ski file directory: ")[1]

            # Get the output directory of the simulation
            outputdir = lines[3].split("output directory: ")[1]

        # Check if the job is in the list of queued or running jobs
        if itemname in jobs:

            # Check the status of this job
            status = jobs[itemname]

            if status == 'Q':

                # This simulation is still queued
                log.info(skifilename + " in " + skifiledir + ": queued")

            elif status == 'R':

                # This simulation is currently running
                log.warning(skifilename + " in " + skifiledir + ": running")

            else:

                # This simulation has an unknown status
                log.failure(skifilename + " in " + skifiledir + ": unknown status")

        else:

            # Check the status from the simulation's log file (finished or crashed)
            simulation = SkirtSimulation(prefix=skifilename+".ski", outpath=outputdir)
            status = simulation.status()

            if status == "Finished":

                # This simulation has sucessfullly finished
                log.success(skifilename + " in " + skifiledir + ": finished")

            elif status == "Crashed":

                # This simulation has crashed with an error
                log.failure(skifilename + " in " + skifiledir + ": crashed")

            else:

                # This simulation has an unkown status
                log.failure(skifilename + " in " + skifiledir + ": unknown status")

# -----------------------------------------------------------------
