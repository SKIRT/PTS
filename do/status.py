#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.status Check the status of running SKIRT simulations on the cluster

# -----------------------------------------------------------------

# Import standard modules
import argparse
import os.path
import subprocess
from datetime import datetime
from operator import itemgetter
from collections import defaultdict
from distutils.spawn import find_executable

# Import the relevant PTS modules
from pts.log import Log
from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------

# This private helper function returns the datetime object corresponding to the time stamp in a line
def _timestamp(line):

    date, time = line.split(" ", 2)
    day, month, year = date.split('/')
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')
    return datetime(year=int(year), month=int(month), day=int(day),
                    hour=int(hour), minute=int(minute), second=int(second), microsecond=int(microsecond))

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--delete', nargs='+', type=int)

# Parse the command line arguments
args = parser.parse_args()
delete = args.delete

# Create a logger
log = Log()

# Determine the full path to the SKIRT executable
skirtpath = find_executable("skirt")

# Determine the path to the SKIRT directory
skirtdir = os.path.dirname(os.path.dirname(os.path.dirname(skirtpath)))

# Determine the path to the SKIRT/run directory
rundir = os.path.join(skirtdir, "run")

# Execute the qstat command and get the output
output = subprocess.check_output("qstat", shell=True, stderr=subprocess.STDOUT)

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
            simulation = SkirtSimulation(prefix=skifilename+".ski", outpath=outputdir)
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

        if delete is not None and str(rank) in delete:

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

# -----------------------------------------------------------------
