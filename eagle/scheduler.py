#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.scheduler Functions to schedule EAGLE SKIRT jobs from the SKIRT-run database
#
# This module offers some functions to schedule EAGLE SKIRT jobs according to the specifications defined in
# the SKIRT-run database.
#
# When the module is used on the Durham Cosma cluster, it creates and submits jobs to the LSF workload management
# platform in accordance with the guidelines described in
# <a href="http://icc.dur.ac.uk/index.php?content=Computing/Batch#section1.1">Job Scheduling on
# the Durham COSMA machines</a>.
#

# -----------------------------------------------------------------

import os.path
import subprocess

from . import config as config
from .database import Database
from .skirtrun import SkirtRun

# -----------------------------------------------------------------

## This function schedules a single SKIRT job to perform the SKIRT-runs corresponding to database records with the
# specified sequence of run-ids, regardless of the original run-status of these records.
# If the function is executed on the Cosma cluster, it creates and submits an actual job to the queuing system
# to perform the run. If the function is executed on a regular computer, its asks the user to perform the runs by hand.
# In both cases, the script updates the run-status of the specified records in the database to 'scheduled',
# and it appropriately updates the value of the 'queue' fields.
#
# The \em runidspec argument can be an integer specifying a single run-id, or a string containing a
# comma-seperated list of run-ids and/or runid ranges expressed as two run-ids with a dash in between.
# Valid arguments include, for example, 6, "6", "6,9", "6-9" and "6-9,17".
#
def schedule(runidspec):
    # deal with the special case of a single integer
    runidspec = str(runidspec)

    # build the sequence of run-ids
    runids = []
    for segment in runidspec.split(","):
        if "-" in segment:
            first,last = map(int,segment.split("-"))
            runids += [ id for id in range(first,last+1) ]
        else:
            if segment!="": runids += [ int(segment) ]

    # open the database
    db = Database()

    # get the records corresponding to the run-id sequence from the database
    records = db.select("runid in (" + ",".join(map(str,runids)) + ")")
    if len(records) != len(runids): raise ValueError("Some of the specified run-ids do not match a database record")
    if len(records) < 1: raise ValueError("A job must have at least one run-id")

    # inform the user and get confirmation for long sequences
    if len(runids)>12:
        print "Scheduling a job for runids", runids
        proceed = raw_input("--> This job has {} runids. Do you want to proceed ? [y/n] ".format(len(runids)))
        if not proceed.lower().startswith("y"): return

    # create the appropriate SKIRT result directories for all runs in the sequence
    for runid in runids:
        skirtrun = SkirtRun(runid, create=True)

    # create and submit a batch job
    if config.queue!=None:
        # create the bash script that will be submitted as a job
        # option -n xx specifies the total number of processes (nodes_per_job * processes_per_node)
        # option -R "span[ptile=xx]" specifies the number of processes on each node (processes_per_node)
        # option -x specifies exclusive access to the nodes so we're free to spawn multiple threads
        firstrun = SkirtRun(runids[0])
        jobscriptname = os.path.join(firstrun.runpath(), "job.sh")
        joblogname = os.path.join(firstrun.runpath(), "job_log.txt")
        jobscript = open(jobscriptname,'w')
        jobscript.write("#!/bin/bash -l\n")
        jobscript.write("# Batch script for running SKIRT on the Cosma cluster\n")
        jobscript.write("#BSUB -L /bin/bash\n")
        jobscript.write("#BSUB -n {}\n".format(config.nodes_per_job*config.processes_per_node))
        jobscript.write("#BSUB -R \"span[ptile={}]\"\n".format(config.processes_per_node))
        jobscript.write("#BSUB -x\n")
        jobscript.write("#BSUB -q {}\n".format(config.queue))
        jobscript.write("#BSUB -P dp004-eagle\n")
        jobscript.write("#BSUB -J \"SKIRT {}\"\n".format(runidspec))
        jobscript.write("#BSUB -cwd {}\n".format(firstrun.runpath()))
        jobscript.write("#BSUB -oo {}\n".format(joblogname))
        jobscript.write("#BSUB -eo {}\n".format(joblogname))
        jobscript.write("#BSUB -W {}:00\n".format(config.maximum_hours))
        for runid in runids:
            skirtrun = SkirtRun(runid)
            jobscript.write("cd {}\n".format(firstrun.runpath()))
            jobscript.write("python -u -m do eagle_run {}\n".format(runid))
        jobscript.write("echo Done with SKIRT-runs {}\n".format(runidspec))
        jobscript.close()

        # update the records to indicate that a run has been scheduled and submit the job;
        # do this within a single transaction context to ensure that the scheduled job sees the updated records
        print "Submitting job to queue", config.queue, "for run-ids", runids
        with db.transaction():
            db.updatestatus(runids, 'scheduled')
            db.updatefield(runids, 'queue', config.queue)
            subprocess.call(("bsub",), stdin=open(jobscriptname))

    # or ask the user to perform the runs
    else:
        print "Please manually perform scheduled run-ids", runids
        with db.transaction():
            db.updatestatus(runids, 'scheduled')
            db.updatefield(runids, 'queue', config.hostname)

    # close the database
    db.close()

# -----------------------------------------------------------------
