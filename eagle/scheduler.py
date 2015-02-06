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
import eagle.config as config
from eagle.database import Database
from eagle.skirtrun import SkirtRun

# -----------------------------------------------------------------

## This function schedules a SKIRT job corresponding to the SKIRT-runs database record with the specified run-id,
# regardless of its original run-status.
# If the function is executed on the Cosma cluster, it creates and submits an actual job to the queuing system
# to perform the run. If the function is executed on a regular computer, its asks the user to perform the job by hand.
# In both cases, the script updates the run-status of the specified record in the database to 'scheduled',
# and it appropriately updates the value of the 'queue' field.
#
def schedule(runid):
    # open the database
    db = Database()

    # get the record from the database
    rows = db.select("runid = ?", (runid,))
    if len(rows) != 1: raise ValueError("The specified run-id does not match a database record: " + str(runid))
    record = rows[0]

    # create the appropriate SKIRT result directories
    skirtrun = SkirtRun(runid, create=True)

    if config.queue!=None:
        # create the bash script that will be submitted as a job
        # option -n xx specifies the total number of processes (nodes_per_job * processes_per_node)
        # option -R "span[ptile=xx]" specifies the number of processes on each node (processes_per_node)
        # option -x specifies exclusive access to the nodes so we're free to spawn multiple threads
        jobscriptname = os.path.join(skirtrun.runpath(), "job.sh")
        joblogname = os.path.join(skirtrun.runpath(), "job_log.txt")
        jobscript = open(jobscriptname,'w')
        jobscript.write("#!/bin/bash -l\n")
        jobscript.write("# Batch script for running SKIRT on the Cosma cluster\n")
        jobscript.write("#BSUB -L /bin/bash\n")
        jobscript.write("#BSUB -n {}\n".format(config.nodes_per_job*config.processes_per_node))
        jobscript.write("#BSUB -R \"span[ptile={}]\"\n".format(config.processes_per_node))
        jobscript.write("#BSUB -x\n")
        jobscript.write("#BSUB -q {}\n".format(config.queue))
        jobscript.write("#BSUB -P dp004-eagle\n")
        jobscript.write("#BSUB -J SKIRT-run-{}\n".format(runid))
        jobscript.write("#BSUB -cwd {}\n".format(skirtrun.runpath()))
        jobscript.write("#BSUB -oo {}\n".format(joblogname))
        jobscript.write("#BSUB -eo {}\n".format(joblogname))
        jobscript.write("#BSUB -W {}:00\n".format(config.maximum_hours))
        jobscript.write("python -u -m do eagle_run {}\n".format(runid))
        jobscript.close()

        # update the record to indicate that a run has been scheduled and submit the job;
        # do this within a single transaction context to ensure that the scheduled job sees the updated record
        print "Submitting job for run-id " + str(runid) + " to queue " + config.queue
        with db.transaction():
            db.updatestatus((runid,), 'scheduled')
            db.updatefield((runid,), 'queue', config.queue)
            subprocess.call(("bsub",), stdin=open(jobscriptname))

    else:
        print "Please manually execute scheduled job for run-id " + str(runid)
        with db.transaction():
            db.updatestatus((runid,), 'scheduled')
            db.updatefield((runid,), 'queue', config.hostname)

    # close the database
    db.close()

## This function schedules a SKIRT job for the specified user's SKIRT-runs database records that have a runstatus
# of 'inserted'. If the username is omitted, it defaults to the current user. After the function returns,
# the records will have a runstatus of 'scheduled' and the 'queue' field will be updated appropriately.
def scheduleall(username=None):
    if username==None: username = config.username;

    # get the eligible records
    db = Database()
    records = db.select("username=? and runstatus='inserted'", (username,))
    db.close()

    # schedule a job for each of them
    for record in records:
        schedule(record['runid'])

# -----------------------------------------------------------------
