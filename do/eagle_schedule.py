#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_schedule Schedule EAGLE SKIRT jobs from the SKIRT-runs database
#
# This script schedules EAGLE SKIRT jobs for selected records in the SKIRT-runs database. If the script is executed
# on the Cosma cluster, it creates and submits an actual job to the queuing system to perform the run.
# If the script is executed on a regular computer, its asks the user to perform the job by hand. In both cases,
# the script updates the run-status of the selected records in the database to 'scheduled'.
#
# If this script is invoked without command line arguments, it schedules a job for each record that has a run-status
# of 'inserted' and a username equal to the current user. After the script exits, these records will have a run-status
# of 'scheduled'.
#
# Alternatively the script accepts a single command line argument specifying the run-id for a database record.
# A job will be scheduled for the specified record regardless of its original run-status, and the record's run-status
# will be updated to 'scheduled'.
#

# -----------------------------------------------------------------

import sys
import eagle.scheduler

# -----------------------------------------------------------------

# get the command-line argument
runid = int(sys.argv[1]) if len(sys.argv) > 1 else 0

if runid>0:
    # schedule a job for the specified run-id
    eagle.scheduler.schedule(runid)

else:
    # schedule a job for all eligible records for the current user
    eagle.scheduler.scheduleall()

# -----------------------------------------------------------------
