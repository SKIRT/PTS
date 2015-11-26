#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_schedule Schedule EAGLE SKIRT jobs from the SKIRT-runs database
#
# This script schedules EAGLE SKIRT jobs for selected records in the SKIRT-runs database. If the script is executed
# on the Cosma cluster, it creates and submits actual batch jobs to the queuing system to perform the runs.
# If the script is executed on a regular computer, its asks the user to perform the skirt-runs by hand. In both cases,
# the script updates the run-status of the selected records in the database to 'scheduled'.
#
# The script schedules a separate batch job for each command line argument. Each command line argument specifies
# one or more run-ids as a comma-seperated list of run-ids and/or runid ranges expressed with a dash (no spaces!).
# For example:
#
#\verbatim
#$ pts schedule 6 16,9 7,10-14
#\endverbatim
#
# would schedule three jobs, respectively executing run-id sequences [6], [9,16] and [7,10,11,12,13,14]
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle.scheduler import schedule

# -----------------------------------------------------------------

# loop over the command-line arguments
for arg in sys.argv[1:]:

    # schedule a job for the specified sequence
    schedule(arg)

# -----------------------------------------------------------------
