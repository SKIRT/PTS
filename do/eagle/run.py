#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.run Run a SKIRT simulation on an EAGLE galaxy according to the SKIRT-runs database.
#
# This script runs a SKIRT simulation on an EAGLE galaxy, according to the specifications defined in a particular
# SKIRT-runs database record. The script expects a single command-line argument specifying the run-id of the
# database record to be processed. The run-status of the record must be 'scheduled'; if not the script fails.
# When the script exits, the runstatus of the database record is updated to either 'completed' or 'failed'.
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle import runner

# -----------------------------------------------------------------

print "Starting eagle_run..."

# get the command-line argument specifying the run-id
if len(sys.argv) != 2: raise ValueError("This script expects a single command-line argument specifying a run-id")
runid = int(sys.argv[1])

# execute for the specified run-id
runner.run(runid)

print "Finished eagle_run."

# -----------------------------------------------------------------
