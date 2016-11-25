#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.status Show a status summary per workflow stage for records in the SKIRT-runs database.
#
# This shows the number of records in the SKIRT-runs database at each processing status for each workflow stage.
#
# The script expects a single command-line argument specifying the label for the records to be considered.
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle import database

# -----------------------------------------------------------------

def showstatus(db, label, stage, status):
    # count the relevant records
    count = len(db.select("label=? and stage=? and status=?", (label, stage, status)))
    if count>0: print label, stage, status, ":", count
 
# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv)!=2: raise ValueError("This script expects a single command-line argument: label")
label = sys.argv[1]

# open the database
db = database.Database()

# show status for each status/stage combination
for stage in ('insert', 'extract', 'simulate', 'observe', 'store'):
    for status in ('scheduled', 'running', 'failed', 'succeeded'):
        showstatus(db, label, stage, status)

# close the database
db.close()


# -----------------------------------------------------------------
