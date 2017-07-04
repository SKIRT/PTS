#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.advance Advance the workflow stage of records in the SKIRT-runs database.
#
# This script advances records in the SKIRT-runs database that succeeded the current workflow stage
# to the next workflow stage, after requesting user confirmation.
#
# The script expects a single command-line argument specifying the label for the records to be considered.
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle import database

# -----------------------------------------------------------------

def advance(db, label, stagefrom, stageto):
    with db.transaction():
        # get the eligible records
        records = db.select("label=? and stage=? and status='succeeded'", (label, stagefrom))
        if len(records) > 0:

            # show some info and ask for confirmation
            print "There are {} successful SKIRT-runs with label {} at stage {}.".format(len(records), label, stagefrom)
            confirm = raw_input("--> Would you like to advance these to stage {}? [y/n]: ".format(stageto))

            # advance stage if user confirmed
            if confirm.lower().startswith('y'):
                db.updatestage(records, stageto)
                print "Advancing {} SKIRT-runs from stage {} to {}".format(len(records), stagefrom, stageto)
            else:
                print "Update was rejected"

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv)!=2: raise ValueError("This script expects a single command-line argument: label")
label = sys.argv[1]

# open the database
db = database.Database()

# advance each stage
advance(db, label, "insert", "extract")
advance(db, label, "extract", "simulate")
advance(db, label, "simulate", "observe")
advance(db, label, "observe", "store")

# close the database
db.close()


# -----------------------------------------------------------------
