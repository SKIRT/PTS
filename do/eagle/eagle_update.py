#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_update Update a field for selected records in the EAGLE SKIRT-runs database.
#
# This script updates the value of a field for selected records in the EAGLE SKIRT-runs database.
# The script expects exactly three command-line arguments specifying respectively:
#  - the selection criteria for the records to be updated in the form of
#    an SQL expression that can be used in a WHERE clause.
#  - the name of the field to be updated
#  - the new value for the field
#
# For example:
#
#\verbatim
#pts eagle_update "runid=403" runstatus scheduled
#pts eagle_update "username='pcamps' and runid>555" skitemplate pan
#\endverbatim

# -----------------------------------------------------------------

import sys
import eagle.database

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 4: raise ValueError("This script expects exactly three command-line arguments")
dummy, selection, fieldname, newvalue  = sys.argv

# open the database
db = eagle.database.Database()

# update the records (but don't yet commit)
records = db.select(selection)
if fieldname == 'runstatus':
    db.updatestatus(records, newvalue)
else:
    db.updatefield(records, fieldname, newvalue)

# show the result and ask for confirmation
db.show(records, refetch=True)
confirm = raw_input("--> Would you like to commit this update to the database? [y/n]: ")

# commit or reject based on user reply
if confirm.lower().startswith('y'):
    db.commit()
    print "Update was committed to the database"
else:
    print "Update was rejected"

# close the database
db.close()

# -----------------------------------------------------------------
