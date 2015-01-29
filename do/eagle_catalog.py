#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_insert Insert selected records into the SKIRT-runs database.
#
# This script cretes the cross-reference and catalog files specified in the eagle.galaxy module
#
# The script expects exactly one command-line argument specifying:
#  - The key for the relevant simulation listed in the eagle.config file
#
# TODO: provide more meaningful selection criteria.

# -----------------------------------------------------------------

import sys
from eagle.galaxy import Snapshot, Galaxy
import eagle.database

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 2: raise ValueError("This script expects exactly one command-line argument")
snapkey = sys.argv[1]

# warn the user that these things take time

while True:
    proceed = raw_input('Are you sure you want to rebuild the catalog for '+ snapkey +'? (y/n)\t')

    if proceed == 'n':
        sys.exit
    elif proceed == 'y':
        break

    print 'Invalid input. Please enter y or n'

# open the snaphot files
print "Opening the snaphot..."
snap = Snapshot(snapkey, redshift=0)
snap.printinfo()

# make crossreference and catalog files
snap.exportcrossreference()
snap.exportcatalog()



# -----------------------------------------------------------------
