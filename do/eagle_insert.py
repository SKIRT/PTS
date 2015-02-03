#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_insert Insert selected records into the SKIRT-runs database.
#
# This script populates the SKIRT-runs database with relevant galaxies selected from the default EAGLE snapshot.
#
# The script expects exactly four command-line arguments specifying respectively:
#  - a label that will identify the set of inserted records in the database
#  - the name of the ski file template (without extension) to be used for these runs
#  - a minimum gas mass (in solar mass units)
#  - a maximum gas mass (in solar mass units)
#

# -----------------------------------------------------------------

import sys
from eagle.galaxy import Snapshot
import eagle.database

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 5: raise ValueError("This script expects exactly four command-line arguments")
label = sys.argv[1]
skiname = sys.argv[2]
mingasmass = float(sys.argv[3])
maxgasmass = float(sys.argv[4])

# backup and open the data base
eagle.database.backup()
db = eagle.database.Database()

# make sure the specified label is new
count = len(db.select("label=?", (label,)))
if count > 0: raise ValueError("This label is already in use for " + str(count) + " records in the database")

# open the snapshot files
print "Opening the snapshot..."
snap = Snapshot()
snap.printinfo()
print ""

# get the list of galaxies within the specified gass mass range
galaxies = snap.galaxies()
galaxies.remove_gasmass_below(mingasmass)
galaxies.remove_gasmass_above(maxgasmass)

# insert a new record into the database for each selected galaxy, without committing
for g in galaxies.galaxies:
    db.insert(label, snap.eaglesim, snap.orig_redshift,
              g.groupnumber, g.subgroupnumber, g.numstarparticles, g.numgasparticles, g.starmass, g.gasmass,
              skiname)

# show the result and ask for confirmation
records = db.select("label=?", (label,))
count = len(records)
db.show(records)
confirm = raw_input("--> Would you like to commit these " + str(count) + " new records to the database? [y/n]: ")

# commit or reject based on user reply
if confirm.lower().startswith('y'):
    db.commit()
    print str(count) + " new records were committed to the database"
else:
    print "New records were rejected"

# close the database
db.close()

# -----------------------------------------------------------------
