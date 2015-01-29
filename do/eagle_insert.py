#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_insert Insert selected records into the SKIRT-runs database.
#
# This script populates the SKIRT-runs database with relevant galaxies selected from a specific EAGLE snapshot.
# It mainly serves as an example that can be adapted to the user's current needs.
#
# The script expects exactly three command-line arguments specifying respectively:
#  - a label that will identify the set of inserted records in the database
#  - a minimum number of particles (for both the stars and the gas)
#  - a maximum number of particles (for both the stars and the gas)
#
# TODO: use James' mechanisms to extract info from the EAGLE snapshot instead of the obsolete galaxy_old package.
#
# TODO: provide more meaningful selection criteria.

# -----------------------------------------------------------------

import sys
from eagle.galaxy import Snapshot, Galaxy
import eagle.database

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 4: raise ValueError("This script expects exactly three command-line arguments")
label = sys.argv[1]
minparticles = int(sys.argv[2])
maxparticles = int(sys.argv[3])

# backup and open the data base
eagle.database.backup()
db = eagle.database.Database()

# make sure the specified label is new
count = len(db.select("label=?", (label,)))
if count > 0: raise ValueError("This label is already in use for " + str(count) + " records in the database")

# open the snaphot files
print "Opening the snaphot..."
snap = Snapshot("Ref12", redshift=0)
snap.printinfo()

# get the list of galaxies
galaxies = snap.galaxies()
galaxies.remove_starmass_above(1000) # these overly large masses (>10^12 solar masses) are probably an anomaly
galaxies.remove_radius_above(2000)   # these overly large radii probably result from lying across the box border

# insert a new record into the database for each selected galaxy, without committing
for g in galaxies.galaxies:
    if g.numgasparticles>=minparticles and g.numstarparticles>=minparticles and \
       g.numgasparticles<=maxparticles and g.numstarparticles<=maxparticles:
        db.insert(label, "Ref12", 0,
              g.groupnumber, g.subgroupnumber, g.numstarparticles, g.numgasparticles, g.starmass(), g.gasmass(),
              "oli")

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
