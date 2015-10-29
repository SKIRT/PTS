#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_insert Insert selected records into the SKIRT-runs database.
#
# This script populates the SKIRT-runs database with relevant galaxies selected from the public EAGLE database.
#
# The script expects exactly six command-line arguments specifying respectively:
#  - label: a label that will identify the set of inserted records in the SKIRT-runs database
#  - skitemplate: the name of the ski file template (without extension) to be used for these runs
#  - eaglesim: the short name identifying the snapshot from which to extract galaxies (e.g. Ref100)
#  - minstarmass: minimum stellar mass within 30kpc aperture (in solar mass units)
#  - maxstarmass: maximum stellar mass within 30kpc aperture (in solar mass units)
#  - fraction: random fraction of the galaxies to select
#

# -----------------------------------------------------------------

import sys
import eagle.config as config
import eagle.database
from eagle.connection import Connection
from eagle.galaxy import Snapshot

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 7: raise ValueError("This script expects exactly six command-line arguments: " \
                + "label skitemplate eaglesim minstarmass maxstarmass fraction")
label = sys.argv[1]
skitemplate = sys.argv[2]
eaglesim = sys.argv[3]
minstarmass = float(sys.argv[4])
maxstarmass = float(sys.argv[5])
fraction = float(sys.argv[6])

# get a list of the requested galaxies from the public EAGLE database
print "Querying public EAGLE database..."
con = Connection(config.public_eagle_database_username, password=config.public_eagle_database_password)
records = con.execute_query('''
    SELECT
        gal.GalaxyID as galaxyid
    FROM
        {0}_SubHalo as gal,
        {0}_Aperture as ape
    WHERE
        ape.Mass_Star > {1} and
        ape.Mass_Star <= {2} and
        gal.RandomNumber <= {3} and
        gal.SnapNum = 28 and
        gal.Spurious = 0 and
        ape.ApertureSize = 30 and
        gal.GalaxyID = ape.GalaxyID
    '''.format(config.public_eagle_database_name[eaglesim], minstarmass, maxstarmass, fraction))

# display some information for verification
print "Selected {} galaxies with stellar masses in range {:.2e} .. {:.2e}".format(len(records),minstarmass,maxstarmass)
print "SKIRT-runs will be inserted for these galaxies with:"
print "  label:", label
print "  skitemplate:", skitemplate
print "  eaglesim:", eaglesim

# ask confirmation
confirm = raw_input("--> Would you like to insert {} new SKIRT-runs? [y/n]: ".format(len(records)))
if not confirm.lower().startswith('y'): exit()

# backup the data base
eagle.database.backup()

# insert a new record into the database for each selected galaxy
db = eagle.database.Database()
with db.transaction():
    for record in records:
        db.insert(label, eaglesim, 0, record["galaxyid"], skitemplate)
db.close()

# confirm completion
print "{} new SKIRT-runs were inserted in the database".format(len(records))

# -----------------------------------------------------------------
