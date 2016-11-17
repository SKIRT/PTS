#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.insert Insert selected records into the SKIRT-runs database.
#
# This script populates the SKIRT-runs database with galaxies selected from the public EAGLE database.
#
# The script expects 5 to 7 command-line arguments specifying respectively:
#  - label: a label that will identify the set of inserted records in the SKIRT-runs database
#  - eaglesim: the name identifying the EAGLE simulation from which to extract galaxies
#  - snaptag: the tag (0-28) identifying the snapshot from which to extract galaxies
#  - minstarmass: minimum stellar mass within 30kpc aperture, in log10 of solar mass units
#  - skitemplate: the name of the ski file template (without extension) to be used for these runs
#  - [numpp]: the number of photon packages to launch -- optional
#  - [deltamax]: the maximum mass fraction in the dust grid -- optional
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle import config
from pts.eagle import database
from pts.eagle.connection import Connection

# -----------------------------------------------------------------

# get the command-line arguments
if not (6 <= len(sys.argv) <= 8): raise ValueError("This script expects 5 to 7 command-line arguments: " \
                                    + "label eaglesim snaptag minstarmass skitemplate [numpp] [deltamax]")
label = sys.argv[1]
eaglesim = sys.argv[2]
snaptag = int(sys.argv[3])
minstarmass = float(sys.argv[4])
skitemplate = sys.argv[5]
numpp = float(sys.argv[6]) if len(sys.argv) > 6 else 5e5
deltamax = float(sys.argv[7]) if len(sys.argv) > 7 else 3e-6

# get a list of the requested galaxies from the public EAGLE database
print "Querying public EAGLE database..."
con = Connection(config.public_eagle_database_username, password=config.public_eagle_database_password)
records = con.execute_query('''
    SELECT
        gal.GalaxyID as galaxyid,
        gal.groupnumber as groupnr,
        gal.subgroupnumber as subgroupnr,
        ape.mass_star as starmass,
        gal.centreofpotential_x as copx,
        gal.centreofpotential_y as copy,
        gal.centreofpotential_z as copz
    FROM
        {0}_SubHalo as gal,
        {0}_Aperture as ape
    WHERE
        gal.SnapNum = {1} and
        gal.Spurious = 0 and
        gal.GalaxyID = ape.GalaxyID and
        ape.ApertureSize = 30 and
        ape.Mass_Star > {2}
    '''.format(eaglesim, snaptag, 10**minstarmass))

# display some information for verification
print "Selected {} galaxies from {}:{} with stellar mass above 10^{}" \
            .format(len(records), eaglesim, snaptag, minstarmass)
print "SKIRT-runs will be inserted with label {} and ski template {}" \
            .format(label, skitemplate)

# ask confirmation
confirm = raw_input("--> Would you like to insert {} new SKIRT-runs? [y/n]: ".format(len(records)))
if not confirm.lower().startswith('y'): exit()

# backup the data base
database.backup()

# insert a new record into the database for each selected galaxy
db = database.Database()
firstid = db.maxrunid()+1
with db.transaction():
    for record in records:
        db.insert(label, eaglesim, snaptag, int(record["galaxyid"]), int(record["groupnr"]), int(record["subgroupnr"]),
                  float(record["starmass"]), float(record["copx"]), float(record["copy"]), float(record["copz"]),
                  skitemplate, numpp, deltamax)
lastid = db.maxrunid()
db.close()

# confirm completion
print "{} new SKIRT-runs were inserted with run-ids {}-{}".format(len(records), firstid, lastid)

# -----------------------------------------------------------------
