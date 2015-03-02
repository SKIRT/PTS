#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_build Build visualization files for completed EAGLE SKIRT-runs.
#
# This script creates compressed zip archives for the \c in and \c out directories of completed EAGLE SKIRT-runs,
# copies these archives to an offsite archive location, and marks the SKIRT-runs as archived in the database.
# As an extra safeguard in case the database would be damaged, the contents of the relevant database record
# is stored in a text file in the SKIRT-run's \c in directory before the zip archive is created.
#
# If the script is invoked without command-line arguments, it will handle all completed SKIRT-runs that have
# not yet been archived. Otherwise, the (single) command-line argument should specify a comma-seperated list
# of run-ids and/or run-id ranges (expressed as two run-ids with a dash in between) for which the archival
# operation should be performed. The specified SKIRT-runs must have a run-status of 'completed' or 'archived'.
#

# -----------------------------------------------------------------

import sys
import os
import os.path
import subprocess
import eagle.config as config
from eagle.database import Database
from eagle.skirtrun import SkirtRun
from eagle.skirtrun import runids_in_range

# -----------------------------------------------------------------
# ==== path definition that may be moved elsewhere over time   ====

# offsite archive path
archive_path = "/Volumes/astro_skirt/pcamps/EAGLE_Archive"

# -----------------------------------------------------------------

# returns a list of run-ids for all currently completed and not-yet archived SKIRT-runs
def completed_runids():
    db = Database()
    runids = sorted([ row['runid'] for row in db.select("runstatus = 'completed'") ])
    db.close()
    return runids

# returns a list of run-ids for the SKIRT-runs in the specified range, restricted to those that
# are currently completed and/or archived
def completed_or_archived_runids_in_range(runidspec):
    runids = sorted(runids_in_range(runidspec))
    db = Database()
    runids = filter(lambda runid: db.select("runid="+str(runid))[0]['runstatus'] in ('completed','archived'), runids)
    db.close()
    return runids

# -----------------------------------------------------------------

# construct a list of SkirtRun objects to be processed, depending on the command line
if len(sys.argv) > 1:
    runids = completed_or_archived_runids_in_range(sys.argv[1])
else:
    runids = completed_runids()

# verify that the offsite archive server is mounted
if len(runids)==0:
    print "No SKIRT-runs to be archived."
    exit()
if not os.path.isdir(archive_path):
    print "Please mount the offsite archive server and try again."
    exit()

# open the database and process each SKIRT-run
db = Database()
for runid in runids:
    print "--> Creating archives for SKIRT-run {}...".format(runid)

    # get the paths for the in and out directories, and the corresponding zip archives
    skirtrun = SkirtRun(runid)
    indir = skirtrun.inpath()
    outdir = skirtrun.outpath()
    inzip = indir + ".zip"
    outzip = outdir + ".zip"

    # export the database record
    info = db.select("runid="+str(runid))[0]
    infofile = open(os.path.join(indir,"database_record.txt"), 'w')
    infofile.write('# SKIRT-run database record\n')
    for key in info.keys():
        infofile.write( ("{} = {}\n").format(key, info[key]) )
    infofile.close()

    # create the archives
    if os.path.isfile(inzip): os.remove(inzip)
    error = subprocess.call(("zip", "-rj", inzip, indir))
    if error: raise ValueError("Error while zipping the in directory: " + str(error))
    if os.path.isfile(outzip): os.remove(outzip)
    error = subprocess.call(("zip", "-rj", outzip, outdir))
    if error: raise ValueError("Error while zipping the out directory: " + str(error))

    # copy the archives offsite
    print "--> Copying archives for SKIRT-run {}...".format(runid)
    timestamp = config.timestamp()
    archdir = os.path.join(archive_path, "g-{0:04}".format(runid//1000))
    archin  = os.path.join(archdir, "r-{0:04}-{1:03}_in_{2}.zip".format(runid//1000, runid%1000, timestamp))
    archout = os.path.join(archdir, "r-{0:04}-{1:03}_out_{2}.zip".format(runid//1000, runid%1000, timestamp))
    if not os.path.isdir(archdir): os.makedirs(archdir)
    error = subprocess.call(("cp", "-v", inzip, archin))
    if error: raise ValueError("Error while copying the in archive: " + str(error))
    error = subprocess.call(("cp", "-v", outzip, archout))
    if error: raise ValueError("Error while copying the out archive: " + str(error))

    # update the runstatus in the database
    db.updatestatus((runid,), 'archived')
    db.commit()

# close the database
db.close()

print "Done..."

# -----------------------------------------------------------------
