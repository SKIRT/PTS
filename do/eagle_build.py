#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_build Build visualization files for completed EAGLE SKIRT-runs.
#
# This script builds RGB images, wavelength movies, or info files from the results of completed EAGLE SKIRT-runs.
# The visualization files are placed in the corresponding \c vis directories, replacing any previous versions.
#
# The script expects exactly two command-line arguments. The first argument specifies what type of files to build.
# It should be one of the following strings (it is sufficient to specify an unambiguous portion of the string):
#  - rgbimages: an optical RGB image for each of the instruments
#  - wavemovie: a wavelength movie for the three instruments
#  - infofile: a text info file with statistics on the simulation results
#
# The second argument further specifies the operation that should be performed. It should be one of the following:
#  - the string "update" to perform the build for all completed SKIRT-runs for which the visualization files do not
#    yet exist.
#  - the string "rebuild" to perform the build for any and all all completed SKIRT-runs, even if the visualization
#    files already exist. Previous versions are replaced.
#  - a comma-seperated list of run-ids and/or run-id ranges (expressed as two run-ids with a dash in between)
#    for which the build should be performed. Visualization files will be rebuilt even if they already exist.
#    The script does not verify whether the specified SKIRT-runs have been completed.
#

# -----------------------------------------------------------------

import sys
import types
import os.path
from eagle.database import Database
from eagle.skirtrun import SkirtRun

# -----------------------------------------------------------------
# ==== utility functions that may be moved elsewhere over time ====

# returns a list of run-ids corresponding to the specified range string, or None in case of syntax error
def runids_in_range(runidspec):
    try:
        runids = []
        for segment in runidspec.split(","):
            if "-" in segment:
                first,last = map(int,segment.split("-"))
                runids += [ id for id in range(first,last+1) ]
            else:
                if segment!="": runids += [ int(segment) ]
        return runids
    except Exception:
        return None

# returns a list of SkirtRun objects corresponding to all completed or archived skirt-runs, in order of run-id,
# optionally omitting any skirt-runs for which all files in the specified sequence exist in the visualization folder
def completed_skirtruns(unless_filenames=None):
    db = Database()
    runids = sorted([ row['runid'] for row in db.select("runstatus='completed' or runstatus='archived'") ])
    runs = [ SkirtRun(runid) for runid in runids ]
    if unless_filenames!=None:
        runs = filter(lambda run: not has_visualization_files(run,unless_filenames), runs)
    db.close()
    return runs

# returns True if all files in the specified sequence exist in the visualization folder of the specified skirt-run
def has_visualization_files(skirtrun, filenames):
    filenames = filenames if isinstance(filenames, (types.TupleType,types.ListType)) else [ filenames ]
    vispath = skirtrun.vispath()
    prefix = skirtrun.prefix()
    for filename in filenames:
        if not os.path.isfile(os.path.join(vispath, prefix+"_"+filename)): return False
    return True

# move any of the files in the specified sequence from the skirt-run output directory to the visualization directory
def move_visualization_files(skirtrun, filenames):
    for visfile in filter(lambda fn: fn.endswith(filenames), os.listdir(skirtrun.outpath())):
        os.rename(os.path.join(skirtrun.outpath(), visfile),
                  os.path.join(skirtrun.vispath(), visfile))

# -----------------------------------------------------------------

# get and parse the command-line arguments
if len(sys.argv) != 3: raise ValueError("This script expects exactly two command-line arguments")

vistype = sys.argv[1].lower()
rgb  = "rgbimages".startswith(vistype)
wave = "wavemovie".startswith(vistype)
info = "infofile".startswith(vistype)
if not rgb and not wave and not info: raise ValueError("Unknown visualization type: " + vistype)

buildrange = sys.argv[2].lower()
update =  "update".startswith(buildrange)
rebuild = "rebuild".startswith(buildrange)
runidrange = runids_in_range(buildrange)
if not update and not rebuild and not runidrange: raise ValueError("Unknown build range: " + buildrange)

# -----------------------------------------------------------------

# construct a list of relevant filename endings depending on the visualization type
filenames = []
if rgb: filenames += [ "xy_total.png", "xz_total.png", "yz_total.png" ]
if wave: filenames += [ "wave.mov" ]
if info: filenames += [ "info.txt" ]
filenames = tuple(filenames)

# construct the list of SKIRT-runs to be processed
if runidrange:
    skirtruns = [ SkirtRun(runid) for runid in runidrange ]
else:
    skirtruns = completed_skirtruns(filenames if update else None)

# build RGB images for each SKIRT-run
if rgb:
    from pts.makergbimages import makergbimages
    print "Building RGB images for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building RGB images for SKIRT-run {}...".format(skirtrun.runid())
        makergbimages(skirtrun.simulation(), wavelength_tuples=((0.75,0.60,0.45),) )
        move_visualization_files(skirtrun, filenames)

# build wavelength movies for each SKIRT-run
if wave:
    from pts.makewavemovie import makewavemovie
    print "Building wavelength movies for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building wavelength movie for SKIRT-run {}...".format(skirtrun.runid())
        makewavemovie(skirtrun.simulation())
        move_visualization_files(skirtrun, filenames)

# build info files for each SKIRT-run
if info:
    # from pts.makeinfofile import makeinfofile
    print "Building info files for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building info file for SKIRT-run {}...".format(skirtrun.runid())
        # makeinfofile(skirtrun)
        move_visualization_files(skirtrun, filenames)

print "Done..."

# -----------------------------------------------------------------
