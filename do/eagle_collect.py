#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_build Build visualization files for completed EAGLE SKIRT-runs.
#
# This script collects the statistics contained in the info files for a set of SKIRT-runs into a single data file,
# which is placed in the "Collections" directory specified in eagle.config.
# The current implementation simply selects the set of SKIRT-runs based on the EAGLE simulation name specified
# as the sole command line argument for the script, for example "Ref25"
#

# -----------------------------------------------------------------

import sys
import types
import os.path
import pickle
import numpy as np
import pts.archive as arch
import eagle.config as config
from eagle.database import Database
from eagle.skirtrun import SkirtRun

# -----------------------------------------------------------------

# returns a list of SkirtRun objects corresponding to all completed or archived skirt-runs
# for the given eagle simulation name, in order of run-id.
def completed_skirtruns_for_simulation(eaglesim):
    db = Database()
    query = "runstatus in ('completed','archived') and eaglesim='{}'".format(eaglesim)
    runids = sorted([ row['runid'] for row in db.select(query) ])
    runs = [ SkirtRun(runid) for runid in runids ]
    db.close()
    return runs

# collects the statistics for a given list of skirt-runs into a single dictionary with a numpy array for each entry
# and dumps the dictionary in pickle format to the specified file.
def collect_info(skirtruns, outfilepath):
    # collect the info in single dict
    collection = { }
    numruns = len(skirtruns)
    for runindex in range(numruns):
        vispath = skirtruns[runindex].vispath()
        infofile = arch.listdir(vispath, "_info.txt")[0]
        for line in arch.opentext(os.path.join(vispath,infofile)):
            if not line.startswith("#"):
                key,dummy,value = line.split(None, 2)
                if not key in collection:
                    collection[key] = np.zeros(numruns)
                collection[key][runindex] = float(value)

    # dump it into file
    outfile = open(outfilepath, 'w')
    pickle.dump(collection, outfile)
    outfile.close()
    print "Created info collection " + outfilepath

# -----------------------------------------------------------------

# get the command-line argument
if len(sys.argv) != 2: raise ValueError("This script expects exactly one command-line argument")
eaglesim = sys.argv[1]

# get a list of skirt-runs for which to collect statistics
skirtruns = completed_skirtruns_for_simulation(eaglesim)

# perform the collection
if len(skirtruns) > 0:
    print "Collecting statistics from {} SKIRT-runs for EAGLE simulation {}...".format(len(skirtruns),eaglesim)
    collect_info(skirtruns, os.path.join(config.collect_path,"{}_info_collection.dat".format(eaglesim)))

else:
    print "No SKIRT-runs for EAGLE simulation {}.".format(eaglesim)

# -----------------------------------------------------------------
