#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.collect Collect statistics contained in the info files for a set of SKIRT-runs.
#
# This script collects the statistics contained in the info files for a set of SKIRT-runs into a single data file,
# which is placed in the "Collections" directory specified in eagle.config.
# The script expects one or more command line arguments, which are used as possible values for the label
# and eaglesim fields in the SKIRT-run database. For example, "pts eagle/collect Recal25 convergence" will collect all
# runs that have a label value of "Recal25" or "convergence" \em and an eaglesim value of "Recal25" or "convergence".
#

# -----------------------------------------------------------------

# Import standard modules
import sys
import os.path
import pickle
import numpy as np

# Import the relevant PTS classes and modules
import pts.core.tools.archive as arch
from pts.eagle import config
from pts.eagle.database import Database
from pts.eagle.skirtrun import SkirtRun

# -----------------------------------------------------------------

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

# chain the command-line arguments into a query list
if len(sys.argv) <= 1: raise ValueError("This script expects one or more command-line arguments")
querylist = "('{}')".format("','".join(sys.argv[1:]))
namelist = "_".join(sys.argv[1:])

# get a list of SkirtRun objects for which to collect statistics, in order of run-id
db = Database()
query = "((stage='observe' and status='succeeded') or (stage in ('store', 'completed')))" \
        " and label in {0} and eaglesim in {0}".format(querylist)
runids = sorted([ row['runid'] for row in db.select(query) ])
skirtruns = [ SkirtRun(runid) for runid in runids ]
db.close()

# perform the collection
if len(skirtruns) > 0:
    print "Collecting statistics from {} SKIRT-runs with label and eaglesim fields in {}...".format(len(skirtruns),querylist)
    collect_info(skirtruns, os.path.join(config.collections_path,"{}_info_collection.dat".format(namelist)))

else:
    print "There are no SKIRT-runs with label and eaglesim fields in {}.".format(querylist)

# -----------------------------------------------------------------
