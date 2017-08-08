#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.plotresources Plot a histogram of the resources used for a selection of SKIRT-run simulations.
#
# This script plots two histograms of the resources used for the SKIRT simulations corresponding to
# a set of SKIRT-run records:
#  - The computation time, determined as the product of the wall-time and the number of parallel processes.
#  - The peak memory usage of the root process
#
# The script expects a single command-line argument specifying the label for the records to be considered.
#

# -----------------------------------------------------------------

# Import standard modules
import os.path
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.eagle import config
from pts.eagle import database
from pts.eagle.skirtrun import SkirtRun

# -----------------------------------------------------------------

def advance(db, label, stagefrom, stageto):
    with db.transaction():
        # get the eligible records
        records = db.select("label=? and stage=? and status='succeeded'", (label, stagefrom))
        if len(records) > 0:

            # show some info and ask for confirmation
            print "There are {} successful SKIRT-runs with label {} at stage {}.".format(len(records), label, stagefrom)
            confirm = raw_input("--> Would you like to advance these to stage {}? [y/n]: ".format(stageto))

            # advance stage if user confirmed
            if confirm.lower().startswith('y'):
                db.updatestage(records, stageto)
                print "Advancing {} SKIRT-runs from stage {} to {}".format(len(records), stagefrom, stageto)
            else:
                print "Update was rejected"

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv)!=2: raise ValueError("This script expects a single command-line argument: label")
label = sys.argv[1]

# get the eligible records
db = database.Database()
records = db.select("((stage='simulate' and status='succeeded') or (stage in ('observe', 'store', 'completed')))" \
                    " and label=?", (label,))
db.close()
size = len(records)
if size==0: raise ValueError("There are no simulated records with label " + label)

# assemble the statistics for all records
log.info("Assembling statistics for {} SKIRT simulations...".format(size))
time = np.zeros(size)
memory = np.zeros(size)
for index in range(size):
    logfilepath = SkirtRun(records[index]["runid"]).simulation().logfilepath()
    for line in open(logfilepath):
        if " Finished simulation " in line:
            segments = line.split()
            processes = float(segments[segments.index("processes")-1])
            timeindex = segments.index("s") if "s" in segments else segments.index("s.")
            walltime = float(segments[timeindex-1])
            time[index] = processes * walltime
        if " Available memory: " in line:
            segments = line.split()
            memory[index] = float(segments[segments.index("usage:")+1])

# construct the time histogram
log.info("Constructing plots...")
figure = plt.figure(figsize=(8,8))
plt.xlabel('Computation time (min)', fontsize='medium')
plt.ylabel('Simulation count', fontsize='medium')
plt.hist(time/60, bins=25)
plt.vlines(np.mean(time)/60, 0, 5, colors='r')
plotfilepath = os.path.join(config.plots_path, label + "_time_hist.pdf")
plt.savefig(plotfilepath)
log.info("Created time histogram " + plotfilepath)

# construct the memory histogram
figure = plt.figure(figsize=(8,8))
plt.xlabel('Peak memory usage (GB)', fontsize='medium')
plt.ylabel('Simulation count', fontsize='medium')
plt.hist(memory, bins=25)
plt.vlines(np.max(memory), 0, 5, colors='r')
plotfilepath = os.path.join(config.plots_path, label + "_memory_hist.pdf")
plt.savefig(plotfilepath)
log.info("Created memory histogram " + plotfilepath)

# -----------------------------------------------------------------
