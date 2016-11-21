#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.perform Perform the process for a given stage on a set of SKIRT-runs records.
#
# This script performs the process for a given stage on a set of SKIRT-runs records, in one of two modes
# depending on the command line arguments:
#
#    pts eagle/perform <stage> loop <runtime>
#    pts eagle/perform <stage> force <runidspec>
#
# where:
#  - <stage> is one of "extract", "simulate", "observe", or "store";
#  - <runtime> is a run time in seconds (given as a floating point expression);
#  - <runidspec> is a comma-seperated list of run-ids and/or run-id ranges (two run-ids with a dash in between).
#
# In loop mode, the script repeatedly performs the specified stage for a SKIRT-run record that is
# at the corresponding workflow stage and has the 'scheduled' status, until there are no more such records,
# or until the specified run time (given in seconds) has been surpassed.
#
# In force mode, the script performs the specified stage for each of the specified SKIRT-run records,
# regardless of the stage and status of those records.
#

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
# (the appropriate module for performing the requested stage is imported conditionally below)
from pts.core.tools.logging import log
from pts.eagle import performer

# -----------------------------------------------------------------

# get the command-line arguments
if len(sys.argv) != 4: raise ValueError("This script expects three command-line arguments: " \
                                        "(stage loop runtime) or (stage force runidspec)")
stage = sys.argv[1]
if not stage in ("extract", "simulate", "observe", "store"):
    raise ValueError("First argument must be 'extract', 'simulate', 'observe', or 'store'")
mode = sys.argv[2]
if not mode in ("loop", "force"):
    raise ValueError("Second argument must be 'loop' or 'force'")
argum = sys.argv[3]

# get appropriate callback for given stage
if stage == "extract":
    from pts.eagle.extractor import extract as callback
if stage == "simulate":
    from pts.eagle.simulator import simulate as callback
if stage == "observe":
    log.error("Sorry - the observe stage is not yet implemented")
    exit()
if stage == "store":
    log.error("Sorry - the store stage is not yet implemented")
    exit()

# perform according to requested mode
log.info("Performing {} in {} mode ...".format(stage, mode))
if mode == "loop":
    performer.loop(callback, stage, float(eval(argum)))
if mode == "force":
    performer.force(callback, stage, argum)
log.info("Finished performing.")

# -----------------------------------------------------------------
