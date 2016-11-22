#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.comparecollections Compare the magnitudes of galaxies in two EAGLE result collections
#
# This script compares the magnitudes of corresponding galaxies in two EAGLE result collections,
# and lists the mean and maximum differences for each property. If the maximum difference for a property
# is smaller than the specified epsilon (optional), the statistics for that property are not listed.
# The collections must already have been created (and placed in the default Collections path) using
# the eagle/collect script.
#
# The script expects exactly three command-line arguments: two collection names (without path nor
# filename extension) and the magnitude epsilon.
#

# ----------------------------------------------------------------------

# Import standard modules
import numpy as np
import os.path
import sys

# Import the relevant PTS classes and modules
import pts.eagle.config as config
from pts.eagle.collection import Collection

# ----------------------------------------------------------------------

# get the command-line arguments
if len(sys.argv)!=4: raise ValueError("This script expects 3 command-line arguments: " \
                                    + "collectionA collectionB magnitude-epsilon")
name1 = sys.argv[1]
name2 = sys.argv[2]
eps = float(sys.argv[3])

# load the collections and print their names
c1 = Collection(name1)
c2 = Collection(name2)
print "Collections:"
print "  A: " + c1.name()
print "  B: " + c2.name()

# find common sets of galaxies and properties and report on differences
ids1 = c1.galaxy_ids()
ids2 = c2.galaxy_ids()
ids = sorted(ids1 & ids2)
props1 = set(filter(lambda key: "_magnitude_" in key, c1.property_names()))
props2 = set(filter(lambda key: "_magnitude_" in key, c2.property_names()))
props = sorted(props1 & props2)
print "Composition:"
print "  {} common galaxies; {} extra galaxies in A; {} extra galaxies in B" \
      .format(len(ids), len(ids1) - len(ids), len(ids2) - len(ids))
print "  {} common properties; {} extra properties in A; {} extra properties in B" \
      .format(len(props), len(props1) - len(props), len(props2) - len(props))

# collection of galaxy ids for which the maximum difference is exceeded
offenders = set()

# report on the differences for each property
maxproplen = max(map(len,props))
print "Differences:" + (" "*(maxproplen-6))+"mean    std    max   nan"
for prop in props:
    # get the values
    values1 = c1.property_values(prop, ids)
    values2 = c2.property_values(prop, ids)
    # locate non-detections (indicated by NaN) and treat them seperately
    nan1 = np.isnan(values1)
    nan2 = np.isnan(values2)
    nandiff = np.count_nonzero(np.logical_xor(nan1,nan2))
    # calculate the numeric differences
    diff = np.abs(values1-values2)
    # build and print the report string
    if nandiff > 0 or np.nanmax(diff) > eps:
        message = "  {:"+str(maxproplen)+"}: {:6.3f} {:6.3f} {:6.3f}"
        if nandiff>0: message += "   {}"
        print message.format(prop, np.nanmean(diff), np.nanstd(diff), np.nanmax(diff), nandiff)

        # remember the galaxy ids for which the maximum difference was exceeded
        diff_notnan = diff[~np.isnan(diff)]
        ids_notnan = np.array(ids)[~np.isnan(diff)]
        offenders.update(ids_notnan[diff_notnan>eps])

# report the galaxy ids for which the maximum difference was exceeded
print "Offender galaxy IDs:"
print sorted(map(int,offenders))

# ----------------------------------------------------------------------
