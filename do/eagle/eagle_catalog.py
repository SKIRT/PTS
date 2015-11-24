#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_catalog Inspect or build the catalog for a given EAGLE snapshot.
#
# This script lists the basic properties of the galaxies in the catalog for the specified EAGLE snapshot.
# If the catalog does not yet exist, the script offers to construct it.
#
# The script expects one or two command-line arguments. The first argument specifies the EAGLE snapshot
# (always at redshift zero) using one of the shorthands listed in eagle.config. If the catalog for this
# snapshot exists, the second argument is not used (i.e. all galaxies in the catalog are listed).
#
# If the catalog does not exist, and the user requests that it be created, the second argument specifies
# the minimum stellar mass within a 30kpc aperture (in solar mass units) for a galaxy to be included in
# the catalog.

# -----------------------------------------------------------------

import sys
import os.path
from eagle.galaxy import Snapshot

# -----------------------------------------------------------------

# parse the arguments
if not len(sys.argv) in (2,3): raise ValueError("This script expects one or two command-line arguments")
snapshotname = sys.argv[1]

# open the snapshot files
print "Opening the snapshot..."
snap = Snapshot(snapshotname)
snap.printinfo()
print ""

# if the catalog exists, show info on its contents
if os.path.isfile(snap.catalogfilepath()):
    galaxies = snap.galaxies()
    galaxies.printinfo()

# otherwise, ask the user whether to construct it
else:
    proceed = raw_input("--> Would you like to build the catalog for {0} at redshift {1}? [y/n] "  \
                        .format(snap.eaglesim, snap.orig_redshift))
    if proceed.lower().startswith("y"):
        if len(sys.argv) != 3: raise ValueError("This script expects a second command-line argument")
        minstarmass = float(sys.argv[2])
        snap.exportcatalog(minstarmass)

# -----------------------------------------------------------------
