#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_catalog Inspect or build the catalog for the default EAGLE snapshot.
#
# This script lists the basic properties of the galaxies in the catalog for the EAGLE snapshot
# that has been configured as the default snapshot in eagle.config.
# If the catalog does not yet exist, the script offers to construct it.
#
# The script expects either zero or exactly two command-line arguments specifying a minimum and a maximum gas mass
# (in solar mass units). If these values are specified, only galaxies are listed with a gas mass inside the range.
#

# -----------------------------------------------------------------

import sys
import os.path
from eagle.galaxy import Snapshot

# -----------------------------------------------------------------

# open the snapshot files
print "Opening the snapshot..."
snap = Snapshot()
snap.printinfo()
print ""

# if the catalog exists, show info on its contents
if os.path.isfile(snap.catalogfilepath()):
    galaxies = snap.galaxies()
    # if appropriate command-line arguments are provided, filter the galaxies
    if not (len(sys.argv) in (1,3)): raise ValueError("This script expects zero or two command-line arguments")
    if len(sys.argv) == 3:
        mingasmass = float(sys.argv[1])
        maxgasmass = float(sys.argv[2])
        print "Restricting gas mass to range [{0:.1e},{1:.1e}]".format(mingasmass,maxgasmass)
        print ""
        galaxies.remove_gasmass_below(mingasmass)
        galaxies.remove_gasmass_above(maxgasmass)
    galaxies.printinfo()

# otherwise, ask the user whether to construct it
else:
    proceed = raw_input("Do you want to build the catalog for {0} at redshift {1}? (y/n) "  \
                        .format(config.default_eaglesim, config.default_redshift))
    if proceed.lower().startswith("y"):
        snap.exportcatalog()

# -----------------------------------------------------------------
