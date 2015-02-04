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
# The script expects up to three command-line arguments specifying a minimum number of particles (for both stars and
# gas) and a minimum and a maximum stellar mass (in solar mass units). If these values are specified, only galaxies
# within these constraints are listed.
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
    if len(sys.argv) > 1:
        minparticles = int(sys.argv[1])
        print "Requiring a minimum of {0} particles for stars and gas".format(minparticles)
        galaxies.remove_starparticles_below(minparticles)
        galaxies.remove_gasparticles_below(minparticles)
    if len(sys.argv) > 2:
        minstarmass = float(sys.argv[2])
        print "Requiring a minimum stellar mass of {0:.1e} solar masses".format(minstarmass)
        galaxies.remove_starmass_below(minstarmass)
    if len(sys.argv) > 3:
        maxstarmass = float(sys.argv[3])
        print "Requiring a maximum stellar mass of {0:.1e} solar masses".format(maxstarmass)
        galaxies.remove_starmass_above(maxstarmass)
    galaxies.printinfo()

# otherwise, ask the user whether to construct it
else:
    proceed = raw_input("--> Would you like to build the catalog for {0} at redshift {1}? [y/n] "  \
                        .format(snap.eaglesim, snap.orig_redshift))
    if proceed.lower().startswith("y"):
        snap.exportcatalog()

# -----------------------------------------------------------------
