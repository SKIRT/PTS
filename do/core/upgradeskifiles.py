#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.upgradeskifiles Upgrade ski files to the latest format version.
#
# This script examines a set of ski files and, if needed, upgrades each file the latest format version.
#
# If the script is invoked without command line arguments, it examines all ski files in the current directory.
# If there is a command line argument, the script examines all ski files in the specified directory.
#
# If the first command line argument is the string "-r" then the script also examines all subdirectories,
# recursively, of either the current directory (if there is no second argument) or the directory specified
# in the second argument.
#
# If one of the examined ski files indeed needs an upgrade, a backup copy of the original file is placed next to it,
# and the original file is overwritten with the upgraded version.
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import sys

# Import the relevant PTS classes and modules
from pts.core.prep.upgradeskifile import upgradeskifile

# -----------------------------------------------------------------

# function to examine the ski files in a single directory
def do_single_dir(target):
    # get a sorted list of the names of all ski files in the target directory
    names = sorted(filter(lambda fn: fn.endswith(".ski"), os.listdir(target)))

    # examine each ski file
    print "-> " + target + "..."
    for name in names:
        changed = upgradeskifile(os.path.join(target,name))
        if changed:
            print "** " + name + " (UPGRADED)."
        else:
            print "   " + name + " (not changed)."

# -----------------------------------------------------------------

# determine if recursive decent is requested
recursive = True if len(sys.argv) > 1 and sys.argv[1]=="-r" else False

# get the target directory
index = 2 if recursive else 1
target = sys.argv[index] if len(sys.argv) > index else ""
target = os.path.realpath(os.path.expanduser(target))

# examine all ski files as requested
if recursive:
    # iterate over all nested directories
    for dirpath, dirs, files in os.walk(target):
        do_single_dir(dirpath)
else:
    do_single_dir(target)

print "All done"

# -----------------------------------------------------------------
