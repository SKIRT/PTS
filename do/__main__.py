#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
# Main module for the do package
# -----------------------------------------------------------------

## \package do.__main__ Execute one of the scripts in this package from any current directory.
#
# Before proceeding, ensure that your login script includes the extra lines as described in \ref InstallSourceConf,
# and that you have logged in again after that change.
#
# To execute a python script named \c example.py residing in this directory, enter "pts example" in a Terminal window.
# This will work regardless of your current directory.
# You can include command line arguments as well, as in "pts example arg1 arg2".
# You can shorten the script name to the first few letters as long as there are no two scripts with matching names.
# For example "pts exa arg1 arg2" would still execute \c example.py, assuming that only one script has a name
# starting with the string "exa".
#
# Use "ipts" rather than "pts" to enter interactive python mode (with the >>> prompt) after executing the script.
# This is useful for testing or experimenting with pts functionality: the script imports the relevant
# pts module(s) and initializes some key objects that then can be used from the interactive python prompt.
#

# -----------------------------------------------------------------

# import standard modules
import sys
import os.path

# get a sorted list of the names (without extension) of all python scripts in the directory containing this file
path = os.path.dirname(sys.argv[0])
names = filter(lambda fn: fn.endswith(".py") and not fn.startswith("_"), os.listdir(path))
names = [ name[:-3] for name in sorted(names) ]

# get a list of the script names that match the first command line argument, if there is one;
# if the argument doesn't match any script names, try again with the portion of the script names after an underscore
if len(sys.argv) > 1:
    matches = filter(lambda name: name.startswith(sys.argv[1]), names)
    if len(matches)==0:
        matches = filter(lambda name: '_' in name and name.split('_',1)[1].startswith(sys.argv[1]), names)
else:
    matches = []

# if there is a unique match, execute the script
# (after adjusting the command line arguments so that it appears that the script was executed directly)
if len(matches)==1:
    target = os.path.join(path, matches[0] + ".py")
    sys.argv[0] = target
    del sys.argv[1]
    print "Executing: " + matches[0] + " " + " ".join(sys.argv[1:])
    exec open(target)

# if not, print the list of available scripts
else:
    print "Available scripts:"
    for name in names: print "  " + name

# -----------------------------------------------------------------
