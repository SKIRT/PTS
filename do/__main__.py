#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
# Main module for the do package
# -----------------------------------------------------------------

## \package pts.do.__main__ Execute one of the scripts in this package from any current directory.
#
# Before proceeding, ensure that your login script includes the extra lines as described in \ref InstallMacSetUp,
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

# Import standard modules
import sys
import os.path

# Import the relevant PTS modules
from pts.core.tools import inspection

# -----------------------------------------------------------------

# Get the name of the script to execute
script_name = sys.argv[1] if len(sys.argv) > 1 else None

# Find matching scripts
match = inspection.find_matching_script(script_name)
if match is None: exit()

# Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
target = os.path.join(inspection.pts_do_dir, match[0], match[1])
sys.argv[0] = target
del sys.argv[1]
print "Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:])
exec open(target)

# -----------------------------------------------------------------
