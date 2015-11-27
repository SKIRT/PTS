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
from operator import itemgetter

# -----------------------------------------------------------------

# Get the path to the 'do' directory
path = os.path.dirname(sys.argv[0])

# Loop over the directories within the 'do' subpackage
scripts = []
for item in os.listdir(path):

    # Get the full path to the item
    item_path = os.path.join(path, item)

    # Skip items that are not directories
    if not os.path.isdir(item_path): continue

    # Loop over the files in the directory
    for name in os.listdir(item_path):

        # Get the full name to the file
        file_path = os.path.join(item_path, name)

        # Skip items that are not files
        if not os.path.isfile(file_path): continue

        # Add the file path
        if file_path.endswith(".py") and "__init__" not in file_path: scripts.append((item, name))

# Sort the script names
scripts = sorted(scripts, key=itemgetter(1))

# Get the name of the script to execute
script_name = sys.argv[1] if len(sys.argv) > 1 else None

# Get a list of the script names that match the first command line argument, if there is one
if "/" in script_name:

    matches = []
    dir_name = script_name.split("/")[0]
    script_name = script_name.split("/")[1]

    # Loop over all found items
    for item in scripts:
        if item[0] == dir_name and item[1].startswith(script_name): matches.append(item)

elif script_name is not None: matches = filter(lambda item: item[1].startswith(script_name), scripts)
else: matches = []

# If there is a unique match, execute the script
# (after adjusting the command line arguments so that it appears that the script was executed directly)
if len(matches) == 1:
    target = os.path.join(path, matches[0][0], matches[0][1])
    sys.argv[0] = target
    del sys.argv[1]
    print "Executing: " + matches[0][0] + "/" + matches[0][1] + " " + " ".join(sys.argv[1:])
    exec open(target)

elif len(matches) == 0:

    print "No match found. Available scripts:"

    # Sort on the 'do' subfolder name
    scripts = sorted(scripts, key=itemgetter(0))

    current_dir = None
    for script in scripts:

        if current_dir == script[0]:
            print " "*len(current_dir) + "/" + script[1]
        else:
            print script[0] + "/" + script[1]
            current_dir = script[0]

else:

    print "The command you provided is ambiguous. Possible matches:"

    # Sort on the 'do' subfolder name
    matches = sorted(matches, key=itemgetter(0))

    current_dir = None
    for script in matches:

        if current_dir == script[0]:
            print " "*len(current_dir) + "/" + script[1]
        else:
            print script[0] + "/" + script[1]
            current_dir = script[0]

# -----------------------------------------------------------------
