#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.status Check the status of SKIRT simulations running remotely

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer

# -----------------------------------------------------------------

# This function returns a list of integer values, based on a string denoting a certain range (e.g. '3-9') or a
# set of integer values seperated by commas ('2,14,20')
def int_list(string):

    # Split the string
    splitted = string.split('-')

    if len(splitted) == 0: raise argparse.ArgumentTypeError("No range given")
    elif len(splitted) == 1:

        splitted = splitted[0].split(",")

        # Check if the values are valid
        for value in splitted:
            if not value.isdigit(): raise argparse.ArgumentTypeError("Argument contains unvalid characters")

        # Only leave unique values
        return list(set([int(value) for value in splitted]))

    elif len(splitted) == 2:

        if not (splitted[0].isdigit() and splitted[1].isdigit()): raise argparse.ArgumentTypeError("Not a valid integer range")
        return range(int(splitted[0]), int(splitted[1])+1)

    else: raise argparse.ArgumentTypeError("Values must be seperated by commas or by a '-' in the case of a range")

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--delete', type=int_list)

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer()

# Run the synchronizer
synchronizer.run(arguments.delete)

# -----------------------------------------------------------------
