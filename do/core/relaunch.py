#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.relaunch Relaunch previous simulations

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer

# -----------------------------------------------------------------

def int_list(string):

    """
    This function returns a list of integer values, based on a string denoting a certain range (e.g. '3-9') or a
    set of integer values seperated by commas ('2,14,20')
    :param string:
    :return:
    """

    # Split the string
    splitted = string.split('-')

    if len(splitted) == 0: raise argparse.ArgumentError("No range given")
    elif len(splitted) == 1:

        splitted = splitted[0].split(",")

        # Check if the values are valid
        for value in splitted:
            if not value.isdigit(): raise argparse.ArgumentError("Argument contains unvalid characters")

        # Only leave unique values
        return list(set([int(value) for value in splitted]))

    elif len(splitted) == 2:

        if not (splitted[0].isdigit() and splitted[1].isdigit()): raise argparse.ArgumentError("Not a valid integer range")
        return range(int(splitted[0]), int(splitted[1])+1)

    else: raise argparse.ArgumentError("Values must be seperated by commas or by a '-' in the case of a range")

# -----------------------------------------------------------------

def simulation_ids(string):

    """
    This function ...
    :param string:
    :return:
    """

    # Initialize a dictionary
    delete = dict()

    # If the string is empty, raise an error
    if not string.strip(): raise argparse.ArgumentError("No input for argument")

    # Split the string by the ';' character, so that each part represents a different remote host
    for entry in string.split(";"):

        # Split again to get the host ID
        splitted = entry.split(":")

        # Get the host ID
        host_id = splitted[0]

        # Get the simulation ID's
        values = int_list(splitted[1])

        # Add the simulation ID's to the dictionary for the correspoding host ID
        delete[host_id] = values

    # Return the dictionary with ID's of simulations that should be deleted
    return delete

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("relaunch", type=simulation_ids, help="the ID's of the simulations to relaunch")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer.from_arguments(arguments)

# Run the synchronizer
synchronizer.run()

# -----------------------------------------------------------------
