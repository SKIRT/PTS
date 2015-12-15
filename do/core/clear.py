#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.clear Clear retreived, crashed, cancelled and aborted simulations

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer

# -----------------------------------------------------------------

def string_list(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    return argument.split(",")

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("remote", nargs='?', default=None, help="(optional) the name of the remote host")
parser.add_argument("--status", type=string_list, help="clear all simulations with these particular statuses")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

if arguments.status is None: arguments.status = ["retreived", "crashed", "aborted", "cancelled"]

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer.from_arguments(arguments)

# Run the synchronizer
synchronizer.run()

# -----------------------------------------------------------------
