#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.clear Clear retrieved, crashed, cancelled and aborted simulations

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer
from pts.core.tools import logging, time, parsing

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("remote", nargs='?', default=None, help="the name of the remote host for which to clear simulations")
parser.add_argument("--status", type=parsing.string_list, help="clear/stop all simulations with these particular statuses")
parser.add_argument("--ids", type=parsing.simulation_ids, help="clear/stop the simulations with these ID's")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

if arguments.status is None: arguments.status = ["retrieved", "crashed", "aborted", "cancelled"]

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer.from_arguments(arguments)

# Run the synchronizer
synchronizer.run()

# -----------------------------------------------------------------
