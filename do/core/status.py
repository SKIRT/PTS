#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.status Check the status of SKIRT simulations running remotely.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.synchronizer import RemoteSynchronizer
from pts.core.tools import logging, time

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('remote', nargs='?', default=None, help="(optional) the name of the remote host")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")
parser.add_argument("--report", action="store_true", help="write a report file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = os.path.join(os.getcwd(), time.unique_name("status") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting status ...")

# -----------------------------------------------------------------

# Create a RemoteSynchronizer instance
synchronizer = RemoteSynchronizer.from_arguments(arguments)

# Run the synchronizer
synchronizer.run()

# -----------------------------------------------------------------
