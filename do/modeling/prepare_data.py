#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.preparedata Do the preparation step of the SKIRT radiative transfer modeling procedure.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.preparation.datapreparation import DataPreparer
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("image", type=str, nargs='?', help="the name of the image for which to run the preparation")
parser.add_argument("--reference", type=str, help="the name of the reference image")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help="write a report file")

parser.add_argument("--steps", action="store_true", help="write the results of intermediate steps")

# Configuration
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Visualisation
parser.add_argument("--visualise", action="store_true", help="make visualisations")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
arguments.path = fs.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(arguments.path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting prepare_data ...")

# -----------------------------------------------------------------

# Create a DataPreparer instance
preparer = DataPreparer.from_arguments(arguments)

# Run the data preparation
preparer.run()

# -----------------------------------------------------------------
