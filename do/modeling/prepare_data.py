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

# Import the relevant PTS classes and modules
from pts.modeling.preparation.datapreparation import DataPreparer
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add required arguments
config.add_required("image", str, "the name of the image for which to run the preparation")

# Add optional arguments
config.add_optional("reference", str, "the name of the reference image")
config.add_flag("steps", "write the results of intermediate steps")
config.add_flag("visualise", "make visualisations")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting prepare_data ...")

# -----------------------------------------------------------------

# Create a DataPreparer instance
preparer = DataPreparer(config.get_settings())

# Run the data preparation
preparer.run()

# -----------------------------------------------------------------
