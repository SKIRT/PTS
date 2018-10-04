#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.complete_data Complete columns that are specified by reference and save as a new file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.data import load_data
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "path of the data file")
definition.add_required("new_filename", "string", "path for the completed data file")
config = parse_arguments("upload_status", definition)

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the data file ...")

# Load
data = load_data(config.filename)

# -----------------------------------------------------------------

# Dereference
data.dereference()

# -----------------------------------------------------------------

# Inform the user
log.info("Writing the completed data file ...")

# Save
data.saveto(config.new_filename)

# -----------------------------------------------------------------
