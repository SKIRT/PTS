#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.history View the history of issued modeling commands.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.modeling.core.component import load_modeling_history

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Get configuration
setter = ArgumentConfigurationSetter("history")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting history ...")

# -----------------------------------------------------------------

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

# Get the history
history = load_modeling_history(fs.cwd())

# Loop over the entries
for i in range(len(history)):
    print(history["Command"][i], history["Start time"][i], history["End time"][i])

# -----------------------------------------------------------------
