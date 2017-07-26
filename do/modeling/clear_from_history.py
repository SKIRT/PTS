#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clear_from_history Remove certain commands from the modeling history.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd, GalaxyModelingEnvironment
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Get config
definition = ConfigurationDefinition()
definition.add_required("commands", "string_list", "commands to remove from the history")
config = parse_arguments("clear_from_history", definition)

# -----------------------------------------------------------------

environment = GalaxyModelingEnvironment(modeling_path)
history = environment.history

# Make backup
#backup_path = fs.appended_filepath(environment.history_file_path, "_backup")
#fs.copy_file(environment.history_file_path, backup_path)

# -----------------------------------------------------------------

# Loop over the commands
for command in config.commands: history.remove_entries(command)

# Save
history.save()

# -----------------------------------------------------------------
