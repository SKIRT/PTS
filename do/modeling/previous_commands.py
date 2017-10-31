#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.previous_commands Show previous modeling commands.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.component.component import load_modeling_commands

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
commands = load_modeling_commands(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_positional_optional("match", "string", "(partial) command name")

# Get configuration
config = parse_arguments("previous_commands", definition)

# -----------------------------------------------------------------

# Loop over the commands
for command, line in zip(commands.as_commands(), commands):

    if config.match is not None and not command.command.startswith(config.match): continue
    print(line)

# -----------------------------------------------------------------
