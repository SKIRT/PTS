#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.previous_config Show previous configuration.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.component.component import get_configuration_file_paths

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
filepaths = get_configuration_file_paths(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_required("match", "string", "(partial) command name")

# Get configuration
config = parse_arguments("previous_config", definition)

# -----------------------------------------------------------------

# TO BE IMPLEMENTED

# -----------------------------------------------------------------
