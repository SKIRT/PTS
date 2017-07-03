#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.configure_skirt_item Test the configuration of a SKIRT simulation item.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema, show_parameters
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.configuration import ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("name", "string", "name of the configurable item")

# Parse
setter = ArgumentConfigurationSetter("configure_skirt_item")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

# Get the configuration parameters interactively
parameters, children = smile.prompt_parameters_for_type(config.name)

# Show the parameters
show_parameters(parameters, children)

# -----------------------------------------------------------------
