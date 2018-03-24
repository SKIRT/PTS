#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.prep.composer import oligo_or_pan, pan_type

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Name and type of simulation
definition.add_positional_optional("name", "string", "name for the model")
definition.add_optional("type", "string", "simulation type", pan_type, choices=oligo_or_pan)

# -----------------------------------------------------------------

# Ski file path
definition.add_optional("from_ski", "file_path", "ski file path")

# -----------------------------------------------------------------

# Commands to be run
definition.add_optional("commands", "string_list", "commands to be run in interactive mode")

# -----------------------------------------------------------------
