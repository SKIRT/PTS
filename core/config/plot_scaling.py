#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.plot.scaling import scaling_properties

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Flags
definition.add_flag("recursive", "look for simulation in directories recursively", True)

#
definition.add_positional_optional("properties", "string_list", "plot the scaling of these properties", choices=scaling_properties, default=scaling_properties)

# -----------------------------------------------------------------
