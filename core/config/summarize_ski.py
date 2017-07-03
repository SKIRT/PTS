#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("ski", "file_path", "the ski file")

# Add flags
definition.add_flag("show", "show the results", True)

# -----------------------------------------------------------------
