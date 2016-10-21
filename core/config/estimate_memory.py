#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("ski", "file_path", "path to the ski file")

# Add optional
definition.add_optional("input", "directory_path", "path to the input directory")

# Add optional
definition.add_optional("ncells", "integer", "number of dust cells (useful only when the ski file includes a tree dust grid)")

# -----------------------------------------------------------------
