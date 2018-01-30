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

# The galaxy name
definition.add_required("galaxy_name", "string", "galaxy name")

# Flags
definition.add_flag("write", "write the components and properties", False)
definition.add_flag("show", "show the parameters", True)

# The output directory
definition.add_optional("output", "directory_path", "output directory", letter="o")

# -----------------------------------------------------------------
