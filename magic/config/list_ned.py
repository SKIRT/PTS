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
definition.add_required("galaxy", "string", "galaxy name")

# The filter
definition.add_positional_optional("filter", "filter", "the filter for which to look for images")

# Add flags
definition.add_flag("unknown", "list images of unknown filters", False)
definition.add_flag("show", "show the list", True)

# -----------------------------------------------------------------
