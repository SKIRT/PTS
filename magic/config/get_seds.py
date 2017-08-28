#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.services.seds import catalog_names

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
definition.add_required("galaxy_name", "string", "galaxy name")

# The filter
definition.add_optional("catalogs", "string_list", "names of the catalogs to query", default=catalog_names, choices=catalog_names)

# Flags
definition.add_flag("list", "list the results", True)
definition.add_flag("write", "write the results", True)

# -----------------------------------------------------------------
