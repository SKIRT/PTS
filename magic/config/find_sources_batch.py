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
definition.add_required("dataset", "file_path", "name of the dataset file")

# Flags
definition.add_flag("find_stars", "find stars in the images")

# Flags
definition.add_flag("find_other_sources", "find other contaminating sources")

# Optional settings
definition.add_optional("galactic_catalog_file", "file_path", "galactic catalog file")
definition.add_optional("stellar_catalog_file", "file_path", "stellar catalog file")

# Regions
definition.add_optional("special_region", "file_path", "region indicating areas that require special attention")
definition.add_optional("ignore_region", "file_path", "region indicating areas that should be ignored")

# Output
definition.add_optional("output", "directory_path", "output directory", letter="o")

# -----------------------------------------------------------------
