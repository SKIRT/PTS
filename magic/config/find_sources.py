#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from .find_point import definition as point_definition
from .find_extended import definition as extended_definition
from .find_other import definition as other_definition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The dataset or image
definition.add_required("dataset", "file_path", "name of the dataset file or image file")

# Number of parallel processes
definition.add_optional("nprocesses", "integer", "number of parallel processes for the preparation", 8)

# Flags
definition.add_flag("find_stars", "find stars in the images", True)

# Flags
definition.add_flag("find_other_sources", "find other contaminating sources", True)

# Optional settings
definition.add_optional("extended_sources_catalog", "file_path", "catalog file for extended sources")
definition.add_optional("point_sources_catalog", "file_path", "catalog file for point sources")

# Regions
definition.add_optional("special_region", "file_path", "region indicating areas that require special attention")
definition.add_optional("ignore_region", "file_path", "region indicating areas that should be ignored")

# Output
definition.add_optional("output", "directory_path", "output directory", letter="o")
definition.add_optional("input", "directory_path", "input directory", letter="i")

# Sections
definition.import_section("extended", "options for extended source finder", extended_definition)
definition.import_section("point", "options for point source finder", point_definition)
definition.import_section("other", "options for finding other contaminating sources", other_definition)

# -----------------------------------------------------------------
