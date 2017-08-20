#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.config.find_sources import definition as sources_definition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Galaxy
definition.add_optional("galaxy_name", "string", "name of the galaxy", "M51")
definition.add_optional("filter", "filter", "image filter", "R")
definition.add_optional("year", "positive_integer", "year of the image")

# Add section for the source finder
definition.import_section("sources", "options for the source finder", sources_definition)

# Flags
default_catalogs = ["2MASS"]
definition.add_optional("catalogs", "string_list", "catalogs for point sources", default_catalogs)
definition.add_flag("catalog_overlapping", "only fetch catalog data in the area where all images are overlapping", True)

# Manual (use SourceMarker)
definition.add_flag("manual", "don't find sources, but mark them from the catalog, and let the selection be done manually", False)

# -----------------------------------------------------------------
