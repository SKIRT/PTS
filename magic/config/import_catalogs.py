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

definition.add_section("galaxies", "galaxy settings")
definition.sections["galaxies"].add_flag("use_catalog_file", "use catalog file")
definition.sections["galaxies"].add_optional("catalog_path", "file_path", "catalog file")

definition.add_section("stars", "stars settings")
definition.sections["stars"].add_flag("use_catalog_file", "use catalog file")
definition.sections["stars"].add_optional("catalog_path", "file_path", "catalog file")

definition.sections["stars"].add_section("fetching", "fetching stellar catalogs")
definition.sections["stars"].sections["fetching"].add_optional("catalogs", "string_list", "catalogs to use", default=["II/246"], choices=["UCAC4", "NOMAD", "PPMXL", "II/246"])

definition.add_section("writing", "writing options")
definition.sections["writing"].add_optional("galactic_catalog_path", "string", "write galactic catalog file")
definition.sections["writing"].add_optional("stellar_catalog_path", "string", "write stellar catalog file")

# -----------------------------------------------------------------
