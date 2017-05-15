#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.magic.config.find_sources import definition as sources_definition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add optional arguments
definition.add_optional("image", "string", "the name of the image for which to run the initialization")

# Add flags
definition.add_flag("visualise", "make visualisations")

# Remote source detection
definition.add_optional("remote", "string", "remote host on which to run the source finder", choices=find_host_ids())
definition.add_flag("attached", "run remotely in attached mode")

# Add section for the source finder
definition.import_section("sources", "options for the source finder", sources_definition)

# Flags
default_catalogs = ["II/246"]
definition.add_optional("catalogs", "string_list", "catalogs for point sources", default_catalogs)
definition.add_flag("catalog_overlapping", "only fetch catalog data in the area where all images are overlapping", True)

definition.add_flag("manual", "don't find sources, but mark them from the catalog, and let the selection be done manually", False)

# -----------------------------------------------------------------
