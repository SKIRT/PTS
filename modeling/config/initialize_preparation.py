#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.magic.config.find_sources import definition as sources_definition
from pts.modeling.component.component import get_cache_host_id
from pts.modeling.core.environment import verify_modeling_cwd
from pts.magic.tools.catalogs import stellar_catalog_descriptions
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

definition = definition.copy()

# Add optional arguments
definition.add_positional_optional("image", "string", "the name of the image for which to run the initialization")

# Add flags
definition.add_flag("visualise", "make visualisations")

# Remote source detection
definition.add_optional("remote", "string", "remote host on which to run the source finder", choices=find_host_ids())
definition.add_flag("attached", "run remotely in attached mode")

# Add section for the source finder
definition.import_section("sources", "options for the source finder", sources_definition)

# Flags
default_catalogs = ["2MASS"]
definition.add_optional("catalogs", "string_list", "catalogs for point sources", default_catalogs, choices=stellar_catalog_descriptions)
definition.add_flag("catalog_overlapping", "only fetch catalog data in the area where all images are overlapping", True)

definition.add_flag("manual", "don't find sources, but mark them from the catalog, and let the selection be done manually", False)

# Cache
cache_host_id = get_cache_host_id(modeling_path)
if cache_host_id is not None: definition.add_flag("cache", "cache image data for which the initialized image has been created", False)
else: definition.add_fixed("cache", "caching not possible since cache host id not defined", False)

# -----------------------------------------------------------------
