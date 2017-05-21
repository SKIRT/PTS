#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.core.tools.parallelization import ncores
from pts.modeling.preparation.preparer import steps
from pts.modeling.component.component import get_cache_host_id
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add required arguments
definition.add_required("image", "string", "name of the image for which to run the preparation")

# Add optional arguments
definition.add_optional("exclude_filters", "string_list", "exclude the data for these filters from the procedure that brings all data to the same resolution and pixelscale")
definition.add_flag("steps", "write the results of intermediate steps")
definition.add_flag("visualise", "make visualisations")

# Remote preparation
definition.add_optional("remote", "string", "remote host on which to run the preparation", choices=find_host_ids())
definition.add_flag("attached", "run remotely in attached mode")

# Parallelization
definition.add_optional("nprocesses", "positive_integer", "number of parallel processes for parallel computations", max(8, ncores()))

# Advanced
definition.add_flag("dustpedia_aperture", "use the DustPedia aperture instead of the galaxy region for the sky subtraction", True)
definition.add_optional("aperture_galaxy_region_factor", "real", "if the galaxy region is used for the aperture, use this scaling factor on the region", 1.5)

# Sky subtraction
definition.add_optional("annulus_inner_factor", "real", "annulus_inner_factor", 1.2)
definition.add_optional("annulus_outer_factor", "real", "annulus_outer_factor", 4.0)

definition.add_optional("saturation_expansion_factor", "real", "saturation expansion factor", 1.5)
definition.add_optional("stars_expansion_factor", "real", "stars expansion factor", 2.)

definition.add_optional("rerun", "string", "rerun a certain step for all images", choices=steps)

# Cache
cache_host_id = get_cache_host_id(modeling_path)
if cache_host_id is not None: definition.add_flag("cache", "cache intermediate image files to the remote host storage", False)
else: definition.add_fixed("cache", "caching not possible since cache host id not defined", False)

# -----------------------------------------------------------------
