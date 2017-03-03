#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.maps.component import get_old_stellar_map_names, get_young_stellar_map_names, get_ionizing_stellar_map_names

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Stellar maps
definition.add_required("old_stars_map", "string", "choice of old stars map", choices=get_old_stellar_map_names(modeling_path))
definition.add_required("young_stars_map", "string", "choice of young stars map", choices=get_young_stellar_map_names(modeling_path))
definition.add_required("ionizing_stars_map", "string", "choice of ionizing stars map", choices=get_ionizing_stellar_map_names(modeling_path))

# Flags
definition.add_flag("bulge", "add bulge", True)
definition.add_flag("old", "add old stars", True)
definition.add_flag("young", "add young stars", True)
definition.add_flag("ionizing", "add ionizing stars", True)
definition.add_flag("additional", "add additional stellar component(s)", True)

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# -----------------------------------------------------------------
