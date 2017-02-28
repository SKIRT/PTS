#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.build.component import get_model_names
from pts.core.tools import filesystem as fs
from pts.modeling.maps.component import get_dust_map_names, get_old_stellar_map_names, get_young_stellar_map_names, get_ionizing_stellar_map_names

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Flags
definition.add_flag("bulge", "add bulge", True)
definition.add_flag("old", "add old stars", True)
definition.add_flag("young", "add young stars", True)
definition.add_flag("ionizing", "add ionizing stars", True)
definition.add_flag("additional", "add additional stellar component(s)", True)

# -----------------------------------------------------------------
