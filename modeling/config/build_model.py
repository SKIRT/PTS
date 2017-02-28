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

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
model_names = get_model_names(fs.cwd())
if len(model_names) == 0: definition.add_positional_optional("name", "string", "name for the model", default="standard")
else: definition.add_required("name", "string", "name for the model")

# Stellar and dust maps
definition.add_required("old_stars", "string", "choice of old stars map", choices=get_old_stellar_map_names(fs.cwd()))
definition.add_required("young_stars", "string", "choice of young stars map", choices=get_young_stellar_map_names(fs.cwd()))
definition.add_required("ionizing_stars", "string", "choice of ionizing stars map", choices=get_ionizing_stellar_map_names(fs.cwd()))
definition.add_required("dust", "string", "choice of dust map", choices=get_dust_map_names(fs.cwd()))

# -----------------------------------------------------------------
