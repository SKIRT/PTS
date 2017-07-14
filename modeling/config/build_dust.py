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
from pts.modeling.maps.component import get_dust_map_names
from pts.core.units.parsing import parse_quantity
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

default_dust_mass = parse_quantity("1.5e7 Msun")

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The name
definition.add_required("name", "string", "name of the model")

# Dust map
definition.add_required("dust_map", "string", "choice of dust map", choices=get_dust_map_names(modeling_path))

# Flags
definition.add_flag("disk", "add dust disk", True)
definition.add_flag("additional", "add additional dust component(s)", True)

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# Dust mix settings
definition.add_optional("default_hydrocarbon_pops", "positive_integer",  "default number of hydrocarbon populations", 25)
definition.add_optional("default_enstatite_pops", "positive_integer", "default number of enstatite populations", 25)
definition.add_optional("forsterite_pops", "positive_integer", "default number of forsterite populations", 25)

## OTHER
definition.add_optional("default_dust_mass", "quantity", "default value for the dust disk mass", default_dust_mass)

# Scaleheight
definition.add_optional("dust_scaleheight_ratio", "real", "ratio of the dust scaleheight to the old stellar scaleheight", 0.5)

# -----------------------------------------------------------------
