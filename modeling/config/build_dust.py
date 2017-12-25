#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.units.parsing import parse_quantity
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.maps.selection import ComponentMapsSelection
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
selection = ComponentMapsSelection.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

default_dust_mass = parse_quantity("1.5e7 Msun")

# -----------------------------------------------------------------

definition = definition.copy()

# The name
definition.add_required("name", "string", "name of the model")

# Dust map
definition.add_required("dust_map", "string", "choice of dust map", choices=selection.dust_map_names)

# Flags
definition.add_flag("disk", "add dust disk", True)
definition.add_flag("additional", "add additional dust component(s)", False)

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# Dust mix settings
definition.add_optional("default_hydrocarbon_pops", "positive_integer",  "default number of hydrocarbon populations", 25)
definition.add_optional("default_silicate_pops", "positive_integer", "default number of silicate populations", 25)

## OTHER
definition.add_optional("default_dust_mass", "quantity", "default value for the dust disk mass", default_dust_mass)

# Scaleheight
definition.add_optional("dust_scaleheight_ratio", "real", "ratio of the dust scaleheight to the old stellar scaleheight", 0.5)

# -----------------------------------------------------------------

# ADVANCED: FLAG TO USE DEFAULTS
definition.add_flag("use_defaults", "use the defaults (do not prompt for component parameters)", False)

# -----------------------------------------------------------------
