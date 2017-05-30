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

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

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

# DUST MASS
#scale_height = 200. * u("pc")  # M51
#dust_mass = 1.5e7 * u("Msun")
#definition.add_optional("default_dust_mass", )

# -----------------------------------------------------------------
