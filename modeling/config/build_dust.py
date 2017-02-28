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
definition.add_required("dust", "string", "choice of dust map", choices=get_dust_map_names(modeling_path))

# Flags
definition.add_flag("additional", "add additional dust component(s)", True)

# -----------------------------------------------------------------
