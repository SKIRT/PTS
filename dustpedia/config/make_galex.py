#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.dustpedia.core.galex import galex_bands

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")
definition.add_optional("bands", "string_list", "the bands (NUV/FUV)", choices=galex_bands, default=galex_bands)

# Output
definition.add_optional("output", "string", "the name of the output directory", default="out")

# Advanced
definition.add_optional("max_nobservations_fuv", "positive_integer", "limit the number of FUV observations")
definition.add_optional("max_nobservations_nuv", "positive_integer", "limit the number of NUV observations")

# -----------------------------------------------------------------
