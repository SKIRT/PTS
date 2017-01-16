#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.dustpedia.core.sdss import sdss_bands

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")
definition.add_optional("bands", "string_list", "the bands (u/g/r/i/z)", choices=sdss_bands, default=sdss_bands)

# Output
definition.add_optional("output", "string", "the name of the output directory", default="out")

# -----------------------------------------------------------------
