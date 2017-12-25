#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

formats = ["pdf", "png"]
default_format = "pdf"

# -----------------------------------------------------------------

types = ["old", "young", "ionizing", "dust"]
#features = ["maps", "contours", "profiles", "negatives"]
features = ["maps", "masks", "deprojected", "deprojected_skirt", "edgeon"]

# -----------------------------------------------------------------

definition = definition.copy()

# Types
definition.add_positional_optional("types", "string_list", "types of maps to be plotted", default=types, choices=types)
definition.add_positional_optional("features", "string_list", "features of the maps to plot", default=features, choices=features)

# -----------------------------------------------------------------

definition.add_optional("format", "string", "plotting format", default=default_format, choices=formats)

# -----------------------------------------------------------------

definition.add_flag("replot", "replot existing plots", False)

# -----------------------------------------------------------------
