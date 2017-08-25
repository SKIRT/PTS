#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# Add optional
#definition.add_optional("ssfr_colour", "string", "SSFR colour to use", default="FUV-H", choices=["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"])

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# -----------------------------------------------------------------
