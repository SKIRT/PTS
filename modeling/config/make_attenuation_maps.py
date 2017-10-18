#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.attenuation import methods

# -----------------------------------------------------------------

# Add optional
#definition.add_optional("ssfr_colour", "string", "SSFR colour to use", default="FUV-H", choices=["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"])

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# -----------------------------------------------------------------

# Select specific input maps
definition.add_flag("select_ssfr", "select specific sSFR maps", False)
definition.add_flag("select_tir", "select specific TIR maps", False)

# -----------------------------------------------------------------

definition.add_flag("plot", "plot Cortese sSFR pixel masks", False)

# -----------------------------------------------------------------

# Methods
definition.add_positional_optional("methods", "string_list", "attenuation map making methods", default=methods, choices=methods)

# -----------------------------------------------------------------
