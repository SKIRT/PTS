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

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Remake?
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------

# Select specific input maps
definition.add_flag("select_ssfr", "select specific sSFR maps", False)
definition.add_flag("select_tir", "select specific TIR maps", False)

# -----------------------------------------------------------------

# Sepcific input maps
definition.add_optional("ssfrs", "string_list", "names of the sSFR maps to use")
definition.add_optional("tirs", "string_list", "names of the TIR maps to use")

# Methods
definition.add_optional("tir_methods", "string_list", "only use TIR maps created with these methods")

# -----------------------------------------------------------------

definition.add_flag("debug_plots", "plot Cortese sSFR pixel masks", False)

# -----------------------------------------------------------------

# Methods
definition.add_positional_optional("methods", "string_list", "attenuation map making methods", default=methods, choices=methods)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------
