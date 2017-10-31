#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.old import methods

# -----------------------------------------------------------------

# Methods
definition.add_positional_optional("methods", "string_list", "dust map making methods", default=methods, choices=methods)

# -----------------------------------------------------------------

# Flags
definition.add_flag("write", "write out the maps", True)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# CLEAR
definition.add_flag("clear", "clear already existing maps (for the methods selected)", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------
