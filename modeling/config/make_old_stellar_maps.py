#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# Remove holes from the cutoff mask
definition.add_flag("remove_holes", "remove holes from the total cutoff mask")

# Flags
definition.add_flag("write", "write out the maps", True)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------
