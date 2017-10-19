#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# The significance level
definition.add_optional("i1_significance", "real", "the significance level of the IRAC I1 image below which to cut-off the stellar map", 3.0)

# Remove holes from the cutoff mask
definition.add_flag("remove_holes", "remove holes from the total cutoff mask")

# Flags
definition.add_flag("write", "write out the maps", True)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------
