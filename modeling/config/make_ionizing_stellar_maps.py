#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# Subtraction factor
definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor that denotes the contribution of the old stellar population to the MIPS 24 micron emission", "0.2,0.7", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 8)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# Clear?
definition.add_flag("clear", "clear maps", False)

# -----------------------------------------------------------------
