#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.dust import methods, default_methods

# -----------------------------------------------------------------

# Methods
definition.add_positional_optional("methods", "string_list", "dust map making methods", default=default_methods, choices=methods)

# -----------------------------------------------------------------

# Old stellar contribution subtraction factor
definition.add_optional("hot_factor_range", "real_range", "range of factor to create the hot dust maps", "0.2>0.7", convert_default=True)
definition.add_optional("factor_nvalues", "positive_integer", "number of factors", 8)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------

# Old stars component
old_components = ["bulge", "disk", "total"]
default_old_component = "total"
definition.add_optional("old_component", "string", "old stellar component to use to subtract diffuse emission by evolved stars", default_old_component, choices=old_components)

# -----------------------------------------------------------------
