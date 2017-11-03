#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.dust import methods, default_methods
from pts.modeling.maps.collection import MapsCollection
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Get the modeling path
modeling_path = verify_modeling_cwd()

# Create the maps collection
collection = MapsCollection.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

# Methods
definition.add_positional_optional("methods", "string_list", "dust map making methods", default=default_methods, choices=methods)

# -----------------------------------------------------------------

# Old stellar contribution subtraction factor
definition.add_optional("hot_factor_range", "real_range", "range of factor to create the hot dust maps", "0.2>0.45", convert_default=True)
definition.add_optional("factor_nvalues", "positive_integer", "number of factors", 8)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot
definition.add_flag("replot", "replot already existing plots", False)

# CLEAR
definition.add_flag("clear", "clear already existing maps (for the methods selected)", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------

# Old stars component
old_components = ["bulge", "disk", "total"]
default_old_component = "total"
definition.add_optional("old_component", "string", "old stellar component to use to subtract diffuse emission by evolved stars", default_old_component, choices=old_components)

# -----------------------------------------------------------------

old_filters = collection.old_stellar_disk_filters
if len(old_filters) == 0: raise ValueError("There are no old stellar disk maps")
elif len(old_filters) == 1: definition.add_fixed("old", "filter for the old stellar disk map to use", old_filters[0])
else: definition.add_required("old", "filter", "filter for the old stellar disk map to use", choices=old_filters)

# -----------------------------------------------------------------
