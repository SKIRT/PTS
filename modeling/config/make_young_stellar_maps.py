#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.collection import MapsCollection
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()

# Create the maps collection
collection = MapsCollection.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

# The significance level
definition.add_optional("fuv_significance", "real", "the significance level of the FUV image below which to cut-off the stellar map", 3.0)

# Old stellar contribution subtraction factors
definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor that denotes the contribution of the old stellar population to the FUV emission", "0.03,0.05", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 4)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# Clear?
definition.add_flag("clear", "clear all maps", False)

# -----------------------------------------------------------------

# Old component
old_components = ["bulge", "disk", "total"]
default_old_component = "total"
definition.add_optional("old_component", "string", "old stellar component to use to subtract diffuse emission by evolved stars", default_old_component, choices=old_components)

# -----------------------------------------------------------------

old_filters = collection.old_stellar_disk_filters
if len(old_filters) == 0: raise ValueError("There are no old stellar disk maps")
elif len(old_filters) == 1: definition.add_fixed("old", "filter for the old stellar disk map to use", old_filters[0])
else: definition.add_required("old", "filter", "filter for the old stellar disk map to use", choices=old_filters)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------

definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# -----------------------------------------------------------------

definition.add_flag("use_cortese", "use Cortese FUV attenuation maps", True)
definition.add_flag("use_buat", "use Buat FUV attenuation maps", True)

# -----------------------------------------------------------------

# Select specific input maps
definition.add_optional("attenuation_maps", "string_list", "names of FUV attenuation maps to use")
definition.add_flag("select_attenuation", "select specific FUV attenuation maps", False)

# -----------------------------------------------------------------
