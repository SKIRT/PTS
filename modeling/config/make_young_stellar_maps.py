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

# THIS STUFF IS NOT USED ANYMORE?
# # Remove holes from the cutoff mask
# definition.add_flag("remove_holes", "remove holes from the total cutoff mask", True)
definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor that denotes the contribution of the old stellar population to the FUV emission", "0.1,0.4", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 4)
# #definition.add_optional("best_factor", "real", "the best estimate for the value of the factor", 0.1) # WAS 0.2, then 0.15. 0.1 SEEMS BEST WHEN LOOKING AT HISTOGRAMS!
# definition.add_optional("histograms_annulus_range", "real_range", "range (min,max) of the radius (relative to the scalelength) of the area for to make histograms of the pixel values of the corrected FUV maps", "0.065,0.28", convert_default=True)
# definition.add_optional("histograms_nbins", "integer", "the number of bins in the histogram plots", 20)

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# Replot?
definition.add_flag("replot", "replot already existing plots", False)

# -----------------------------------------------------------------

old_filters = collection.old_stellar_disk_filters
if len(old_filters) == 0: raise ValueError("There are no old stellar disk maps")
elif len(old_filters) == 1: definition.add_fixed("old", "filter for the old stellar disk map to use", old_filters[0])
else: definition.add_required("old", "filter", "filter for the old stellar disk map to use", choices=old_filters[0])

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------
