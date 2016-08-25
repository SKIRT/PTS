#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The significance level
definition.add_optional("mips24_significance", "real", "significance level of the MIPS 24 micron image below which to cut-off the stellar map", 3.0)
definition.add_optional("halpha_significance", "real", "significance level of the H-alpha image below which to cut-off the stellar map", 3.0)

# Remove holes from the cutoff mask
definition.add_flag("remove_holes", "remove holes from the total cutoff mask")

# Map creation parameters
definition.add_optional("factor_range", "real_range", "range (min,max) of values for the factor that denotes the contribution of the old stellar population to the MIPS 24 micron emission", "0.2,0.7", convert_default=True)
definition.add_optional("factor_nvalues", "integer", "the number of values for the factor", 8)
definition.add_optional("best_factor", "real", "the best estimate for the value of the factor", 0.48)

# Histogram options
definition.add_optional("histograms_annulus_range", "real_range", "range (min,max) of the radius (relative to the scalelength) of the area for to make histograms of the pixel values of the corrected 24 micron maps", "0.065,0.28", convert_default=True)
definition.add_optional("histograms_nbins", "integer", "the number of bins in the histogram plots", 20)

# -----------------------------------------------------------------
