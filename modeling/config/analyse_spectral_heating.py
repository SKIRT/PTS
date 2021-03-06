#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import warnings

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# ANALYSIS RUNS
if runs.empty: warnings.warn("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run for which to analyse the projected heating", runs.last_name, runs.names)

# -----------------------------------------------------------------

definition.add_optional("emission_filters", "lazy_broad_band_filter_list", "filters for which to plot a map of the heating fraction by dust emission", "W3,W4,MIPS 24mu,Herschel", convert_default=True)
definition.add_optional("absorption_filters", "lazy_broad_band_filter_list", "filters for which to plot a map of the heating fraction by dust absorption", "GALEX,SDSS", convert_default=True)

# -----------------------------------------------------------------

# Convolution?
definition.add_flag("spectral_convolution", "use spectral convolution for maps from the cubes", False)

# -----------------------------------------------------------------

# For creating interpolated maps
definition.add_optional("min_ncells", "positive_integer", "minimum number of cells for each pixel of the interpolated maps of the heating fraction", 5)
definition.add_optional("not_nans_dilation_radius", "positive_real", "radius for dilating not-nans", 3)

# -----------------------------------------------------------------
