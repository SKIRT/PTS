#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.decomposition.decomposition import scalelength_scaleheight_ratios, degeyter_ratio
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load the modeling environment
environment = load_modeling_environment_cwd()

# Set the default filter
default_filter_name = "IRAC I1" if "IRAC I1" in environment.nir_filter_names else None

# -----------------------------------------------------------------

definition = definition.copy()

# The filter for decomposition: choose from NIR filters
definition.add_optional("filter", "filter", "filter for which to use the data for decomposition", default=default_filter_name, convert_default=True, choices=environment.nir_filters)

# The method
definition.add_optional("method", "string", "method for decomposition", default="s4g", choices=["s4g", "fit", "imfit"])

# Add optional
choices = dict()
choices["azimuth"] = "azimuth angle and y flattening"
choices["tilt"] = "tilt angle and z flattening"
definition.add_optional("bulge_deprojection_method", "string", "method of deprojecting a 2D bulge with position angle difference w.r.t. the disk", choices=choices, default="tilt")

# -----------------------------------------------------------------

definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=degeyter_ratio, suggestions=scalelength_scaleheight_ratios)

# -----------------------------------------------------------------
