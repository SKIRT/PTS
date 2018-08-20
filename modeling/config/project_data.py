#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.projection.projector import orientations

# -----------------------------------------------------------------

default_spacing_measure = "mean"
spacing_measures = ["min", "max", "mean", "median"]
default_spacing_factor = 10.

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# 3D data file
definition.add_required("filename", "file_path", "path of the 3D data file")

# Orientations
definition.add_required("orientations", "string_list", "orientations from which to project the data", choices=orientations)

# -----------------------------------------------------------------

# Width/height
definition.add_optional("height", "length_quantity", "maximum height above/below midplane for the faceon projection")
definition.add_optional("width", "length_quantity", "maximum width before/behind center vertical plane for the edgeon projection")

# Output directory
definition.add_optional("output", "directory_path", "output directory")

# -----------------------------------------------------------------

# Spacing of the maps
definition.add_optional("spacing", "string", "method of determining the grid cell spacing", default_spacing_measure, choices=spacing_measures)
definition.add_optional("spacing_factor", "positive_real", "factor by which to multiply the grid cells spacing measure to become the actual map spacing", default_spacing_factor)

# -----------------------------------------------------------------

# Do interpolation on the maps
definition.add_flag("interpolate", "do interpolation", False)

# -----------------------------------------------------------------

# Do plot
definition.add_flag("plot", "do plotting", True)

# Plotting sections
definition.add_section("plotting", "plotting options")
definition.sections["plotting"].add_optional("interval", "string", "interval", "minmax")
definition.sections["plotting"].add_flag("contours", "show contours", False)
definition.sections["plotting"].add_optional("ncontours", "positive_integer", "number of contour levels", 5)
definition.sections["plotting"].add_optional("contours_color", "string", "color for the contour lines", "white")
definition.sections["plotting"].add_optional("minmax", "real_pair", "plotting minimum and maximum value")

# -----------------------------------------------------------------
