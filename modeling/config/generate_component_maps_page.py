#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition
from pts.modeling.maps.component import default_central_ellipse_factor

# -----------------------------------------------------------------

# Show
definition.add_flag("show", "show the page", False)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# Flags
definition.add_flag("thumbnails", "add map thumbnails", True)
definition.add_optional("thumbnail_height", "positive_integer", "height of the thumbnails (in pixels)", 50)
definition.add_flag("previews", "add previews of the maps when hovering over the thumbnails", True)
definition.add_flag("methods", "make a separate table for each method", True)

# Filtering maps
definition.add_flag("filter", "filter the maps", True)
definition.add_optional("central_ellipse_factor", "real", "central ellipse factor for filtering", default_central_ellipse_factor)
definition.add_optional("ninvalid_pixels_tolerance", "percentage", "number of invalid pixels tolerated within central ellipse", "10", convert_default=True)
definition.add_optional("nzero_pixels_tolerance", "percentage", "number of zero pixels tolerated within central ellipse", "40", convert_default=True)
definition.add_optional("nnegative_pixels_tolerance", "percentage", "number of negative pixels tolerated within central ellipse", "20", convert_default=True)

definition.add_flag("hide_invalid", "hide invalid maps", False)
definition.add_flag("hide_zero", "hide maps with too many zeros", False)
definition.add_flag("hide_constant", "hide constant maps", False)
definition.add_flag("hide_negative", "hide maps with too many negative values", False)

# Group maps with different factors but same origins
definition.add_flag("group_factors", "group maps that only differ by a factor in their name", True)

# Resolution
definition.add_optional("max_pixelscale", "angle", "maximum pixelscale")
definition.add_optional("max_fwhm", "angle", "maximum FWHM")

# -----------------------------------------------------------------
