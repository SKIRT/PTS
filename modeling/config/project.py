#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

default_downsample_factor = 2.

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Downsample
definition.add_optional("downsample_factor", "positive_real", "downsample factor of the deprojected map w.r.t the original map", default_downsample_factor)

# Writing
definition.add_section("writing", "writing options")

# -----------------------------------------------------------------

# Faceon, edgeon?
definition.add_flag("faceon", "create faceon projections", True)
definition.add_flag("edgeon", "create edgeon projections", True)

# -----------------------------------------------------------------

# Properties of the projection
definition.add_required("distance", "length_quantity", "distance to the object")
definition.add_required("center", "pixelcoordinate", "pixel coordinate of the center of the object")
definition.add_optional("npixels", "pixel_shape", "number of pixels")
definition.add_optional("field", "physical_extent", "field of view")

# Further galaxy projection properties
definition.add_optional("inclination", "angle", "inclination angle")
definition.add_optional("azimuth", "angle", "azimuth angle", "0 deg", convert_default=True)
definition.add_optional("position_angle", "angle", "position angle")

# -----------------------------------------------------------------
