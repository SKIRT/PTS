#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

default_downsample_factor = 2.
default_npackages = 1e7
default_parallelization = "2:1:2" # 2 cores, 1 process, 2 threads per core

# -----------------------------------------------------------------

definition = definition.copy()

# Writing
definition.add_section("writing", "writing options")
definition.sections["writing"].add_flag("projections", "write the projections", True)

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

# SKIRT options
definition.add_optional("npackages", "positive_integer", "number of photon packages", default_npackages)
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations", default_parallelization, convert_default=True)

# -----------------------------------------------------------------
