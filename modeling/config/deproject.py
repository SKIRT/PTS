#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

methods = ["skirt", "pts"]
default_method = "pts"
default_downsample_factor = 2.
default_parallelization = "2:1:2" # 2 cores, 1 process, 2 threads per core

# -----------------------------------------------------------------

definition = definition.copy()

# Method
definition.add_optional("method", "string", "method for deprojection", default_method, choices=methods)

# Downsample
definition.add_optional("downsample_factor", "positive_real", "downsample factor of the deprojected map w.r.t the original map", default_downsample_factor)

# Writing
definition.add_section("writing", "writing options")
definition.sections["writing"].add_flag("deprojections", "write the deprojections", True)
definition.sections["writing"].add_flag("maps", "write the maps (or don't clear them)", True)

# -----------------------------------------------------------------

# DUST GRIDS
# Settings for the dust grid
definition.add_section("dg", "options for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("rel_scale", "real", "the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 10.) #1.)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "the maximum mass fraction per cell", 1e-5) #1e-6)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 10.)
definition.sections["dg"].add_optional("bintree_min_level", "integer", "minimum depth level for binary trees", 9)
definition.sections["dg"].add_optional("octtree_min_level", "integer", "minimum depth level for octrees", 3)

# -----------------------------------------------------------------

# Simulation options
definition.add_optional("parallelization", "parallelization", "parallelization scheme for SKIRT", default_parallelization, convert_default=True)

# -----------------------------------------------------------------
