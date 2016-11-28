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
definition = ConfigurationDefinition()

# Image path
definition.add_optional("image", "file_path", "name/path of the input image")

# The path of a region file for the sky estimation
definition.add_optional("sky_region", "file_path", "region file for the sky estimation")

# Perform sigma-clipping step
definition.add_flag("sigma_clip_mask", "sigma-clippin", True)

# Estimate the sky (obviously)
definition.add_flag("estimate", "estimate the sky", True)

# Set zero outside of the principal galaxy
definition.add_flag("set_zero_outside", "set zero outside of principal galaxy", False)

# Eliminate negative values (replace them by zero)
definition.add_flag("eliminate_negatives", "replace negative pixels by zero", False)

# Creation of the sky mask
definition.add_section("mask", "creation of sky mask")
definition.sections["mask"].add_optional("annulus_inner_factor", "real", "sky annulus inner factor (based on principal galaxy ellipse", 1.6)
definition.sections["mask"].add_optional("annulus_outer_factor", "real", "sky annulus outer factor (based on principal galaxy ellipse", 4.0)

# Sigma clipping
definition.add_section("sigma_clipping", "sigma clipping")
definition.sections["sigma_clipping"].add_optional("sigma_level", "real", "sigma level", 3.0)

# Histogram
definition.add_section("histogram", "histogram")
definition.sections["histogram"].add_flag("log_scale", "log scale", True)

# Estimation
definition.add_section("estimation", "sky estimation")
definition.sections["estimation"].add_optional("method", "string", "method used for sky estimation", "pts")
definition.sections["estimation"].add_optional("finishing_step", "string", "finishing step", choices=["polynomial", "interpolation"])

# Setting zero outside
definition.add_section("zero_outside", "setting zero outside")
definition.sections["zero_outside"].add_optional("factor", "real", "factor", 2.0)

# -----------------------------------------------------------------
