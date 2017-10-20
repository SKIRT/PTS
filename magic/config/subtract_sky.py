#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.sky.skysubtractor import estimation_methods, finishing_steps, interpolation_methods, estimators, noise_estimators

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Image path
definition.add_positional_optional("image", "file_path", "name/path of the input image")
definition.add_optional("sources_mask_plane", "string", "name of the plane in the input image which contains masks of the sources")
definition.add_optional("principal_shape_region", "file_path", "region file with one shape that indicates the contour of the principal galaxy")

# The path of a region file for the sky estimation
definition.add_optional("sky_region", "file_path", "region file for the sky estimation")

# Perform sigma-clipping step
definition.add_flag("sigma_clip_mask", "sigma-clipping", True)

# Estimate the sky (obviously)
definition.add_flag("estimate", "estimate the sky", True)

# Set zero outside of the principal galaxy
definition.add_flag("set_zero_outside", "set zero outside of principal galaxy", True)

# Eliminate negative values (replace them by zero)
definition.add_flag("eliminate_negatives", "replace negative pixels by zero", False)

# Creation of the sky mask
definition.add_section("mask", "creation of sky mask")
definition.sections["mask"].add_optional("annulus_inner_factor", "real", "sky annulus inner factor (based on principal galaxy ellipse", 1.2)
definition.sections["mask"].add_optional("annulus_outer_factor", "real", "sky annulus outer factor (based on principal galaxy ellipse", 4.0)
definition.sections["mask"].add_optional("saturation_expansion_factor", "real", "expansion factor for saturation regions", 1.5)
definition.sections["mask"].add_optional("stars_expansion_factor", "real", "expansion factor for star regions", 1.5)

# Sigma clipping
definition.add_section("sigma_clipping", "sigma clipping")
definition.sections["sigma_clipping"].add_optional("sigma_level", "positive_real", "sigma level", 3.0)
definition.sections["sigma_clipping"].add_optional("niterations", "positive_integer", "number of iterations", 5)

# Histogram
definition.add_section("histogram", "histogram")
definition.sections["histogram"].add_flag("log_scale", "log scale", True)

# Estimation
definition.add_section("estimation", "sky estimation")
definition.sections["estimation"].add_optional("method", "string", "method used for sky estimation", "photutils", choices=estimation_methods)
definition.sections["estimation"].add_optional("finishing_step", "string", "finishing step", choices=finishing_steps, default="interpolation")
definition.sections["estimation"].add_optional("interpolation_method", "string", "method of interpolation (finishing step)", choices=interpolation_methods, default="zoom")
definition.sections["estimation"].add_optional("aperture_radius", "positive_real", "aperture radius in pixel coordinates (if not defined, aperture_fwhm_factor * fwhm of the frame will be used)")
definition.sections["estimation"].add_optional("aperture_fwhm_factor", "positive_real", "aperture radius = aperture_fwhm_factor * frame FWHM", 3.0)
definition.sections["estimation"].add_optional("relative_napertures_max", "positive_real", "fraction of the theoretical maximal number of apertures to be actually used", 0.5)
definition.sections["estimation"].add_optional("min_napertures", "positive_integer", "minimum number of sky apertures", 50)
definition.sections["estimation"].add_optional("polynomial_degree", "positive_integer", "degree of the polynomial for the finishing step", 2)
definition.sections["estimation"].add_optional("estimator", "string", "estimator for the sky in each aperture", choices=estimators, default="sextractor")
definition.sections["estimation"].add_optional("noise_estimator", "string", "estimator for the noise in each aperture", choices=noise_estimators, default="stddev")
definition.sections["estimation"].add_optional("photutils_fixed_width", "positive_integer", "fixed value for the width of the grid meshes (otherwise 2 * aperture_fwhm_factor * fwhm is used or 2 * aperture_radius)", suggestions=[50])
definition.sections["estimation"].add_optional("photutils_filter_size", "positive_integer", "filter size", 3)
definition.sections["estimation"].add_optional("photutils_global", "string", "use photutils method but only use the global mean/median value and noise (if None, actual frames are used)", choices=["mean", "median"]) #suggestions=["median"])
definition.sections["estimation"].add_optional("relative_mesh_scale", "positive_real", "scale of the meshes relative to the aperture width = 2 x aperture radius", 0.5)
definition.sections["estimation"].add_flag("fill_before_mesh", "fill frame with aperture values before mesh (instead of just creating mesh from aperture values)", False)
definition.sections["estimation"].add_optional("fill_method", "string", "fill method", "NN", choices=["NN", "linear", "cubic", "quintic", "zoom"])
#definition.sections["estimation"].add_flag("do_mesh", "grid into mesh before final interpolation", True)

# Setting zero outside
definition.add_section("zero_outside", "setting zero outside")
definition.sections["zero_outside"].add_optional("factor", "real", "factor", 2.0)

# Flags
definition.add_flag("write", "writing")
definition.add_flag("plot", "plotting")

# Add extra mask to output
definition.add_flag("add_extra_mask", "add the extra mask to the output maps", True)
definition.add_optional("mask_value", "real", "value to replace masked pixels", 0.0)

# -----------------------------------------------------------------
