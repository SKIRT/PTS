#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

default_min_points_per_filter = 8
default_min_points_per_fwhm = 5

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Filters and instruments for which to create the images
definition.add_positional_optional("instruments", "string_list", "instruments for which to create the images")
definition.add_positional_optional("filters", "filter_list", "filters for which to create the images")

# -----------------------------------------------------------------

# Add optional
definition.add_optional("output", "string", "output directory")

# Intermediate results
definition.add_flag("write_intermediate", "write intermediate results", False)
definition.add_flag("keep_intermediate", "keep intermediate results", False)

# Save kernels
definition.add_flag("write_kernels", "write the prepared convolution kernels", False)

# -----------------------------------------------------------------

# Add flags
definition.add_flag("spectral_convolution", "convolve over the wavelengths to get the most accurate images", True)
definition.add_flag("group", "group the images per instrument", False)

# Number of parallel processes to use to create the images
definition.add_optional("nprocesses_local", "positive_integer", "number of parallel processes to use for local calculation", 2)
definition.add_optional("nprocesses_remote", "positive_integer", "number of parallel processes to use for remote calculation", 8)

# -----------------------------------------------------------------

# LIMIT TO THE RATIO BETWEEN CONVOLUTION KERNEL FWHM AND IMAGE PIXELSCALE
definition.add_optional("max_fwhm_pixelscale_ratio", "positive_real", "maximum ratio between the convolution kernel FWHM and the image pixelscale before downsampling is applied", 50)

# REBIN TO SMALLER PIXELSCALES?
definition.add_flag("upsample", "rebin images where upsampling is required?", False)

# -----------------------------------------------------------------

# ADVANCED: regenerate?
definition.add_flag("regenerate", "regenerate images that are already present", False)

# Update PTS on the remote
definition.add_flag("deploy_pts", "deply (install or update) PTS on the remote host", True)
definition.add_flag("update_dependencies", "update PTS dependencies (use with care!)", False)
definition.add_flag("deploy_clean", "perform a clean install", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------

# Plot?
definition.add_flag("plot", "make plots", True)
definition.add_flag("plot_images", "plot the mock obseved images", False)

# -----------------------------------------------------------------

# Special options
definition.add_flag("check_wavelengths", "check the sampling of the wavelength grid", True)
definition.add_flag("ignore_bad", "ignore bad sampling of wavelength grid (just give warning)", False)
definition.add_flag("skip_ignored_bad_convolution", "skip filters that are ignored because of bad sampling (for convolution)", True)
definition.add_flag("skip_ignored_bad_closest", "skip filters that are ignored because the closest wavelength is outside of the inner region of the filter wavelength range", True)
definition.add_optional("min_npoints", "positive_integer", "minimum number of points required in filter wavelength range", default_min_points_per_filter)
definition.add_optional("min_npoints_fwhm", "positive_integer", "minimum number of points required in FWHM filter wavelength range", default_min_points_per_fwhm)

# -----------------------------------------------------------------

# Convolve
definition.add_flag("convolve", "apply spatial convolution by automatically creating PSFs for each filter")

# -----------------------------------------------------------------

# Artifical sky and stars
definition.add_flag("sky", "add artificial sky", False)
definition.add_flag("stars", "add artificial stars", False)

# -----------------------------------------------------------------

# The target unit of the observed
definition.add_optional("unit", "photometric_unit", "target unit of the images")

# -----------------------------------------------------------------
