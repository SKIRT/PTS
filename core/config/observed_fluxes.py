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

# Add optional
definition.add_optional("output", "string", "output directory")

# Add flags
definition.add_flag("spectral_convolution", "convolve over the wavelengths to get the most accurate fluxes", True)

# -----------------------------------------------------------------

definition.add_flag("from_images", "calculate observed fluxes from images created from the output datacubes", False)
definition.add_flag("write_images", "write the images created from the output datacubes", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "make plots", True)
definition.add_flag("plot_seds", "plot the mock observed SEDs", False)
definition.add_flag("plot_images", "plot the images created from the output datacubes", False)

# -----------------------------------------------------------------

definition.add_optional("images_nprocesses", "positive_integer", "number of parallel processes for creating images", 1)

# -----------------------------------------------------------------

# Update PTS on the remote
definition.add_flag("deploy_pts", "deply (install or update) PTS on the remote host", True)
definition.add_flag("update_dependencies", "update PTS dependencies (use with care!)", False)
definition.add_flag("deploy_clean", "perform a clean install", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------

# Special options
definition.add_flag("check_wavelengths", "check the sampling of the wavelength grid", True)
definition.add_flag("ignore_bad", "ignore bad sampling of wavelength grid (just give warning)", False)
definition.add_flag("skip_ignored_bad_convolution", "skip filters that are ignored because of bad sampling (for convolution)", True)
definition.add_flag("skip_ignored_bad_closest", "skip filters that are ignored because the closest wavelength is outside of the inner region of the filter wavelength range", True)
definition.add_optional("min_npoints", "positive_integer", "minimum number of points required in filter wavelength range", default_min_points_per_filter)
definition.add_optional("min_npoints_fwhm", "positive_integer", "minimum number of points required in FWHM filter wavelength range", default_min_points_per_fwhm)

# -----------------------------------------------------------------

# Unit for the fluxes
definition.add_optional("unit", "photometric_unit", "unit for the fluxes (and images)", "Jy")

# -----------------------------------------------------------------
