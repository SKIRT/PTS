#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.core.config.observed_fluxes import default_min_points_per_filter, default_min_points_per_fwhm

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_positional_optional("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "integer", "ID of the simulation")

# Add flags
definition.add_flag("ignore_missing_data", "ignore missing data when analysing the simulations", False)

# -----------------------------------------------------------------

# Special options for misc/flux calculation
definition.add_flag("check_wavelengths", "check the sampling of the wavelength grid", True)
definition.add_flag("ignore_bad", "ignore bad sampling of wavelength grid (just give warning)", False)
definition.add_flag("skip_ignored_bad_convolution", "skip filters that are ignored because of bad sampling (for convolution)", True)
definition.add_flag("skip_ignored_bad_closest", "skip filters that are ignored because the closest wavelength is outside of the inner region of the filter wavelength range", True)
definition.add_optional("min_npoints", "positive_integer", "minimum number of points required in filter wavelength range", default_min_points_per_filter)
definition.add_optional("min_npoints_fwhm", "positive_integer", "minimum number of points required in FWHM filter wavelength range", default_min_points_per_fwhm)

# -----------------------------------------------------------------

# Update PTS on the remote
definition.add_flag("deploy_pts", "deply (install or update) PTS on the remote host", True)
definition.add_flag("update_dependencies", "update PTS dependencies (use with care!)", False)
definition.add_optional("pubkey_password", "string", "pubkey password for accessing the repo URL")

# -----------------------------------------------------------------

# Local: override analysis options and perform everything locally
definition.add_flag("local", "override analysis options and perform everything locally", False)

# -----------------------------------------------------------------

# Flags
definition.add_flag("extract", "do extraction", True)
definition.add_flag("plot", "do plotting", True)
definition.add_flag("misc", "do misc", True)

# -----------------------------------------------------------------
