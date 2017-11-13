#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add option for the remote host ID
definition.add_positional_optional("remote", "string", "remote host to use for the aperture correction calculation", choices=find_host_ids(schedulers=False))

# Add option
noise_methods = ["caapr", "pts"]
definition.add_optional("noise_method", "string", "method to use for the aperture noise calculation", choices=noise_methods, default="caapr")

# -----------------------------------------------------------------

definition.add_flag("plot", "plot SEDs", True)

# -----------------------------------------------------------------

definition.add_optional("reprocess", "lazy_broad_band_filter_list", "reprocess images and calculate fluxes for certain filters")

# -----------------------------------------------------------------

definition.add_optional("plot_images", "lazy_broad_band_filter_list", "plot the image frame and masks for certain filters")

# -----------------------------------------------------------------
