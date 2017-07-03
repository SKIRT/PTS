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
definition.add_required("remote", "string", "remote host to use for the aperture correction calculation", find_host_ids())

# Add option
definition.add_optional("noise_method", "string", "method to use for the aperture noise calculation", choices=["caapr", "pts"], default="caapr")

# -----------------------------------------------------------------
