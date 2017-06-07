#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.config.generate_genetic import default_scale, scales

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
definition.add_required("name", "string", "fitting run name")

# Scale
definition.add_optional("default_scale", "string", "default parameter scale (also for generating the initial generation manually)", default=default_scale, choices=scales)

# -----------------------------------------------------------------
