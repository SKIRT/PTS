#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# -----------------------------------------------------------------

# Flags
definition.add_flag("use_session", "use remote python session to create the images", False)
definition.add_flag("show", "show the page", False)

# -----------------------------------------------------------------
