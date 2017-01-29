#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.fitting.component import get_generation_names
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
definition.add_positional_optional("features", "string_list", "features to be plotted")
definition.add_positional_optional("generation", "string", "generation for which to plot the features", choices=get_generation_names(fs.cwd()))

# -----------------------------------------------------------------
