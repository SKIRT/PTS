#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# Add settings

# Extraction
definition.add_flag("extract_progress", "extract progress information", False)
definition.add_flag("extract_timeline", "extract timeline information", False)
definition.add_flag("extract_memory", "extract memory information", False)

# Plotting
default_format = "pdf"
formats = ["pdf", "png"]
definition.add_flag("plot_progress", "plot progress information", False)
definition.add_flag("plot_timeline", "plot simulation timeline", False)
definition.add_flag("plot_memory", "plot memory information", False)
definition.add_flag("plot_seds", "plot the SEDs of individual simulations", True)
definition.add_optional("plotting_format", "string", "plotting format", default_format, choices=formats)

# -----------------------------------------------------------------
