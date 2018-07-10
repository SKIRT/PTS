#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create
definition = ConfigurationDefinition()

# Commands to be run
definition.add_optional("commands", "string_list", "commands to be run in interactive mode")

# Interactive mode
definition.add_flag("interactive", "use interactive mode", default=None)

# -----------------------------------------------------------------

# Showing
definition.add_flag("show", "show stuff", True)
definition.add_flag("show_components", "show the model components", False)

# Plotting
definition.add_flag("plot", "do plotting", True)

# Writing
definition.add_flag("write", "do writing", True)

# -----------------------------------------------------------------
