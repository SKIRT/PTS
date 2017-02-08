#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Flags
definition.add_flag("recursive", "look for images in directories recursively", True)
definition.add_flag("list", "list the found images", True)

definition.add_flag("pixelscale", "arrange by increasing pixelscale")
definition.add_flag("fwhm", "arrange by increasing FWHM")

# -----------------------------------------------------------------
