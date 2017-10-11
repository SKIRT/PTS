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

# Add optional
definition.add_optional("output", "string", "output directory")

# -----------------------------------------------------------------

# Define the wavelengths for the RGB
#wavelengths = dict()
#wavelengths["optical"] = (0.77, 0.55, 0.33)
#wavelengths["infrared"] = (333, 100, 24)

# -----------------------------------------------------------------
