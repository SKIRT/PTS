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

# Add flags
definition.add_flag("spectral_convolution", "convolve over the wavelengths to get the most accurate fluxes", True)

# -----------------------------------------------------------------
