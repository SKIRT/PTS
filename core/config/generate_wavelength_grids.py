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
definition.add_optional("min_wavelengths_in_filter", "positive_integer", "minimum number of wavelength points to sample a filter", 5)
definition.add_optional("min_wavelengths_in_fwhm", "positive_integer", "minimum number of wavelength points to sample within inner range of filter", 3)

# Add flags
definition.add_flag("show", "show", False)
definition.add_flag("plot", "plot", False)
definition.add_flag("write", "write", False)

# -----------------------------------------------------------------
