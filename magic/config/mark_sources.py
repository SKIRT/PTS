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

# Optional
definition.add_optional("default_fwhm", "quantity", "default FWHM", "2.0 arcsec", convert_default=True)
definition.add_optional("sigma_level", "real", "sigma level", 4.0)

# -----------------------------------------------------------------
