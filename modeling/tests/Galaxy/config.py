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
definition = ConfigurationDefinition(write_config=False)

definition.add_optional("nxpixels", "positive_integer", "number of x pixels for reference simulation", 1000)
definition.add_optional("nypixels", "positive_integer", "number of y pixels for reference simulation", 1000)

# -----------------------------------------------------------------
