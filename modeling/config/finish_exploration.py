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

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
definition.add_required("name", "string", "fitting run name")

# N significant digits: NO NOT ANYMORE: NOW DEFINED IN THE FITTING RUN DEFINITION
#definition.add_optional("ndigits", "positive_integer", "number of digits (in base 10) for binary genomes", 4)

# -----------------------------------------------------------------
