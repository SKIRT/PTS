#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Number of processes to use
definition.add_optional("nprocesses", "positive_integer", "number of parallel processes", 8)

# Flags
definition.add_flag("write", "write results", True)

# -----------------------------------------------------------------
