#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

default_npackages = 1e7
default_parallelization = "2:1:2" # 2 cores, 1 process, 2 threads per core

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# SKIRT options
definition.add_optional("npackages", "positive_integer", "number of photon packages", default_npackages)
definition.add_optional("parallelization", "parallelization", "parallelization scheme for the simulations", default_parallelization, convert_default=True)

# -----------------------------------------------------------------
