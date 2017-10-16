#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.check_simulation_data Check whether simulation data is valid.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.data import SimulationData

# -----------------------------------------------------------------

# Parse arguments
definition = ConfigurationDefinition(write_config=False)
definition.add_positional_optional("prefix", "string", "simulation prefix")
config = parse_arguments("check_simulation_data", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

# Get the simulation data
data = SimulationData.from_cwd(prefix=config.prefix)

# -----------------------------------------------------------------

# Show
data.show()

# -----------------------------------------------------------------
