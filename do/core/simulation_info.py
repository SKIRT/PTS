#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.simulation_info Show information about a certain remote simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_required("simulation_id", "integer", "ID of the simulation")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("simulation_info", definition, description="Show information about a certain remote simulation")

# -----------------------------------------------------------------

print("Not implemented yet")

# -----------------------------------------------------------------
