#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_analysis_options Show information about a certain remote simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.simulation.remote import get_simulation_for_host

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_required("id", "integer", "ID of the simulation")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("show_analysis_options", definition, description="Show information about a certain remote simulation")

# -----------------------------------------------------------------

# Load the simulation
simulation = get_simulation_for_host(config.remote, config.id)

# -----------------------------------------------------------------

print(simulation.analysis)

# -----------------------------------------------------------------
