#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.create_assignment Create simulation assignment file from simulation names.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids
from pts.core.launch.manager import SimulationManager

# -----------------------------------------------------------------

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Simulations
definition.add_positional_optional("names", "string_list", "simulation names")
definition.add_optional("remotes", "string_list", "remote hosts for which to look for matching simulations", default=host_ids, choices=host_ids)
definition.add_flag("from_directories", "use directory names as simulation names")

# Options
definition.add_flag("local", "treat simulations without a match as local simulations", False)
definition.add_flag("success", "success flag to fill in into the assignment table for all simulations", True)
definition.add_flag("show", "show the assignment", False)

# Create configuration
config = parse_arguments("create_assignment", definition, "Create simulation assignment file from simulation names")

# -----------------------------------------------------------------

# Create simulation manager
manager = SimulationManager()

# Set options
manager.config.simulation_names = config.names
manager.config.remotes = config.remotes
manager.config.from_directories = config.from_directories
manager.config.local = config.local
manager.config.success = config.success
manager.config.show_assignment = config.show
manager.config.write_assignment = True

# Run the manager
manager.run()

# -----------------------------------------------------------------
