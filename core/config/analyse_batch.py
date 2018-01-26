#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_positional_optional("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "integer", "ID of the simulation")

# Add flags
definition.add_flag("ignore_missing_data", "ignore missing data when analysing the simulations", False)

# Replace already existing simulations
definition.add_flag("replace", "replace the timing and memory information when a simulation is already in the table (otherwise, an error is raised)", False)

# -----------------------------------------------------------------
