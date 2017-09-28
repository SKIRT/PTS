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

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add positional optional
definition.add_positional_optional("host_ids", "string_list", "name of the remote host for which to show/retrieve simulations and tasks", choices=find_host_ids())

# Add flag
definition.add_flag("show_progress", "show the progress of the simulation that is still running (only if there is just one)", False)
definition.add_flag("debug_output", "show all simulation output in debug mode")

# -----------------------------------------------------------------
