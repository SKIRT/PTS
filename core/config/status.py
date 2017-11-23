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
definition.add_positional_optional("host_ids", "string_list", "name of the remote host(s) for which to show/retrieve simulations and tasks", choices=find_host_ids())
definition.add_positional_optional("ids", "string_integer_list_dictionary", "simulation IDs for the different remote hosts")

# Add flag
definition.add_flag("show_progress", "show the progress of the simulation that is still running (only if there is just one)", False)
definition.add_flag("debug_output", "show all simulation output in debug mode")

# -----------------------------------------------------------------

definition.add_flag("retrieve", "retrieve finished simulations", True)
definition.add_flag("analyse", "analyse retrieved simulations", True)

# -----------------------------------------------------------------

# Crashed
definition.add_flag("check_crashed", "check whether crashed simulations have the necessary output, if so, retrieve them")
definition.add_optional("retrieve_crashed", "string_integer_list_dictionary", "retrieve crashed simulations for these hosts and simulation IDs")
definition.add_flag("check_data", "for crashed simulations, check whether the simulation data is valid")

# -----------------------------------------------------------------

definition.add_flag("offline", "run in offline mode: only analyse already retrieved simulations and tasks, don't try to connect to remotes")

# -----------------------------------------------------------------
