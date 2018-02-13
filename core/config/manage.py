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

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Assignment file
definition.add_positional_optional("assignment", "file_path", "path of assignment file")

# Remotes for which to find corresponding simulations
definition.add_optional("remotes", "string_list", "remote hosts for which to look for matching simulations (not necessary when assignment is specified)", default=host_ids, choices=host_ids)

# Timing and memory table
definition.add_optional("timing", "file_path", "timing table path")
definition.add_optional("memory", "file_path", "memory table path")

# -----------------------------------------------------------------

# Flags
definition.add_flag("local", "treat simulations without a match as local simulations", False)
definition.add_flag("warn_local", "give a warning for each simulation that didn't have a match for any host ID")
definition.add_flag("success", "success flag to fill in into the assignment table for all simulations", True)

# -----------------------------------------------------------------

# Showing
definition.add_flag("show", "showing", True)
definition.add_flag("show_assignment", "show the assignment scheme")
definition.add_flag("show_status", "show the simulation status")
definition.add_flag("show_runtimes", "show runtimes")
definition.add_flag("show_memory", "show memory")

# Plotting
definition.add_flag("plot", "plotting", False)
definition.add_flag("plot_runtimes", "plot runtimes")
definition.add_flag("plot_memory", "plot memory usage")

# Writing
definition.add_flag("write", "writing", False)

# Analysis
definition.add_flag("analyse", "analysis", False)

# -----------------------------------------------------------------
