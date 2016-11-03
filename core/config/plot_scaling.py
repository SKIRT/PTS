#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.plot.scaling import scaling_properties, simulation_phases

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Flags
definition.add_flag("recursive", "look for simulation in directories recursively", True)

# Add optional
definition.add_positional_optional("properties", "string_list", "plot the scaling of these properties", choices=scaling_properties, default=scaling_properties)
definition.add_positional_optional("phases", "string_list", "the simulation phases for which to do the plotting", choices=simulation_phases, default=["total"])

definition.add_flag("hybridisation", "plot as a function of number of processes for constant number of cores")

definition.add_optional("output", "directory_path", "output directory", letter="o")
definition.add_optional("figsize", "integer_tuple", "figure size", default=(12,8))

definition.add_optional("sigma_level", "real", "sigma level for plotting error bars", 1.0)

definition.add_flag("fit", "fit theoretical curves to timing and memory data", True)
definition.add_flag("plot_fit", "plot the fitted relations", True)

# -----------------------------------------------------------------
