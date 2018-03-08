#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Flags
definition.add_flag("recursive", "look for simulation in directories recursively", True)

# Output
definition.add_optional("output", "directory_path", "output directory", letter="o")

# Flags
definition.add_flag("other", "also plot the 'other' phases")
definition.add_flag("group", "group timelines for the same number of processes in the same plot")

# Plot
definition.add_optional("label_fontsize", "positive_integer", "fontsize for the axes labels", default=18)
definition.add_optional("figsize", "integer_pair", "size of the figure", default=(12,8))
definition.add_flag("percentages", "show percentages")
definition.add_flag("totals", "show totals")
definition.add_optional("title", "string", "plot title")
definition.add_flag("add_border", "add plot border")
definition.add_flag("show_ranks", "show process ranks")
definition.add_optional("ticks_fontsize", "positive_integer", "fontsize of the axes ticks", default=12)

# -----------------------------------------------------------------

definition.add_flag("write_timelines", "write the timelines", True)
definition.add_flag("write_data", "write the plotting data", True)

# -----------------------------------------------------------------
