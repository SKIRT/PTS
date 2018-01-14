#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.config.plot import definition as plot_definition
from pts.core.filter.broad import identifiers, BroadBandFilter

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add positional optional
definition.add_positional_optional("filters", "lazy_filter_list", "broad band filters for which to plot the transmission and narrow band filters to plot as lines", default=map(BroadBandFilter, identifiers.keys()))
definition.add_positional_optional("wavelengths", "quantity_list", "wavelengths to plot as lines on the transmission curves")

# Add flags
definition.add_flag("emission", "add emission lines")

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)

# Change defaults
definition.sections["plot"].optional["xsize"].default = 17
definition.sections["plot"].optional["ysize"].default = 4

# Add optional
definition.add_optional("title", "string", "plot title")
definition.add_optional("output", "string", "output plot file path")

# Add flags
definition.add_flag("write", "write", False)

# -----------------------------------------------------------------
