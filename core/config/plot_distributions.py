#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.config.plot import definition as plot_definition
from pts.core.basics.plot import plotting_libraries, mpl

# -----------------------------------------------------------------

formats = ["pdf", "png"]
default_format = "pdf"

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Distributions from file
definition.add_positional_optional("distributions", "filepath_list", "distribution files to be plotted")
definition.add_flag("panels", "plot the distributions in separate panels")

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)
definition.sections["plot"].optional["figsize"].default = (8,4)
definition.add_flag("logscale", "use value log scale")
definition.add_flag("logfrequency", "use log scale for frequency")

# Add features
definition.add_flag("smooth", "add smooth curves to plot")
definition.add_flag("statistics", "add statistics to plot")
definition.add_flag("extrema", "add extrema to plot")

# Show
definition.add_flag("edges", "show edges")
definition.add_flag("frequencies", "show frequencies")

# Add optional
definition.add_optional("output", "string", "output directory")
definition.add_optional("format", "string", "plotting format", default=default_format, choices=formats)

# -----------------------------------------------------------------

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------
