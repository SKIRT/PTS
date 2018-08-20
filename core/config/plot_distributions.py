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
normalizations = ["max", "sum"]

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Distributions
definition.add_positional_optional("distributions", "filepath_list", "distribution files to be plotted")
definition.add_flag("panels", "plot the distributions in separate panels")
definition.add_flag("recursive", "search distribution files recursively in the working directory")
definition.add_optional("normalize", "string", "normalize all distributions by a certain method", choices=normalizations)
definition.add_optional("normalization_value", "real", "value for normalization", 1.)

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)
definition.sections["plot"].optional["xsize"].default = 8
definition.sections["plot"].optional["ysize"].default = 4
definition.add_flag("logscale", "use value log scale")
definition.add_flag("logfrequency", "use log scale for frequency")
definition.add_optional("bar_width", "positive_real", "relative width of the bars (1 means edges touch)", 1.)
definition.add_flag("use_name_xlabel", "use the distribution name(s) for the x labels of the panels")
definition.add_flag("colours_per_panel", "reuse the same colours for each panel")

# Add features
definition.add_flag("smooth", "add smooth curves to plot")
definition.add_flag("statistics", "add statistics to plot")
definition.add_flag("extrema", "add extrema to plot")
definition.add_flag("minima", "add minima", None)
definition.add_flag("maxima", "add maxima", None)
definition.add_flag("hatches", "add hatches", False)
definition.add_optional("y_label", "string", "label for vertical axes")
definition.add_flag("distribution_ticks", "use the distribution values as the horizontal axes ticks")
definition.add_flag("y_ticks", "show the y ticks", True)
definition.add_optional("lower_than", "positive_real", "mark frequencies/counts lower than this amount")
definition.add_optional("higher_than", "positive_real", "mark frequencies/counts higher than this amount")
definition.add_flag("add_lower_frequencies", "add frequencies for marked lower")
definition.add_flag("add_higher_frequencies", "add frequencies for marked higher")
definition.add_flag("lower_frequency_percentages", "display frequencies as percentages")
definition.add_flag("higher_frequency_percentages", "display frequencies as percentages")

# Show
definition.add_flag("edges", "show edges")
definition.add_flag("frequencies", "show frequencies")

# Add optional
definition.add_optional("output", "string", "output directory")
definition.add_optional("format", "string", "plotting format", default=default_format, choices=formats)

# Alpha
definition.add_optional("alpha", "real", "alpha")

# -----------------------------------------------------------------

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------
