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

default_residual_reference = "models"
residual_references = ["models", "observations"]

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# SEDs from file
definition.add_positional_optional("seds", "filepath_list", "SED files to be plotted")
definition.add_flag("multi", "look for multiple SEDs per file")
definition.add_optional("wavelength_unit_file", "length_unit", "wavelength unit in SED file")
definition.add_optional("unit_file", "photometric_unit", "photometric unit in SED file")
definition.add_optional("distance", "length_quantity", "object distance (for flux <> luminosity unit conversion")

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)

# Add optional
definition.add_optional("output", "string", "output directory")
definition.add_optional("format", "string", "plotting format", default=default_format, choices=formats)

# -----------------------------------------------------------------

# The unit in which to plot
definition.add_optional("wavelength_unit", "length_unit", "unit of wavelength", "micron", convert_default=True)
definition.add_optional("unit", "photometric_unit", "photometric unit", "Jy", convert_default=True)

# Residual reference
definition.add_optional("residual_reference", "string", "reference for the residuals", default_residual_reference, choices=residual_references)

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

# Ignore these filters
definition.add_optional("ignore_filters", "lazy_filter_list", "ignore these filters from the observed SEDs")

# -----------------------------------------------------------------

definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------

definition.add_flag("ignore_upper", "ignore upper limit points", False)

# -----------------------------------------------------------------
