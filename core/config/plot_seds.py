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
from pts.core.plot.sed import residual_references

# -----------------------------------------------------------------

formats = ["pdf", "png"]
default_format = "pdf"

# -----------------------------------------------------------------

default_residual_reference = "models"

# -----------------------------------------------------------------

# SEDs
seds = ["mappings", "bruzual_charlot"]
default_ages = "8.0 Gyr,0.1 Gyr"

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
definition.add_flag("models_residuals", "plot residuals for only model SEDs", False)

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

# Add SEDs
definition.add_flag("add_templates", "add template SEDs")
definition.add_optional("templates", "string_list", "template SEDs to add", default=seds, choices=seds)
definition.add_optional("metallicity", "positive_real", "metallicity for templates", 0.02)

# Mappings template
definition.add_optional("compactness", "positive_real", "compactness for MAPPINGS template", 6.)
definition.add_optional("pressure", "quantity", "pressure for MAPPINGS template", "1e12 K/m3", convert_default=True)
definition.add_optional("covering_factor", "positive_real", "covering factor for MAPPINGS template", 0.2)

# Bruzual-Charlot
definition.add_optional("ages", "time_quantity_list", "ages for Bruzual-Charlot templates", default_ages, convert_default=True)

# -----------------------------------------------------------------

definition.add_flag("interpolate_models_for_residuals", "interpolate models to get the photometry at specific wavelengths", True)

# -----------------------------------------------------------------

# Load only certain files
definition.add_optional("contains", "string", "only load SED files containing this string in their name")
definition.add_optional("not_contains", "string", "don't load SED files containing this string in their name")
definition.add_optional("exact_name", "string", "only load SED files with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't load SED files with this exact string as their name")
definition.add_optional("startswith", "string", "only load SED files whose name starts with this string")
definition.add_optional("endswith", "string", "only load SED files whose name starts with this string")

# -----------------------------------------------------------------

# Additional relative error
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")
definition.add_optional("additional_for", "string_list", "observed SED labels on which to apply the additional error")
definition.add_optional("additional_not_for", "string_list", "observed SED labels on which to not apply the additional error")

# -----------------------------------------------------------------
