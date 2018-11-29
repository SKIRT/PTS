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
from pts.core.tools import sequences

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

# Legends
default_instruments_loc = "lower center"
default_observations_loc = "upper left"
default_models_loc = "upper right"

# Locations
locations = sequences.get_lists_combinations(["upper", "center", "lower"], ["left", "center", "right"], combine=lambda x: " ".join(x) if x[0] != x[1] else x[0])

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

#definition.add_flag("old", "old plotting methods", False)

# SEDs from file
definition.add_positional_optional("seds", "filepath_list_or_string_filepath_dictionary", "SED files to be plotted")
definition.add_flag("multi", "look for multiple SEDs per file", False)
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
definition.add_flag("observations_residuals", "plot residuals between observed SEDs", True)
definition.add_flag("models_residuals", "plot residuals between model SEDs", False)

# Minimum and maximum reference
definition.add_optional("minmax_wavelength_reference", "string", "reference for determining the minimum and maximum of the wavelength axis (default is both models and observations)", choices=residual_references)
definition.add_optional("minmax_photometry_reference", "string", "reference for determining the minimum and maximum of the photometry axis (default is both models and observations)", choices=residual_references)

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

# Flags for residuals
definition.add_flag("interpolate_models_for_residuals", "interpolate models to get the photometry at specific wavelengths", True)
definition.add_flag("smooth_residuals", "splot a spline interpolation of the residual of the models when they are plotted against observations as reference")
definition.add_flag("show_smooth", "plot the spline interpolation of the observation")
definition.add_optional("smooth_residuals_closeness_limit", "real", "closeness of observed bands", default=0.01)
definition.add_flag("smooth_residuals_ignore_close", "ignore close observed bands for smooth interpolation")

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

# Min and max
definition.add_optional("min_wavelength", "length_quantity", "minimum wavelength")
definition.add_optional("max_wavelength", "length_quantity", "maximum wavelength")
definition.add_optional("min_flux", "photometric_quantity", "minimum flux")
definition.add_optional("max_flux", "photometric_quantity", "maximum flux")

# -----------------------------------------------------------------

# TeX rendering of text
definition.add_flag("tex", "enable TeX rendering", True)

# -----------------------------------------------------------------

# Axis positions
definition.add_optional("xaxis_position", "string", "position of x axis ticks and label", "bottom", choices=["bottom", "top"])
definition.add_optional("yaxis_position", "string", "position of y axis ticks and label", "left", choices=["left", "right"])

# -----------------------------------------------------------------

# Legends positions
definition.add_section("legends", "legend options")

# Add?
definition.sections["legends"].add_flag("instruments", "add instruments legend", True)
definition.sections["legends"].add_flag("observations", "add observations legend", True)
definition.sections["legends"].add_flag("models", "add models legend", True)
definition.sections["legends"].add_flag("residuals", "add residuals legends", True)

# LOC
definition.sections["legends"].add_optional("instruments_location", "string", "location of instruments legend", default_instruments_loc, choices=locations)
definition.sections["legends"].add_optional("observations_location", "string", "location of observations legend", default_observations_loc, choices=locations)
definition.sections["legends"].add_optional("models_location", "string", "location of models legend", default_models_loc, choices=locations)

# Residuals
definition.sections["legends"].add_optional("observations_residuals_location", "string", "location of observations legend on residuals panels", default_observations_loc, choices=locations)
definition.sections["legends"].add_optional("models_residuals_location", "string", "location of models legend on residuals panels", default_models_loc, choices=locations)

# NCOL
definition.sections["legends"].add_optional("instruments_ncols", "positive_integer", "number of columns for instruments legend", 2)
definition.sections["legends"].add_optional("observations_ncols", "positive_integer", "number of columns for observations legend", 2)
definition.sections["legends"].add_optional("models_ncols", "positive_integer", "number of columns for models legend", 2)

# Residuals
definition.sections["legends"].add_optional("observations_residuals_ncols", "positive_integer", "number of columns for observations legend on residuals panels", 2)
definition.sections["legends"].add_optional("models_residuals_ncols", "positive_integer", "number of columns for models legend on residuals panels", 2)

# -----------------------------------------------------------------

# For
definition.add_flag("only_residuals_legend", "for SEDs for which only_residuals is enabled, show only in a seperate legend on the appropriate residuals panel, instead of in the list on the main panel")

# -----------------------------------------------------------------

definition.add_optional("single_observation_cmap", "string", "colormap for if only one observation", default="rainbow")

# -----------------------------------------------------------------
