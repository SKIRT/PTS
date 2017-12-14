#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import parsing
from pts.core.basics.emissionlines import get_identifiers, important_lines, strong_lines
from pts.magic.tools import wavelengths

# -----------------------------------------------------------------

# Npoints
default_nwavelengths_range = "50>400"

# Wavelength range
default_min_wavelength = "0.019 micron"
default_max_wavelength = "2050. micron"

# Set default wavelength range
default_wavelength_range = "0.019 micron > 2050. micron"

# -----------------------------------------------------------------

# Filters
default_filter_names = "GALEX,SDSS,2MASS,IRAC,MIPS 24mu,Herschel"
default_filters = parsing.lazy_filter_list(default_filter_names)

# -----------------------------------------------------------------

filter_wavelengths = [fltr.wavelength for fltr in default_filters]

# Emission lines
all_lines = get_identifiers()
#default_lines = important_lines
default_lines = strong_lines

# -----------------------------------------------------------------

# FOR PLOTTING:

# SEDs
seds = ["mappings", "bruzual_charlot"]
default_ages = "8.0 Gyr,0.1 Gyr"

# Wavelength regimes
regimes = wavelengths.all_regimes()
default_regimes = ["euv", "fuv", "muv", "nuv", "optical", "nir", "mir", "fir", "submm", "microwave"]

# -----------------------------------------------------------------

default_min_points_per_filter = 8
default_min_points_per_fwhm = 5

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("ngrids", "positive_integer", "number of wavelength grids to generate")
definition.add_positional_optional("npoints_range", "integer_range", "range of the number of wavelength points for the grids", default=default_nwavelengths_range)
definition.add_positional_optional("range", "length_quantity_range", "wavelength range", default_wavelength_range, convert_default=True)

# Emission lines
definition.add_flag("add_emission_lines", "add emission lines", True)
definition.add_optional("emission_lines", "string_list", "emission lines to use for building the grid", default=default_lines, choices=all_lines)

# Filters
definition.add_optional("check_filters", "filter_list", "check coverage for these filters")
definition.add_flag("adjust_minmax", "adjust min or max based on the filters to be checked", False)
definition.add_optional("filters", "filter_list", "resample for these filters", default_filters)
definition.add_optional("adjust_to", "length_quantity_list", "adjust wavelength points to these exact wavelengths", filter_wavelengths)
definition.add_optional("fixed", "length_quantity_list", "fixed wavelengths to add to the grid")

# Add optional
definition.add_optional("min_wavelengths_in_filter", "positive_integer", "minimum number of wavelength points to sample a filter", default_min_points_per_filter)
definition.add_optional("min_wavelengths_in_fwhm", "positive_integer", "minimum number of wavelength points to sample within inner range of filter", default_min_points_per_fwhm)

# Add flags
definition.add_flag("show", "show", False)
definition.add_flag("plot", "plot", False)
definition.add_flag("write", "write", True)
definition.add_flag("write_grids", "write grids", True)
definition.add_flag("write_table", "write table", True)
definition.add_flag("write_elements", "write elements", True)

# Table
definition.add_flag("table", "create a table if not passed as input", True)

# Output and plot directories
definition.add_optional("output", "directory_path", "output directory")
definition.add_optional("plot_path", "directory_path", "plotting directory")

# -----------------------------------------------------------------

# Plotting

definition.add_flag("plot_reference", "plot reference grid and differences", False)

# Regimes?
definition.add_flag("plot_regimes", "show wavelength regimes", True)
definition.add_optional("regimes", "string_list", "regimes to show", default=default_regimes, choices=regimes)
definition.add_flag("only_subregimes", "only show subregimes", True)

# Filters
definition.add_flag("plot_filters", "add filters", True)
definition.add_flag("categorize_filters", "categorize filters per instrument", True)
definition.add_optional("plotting_filters", "broad_band_filter_list", "filters to show on the plot (default is filters used for resampling the grid)")

# Emission lines
definition.add_flag("plot_lines", "add emission lines", True)

# Add SEDs
definition.add_flag("plot_seds", "add template SEDs", False)
definition.add_optional("seds", "string_list", "template SEDs to add", default=seds, choices=seds)
definition.add_optional("metallicity", "positive_real", "metallicity for templates", 0.02)

# Mappings template
definition.add_optional("compactness", "positive_real", "compactness for MAPPINGS template", 6.)
definition.add_optional("pressure", "quantity", "pressure for MAPPINGS template", "1e12 K/m3", convert_default=True)
definition.add_optional("covering_factor", "positive_real", "covering factor for MAPPINGS template", 0.2)

# Bruzual-Charlot
definition.add_optional("ages", "time_quantity_list", "ages for Bruzual-Charlot templates", default_ages, convert_default=True)

# Plot residuals and resampled SEDs
definition.add_flag("plot_resampled", "plot the SEDs resampled on the complete wavelength grid", False)
definition.add_flag("plot_interpolated", "plot the SEDs interpolated", False)
definition.add_flag("plot_residuals", "plot the residuals between the original SEDs and the resampled SEDs", False)

# -----------------------------------------------------------------
