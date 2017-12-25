#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.emissionlines import get_identifiers
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Emission lines
all_lines = get_identifiers()
default_lines = None # no default: means all

# -----------------------------------------------------------------

default_min_points_per_filter = 8
default_min_points_per_fwhm = 5

# -----------------------------------------------------------------

definition = definition.copy()

# Output directory
definition.add_optional("output", "directory_path", "output directory", letter="o")

# Flags
definition.add_flag("write", "do writing", True)
definition.add_flag("plot", "do plotting", True)
definition.add_flag("write_grid", True)

# -----------------------------------------------------------------

definition.add_required("npoints", "positive_integer", "number of aimed points in the grid")

# -----------------------------------------------------------------

definition.add_flag("add_emission_lines", "add emission lines", False)
definition.add_optional("emission_lines", "string_list", "emission lines to use for building the grid", default=default_lines, choices=all_lines)
definition.add_optional("min_wavelength", "length_quantity", "minimum wavelength", "0.05 micron", convert_default=True)
definition.add_optional("max_wavelength", "length_quantity", "maximum wavelength", "2000 micron", convert_default=True)
definition.add_optional("check_filters", "filter_list", "check coverage for these filters")
definition.add_flag("adjust_minmax", "adjust min or max based on the filters to be checked", False)
definition.add_optional("filters", "filter_list", "resample for these filters")
definition.add_optional("adjust_to", "length_quantity_list", "adjust wavelength points to these exact wavelengths")
definition.add_optional("fixed", "length_quantity_list", "fixed wavelengths to add to the grid")

# -----------------------------------------------------------------

# Add optional
definition.add_optional("min_wavelengths_in_filter", "positive_integer", "minimum number of wavelength points to sample a filter", default_min_points_per_filter)
definition.add_optional("min_wavelengths_in_fwhm", "positive_integer", "minimum number of wavelength points to sample within inner range of filter", default_min_points_per_fwhm)

# -----------------------------------------------------------------
