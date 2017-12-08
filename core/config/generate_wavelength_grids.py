#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import parsing
from pts.core.basics.emissionlines import get_id_strings, important_lines
from pts.magic.tools import wavelengths

# -----------------------------------------------------------------

# Npoints
default_nwavelengths_range = "100>500"

# Wavelength range
default_min_wavelength = "0.019 micron"
default_max_wavelength = "2050. micron"

# Set default wavelength range
default_wavelength_range = "0.019 micron > 2050. micron"

# -----------------------------------------------------------------

# Filters
default_filter_names = "GALEX,SDSS,2MASS,Spitzer,Herschel"
default_filters = parsing.lazy_filter_list(default_filter_names)

# Emission lines
all_lines = get_id_strings()
default_lines = important_lines

# -----------------------------------------------------------------

# FOR PLOTTING:

# SEDs
seds = ["mappings", "bruzual_charlot"]
default_ages = "8.0 Gyr,0.1 Gyr"

# Wavelength regimes
regimes = wavelengths.all_regimes()
default_regimes = ["euv", "fuv", "muv", "nuv", "optical", "nir", "mir", "fir", "submm", "microwave"]

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("ngrids", "positive_integer", "number of wavelength grids to generate")
definition.add_positional_optional("npoints_range", "integer_range", "range of the number of wavelength points for the grids", default=default_nwavelengths_range)
definition.add_positional_optional("range", "length_quantity_range", "wavelength range", default_wavelength_range, convert_default=True)

# Add optional
definition.add_optional("min_wavelengths_in_filter", "positive_integer", "minimum number of wavelength points to sample a filter", 5)
definition.add_optional("min_wavelengths_in_fwhm", "positive_integer", "minimum number of wavelength points to sample within inner range of filter", 3)

# Add flags
definition.add_flag("show", "show", False)
definition.add_flag("plot", "plot", False)
definition.add_flag("write", "write", False)

# -----------------------------------------------------------------

# Plotting

# Regimes?
definition.add_flag("add_regimes", "show wavelength regimes")

definition.add_optional("regimes", "string_list", "regimes to show", default=default_regimes, choices=regimes)
definition.add_flag("only_subregimes", "only show subregimes", True)

# Filters
definition.add_flag("add_filters", "add filters")
definition.add_optional("filters", "lazy_filter_list", "filters to add", default=default_filters)
definition.add_flag("categorize_filters", "categorize filters per instrument", True)

# Emission lines
definition.add_flag("add_lines", "add emission lines")
definition.add_optional("lines", "string_list", "emission line IDs to plot", default=default_lines, choices=all_lines)

# Add SEDs
definition.add_flag("add_seds", "add template SEDs")
definition.add_optional("seds", "string_list", "template SEDs to add", default=seds, choices=seds)
definition.add_optional("metallicity", "positive_real", "metallicity for templates", 0.02)

# Mappings template
definition.add_optional("compactness", "positive_real", "compactness for MAPPINGS template", 6.)
definition.add_optional("pressure", "quantity", "pressure for MAPPINGS template", "1e12 K/m3", convert_default=True)
definition.add_optional("covering_factor", "positive_real", "covering factor for MAPPINGS template", 0.2)

# Bruzual-Charlot
definition.add_optional("ages", "time_quantity_list", "ages for Bruzual-Charlot templates", default_ages, convert_default=True)

# -----------------------------------------------------------------
