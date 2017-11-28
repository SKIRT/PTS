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
from pts.core.prep.wavelengthgrids import subgrids
from pts.magic.tools import wavelengths

# -----------------------------------------------------------------

#seds = ["mappings", "bruzual_charlot"]
seds = ["mappings"]
default_ages = "8.0 Gyr,0.1 Gyr"

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Grids
definition.add_positional_optional("grids", "filepath_list", "wavelength grids to be plotted")
definition.add_flag("load_subgrids", "plot from the output of subgrid wavelength grid generation in the working directory")
definition.add_flag("create_subgrids", "generate subgrid wavelength grid")
definition.add_optional("subgrids", "string_list", "subgrids to load", default=subgrids, choices=subgrids)

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints", "positive_integer", "range of the wavelength grid size", 100)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", True)
definition.sections["wg"].add_optional("range", "quantity_range", "range of wavelengths", "0.02 micron > 2000 micron", convert_default=True)
definition.sections["wg"].add_optional("fixed", "length_quantity_list", "fixed wavelengths")
definition.sections["wg"].add_optional("filters", "filter_list", "resample for these filters")
definition.sections["wg"].add_optional("adjust_to", "length_quantity_list", "adjust wavelength points to these exact wavelengths")
definition.sections["wg"].add_optional("check_filters", "filter_list", "check coverage for these filters")

# -----------------------------------------------------------------

default_min_wavelength = "0.019 micron"
default_max_wavelength = "2050. micron"
definition.add_optional("min_wavelength", "length_quantity", "minimum plot wavelength", default_min_wavelength, convert_default=True)
definition.add_optional("max_wavelength", "length_quantity", "maximum plot wavelength", default_max_wavelength, convert_default=True)

# -----------------------------------------------------------------

# Regimes?
regimes = wavelengths.all_regimes()
default_regimes = ["euv", "fuv", "muv", "nuv", "optical", "nir", "mir", "fir", "submm", "microwave"]
definition.add_flag("add_regimes", "show wavelength regimes")
definition.add_optional("regimes", "string_list", "regimes to show", default=default_regimes, choices=regimes)

# Filters
default_filter_names = "GALEX,SDSS,2MASS,IRAC,Herschel"
definition.add_flag("add_filters", "add filters")
definition.add_optional("filters", "lazy_filter_list", "filters to add", default=default_filter_names, convert_default=True)
definition.add_flag("categorize_filters", "categorize filters per instrument", False)

# Emission lines
definition.add_flag("add_lines", "add emission lines")

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

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)
definition.sections["plot"].optional["figsize"].default = (11,5)

# The unit in which to plot
definition.add_optional("wavelength_unit", "length_unit", "unit of wavelength", "micron", convert_default=True)

# Residual reference
#definition.add_optional("residual_reference", "string", "reference for the residuals", default_residual_reference, choices=residual_references)

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

# Ignore these filters
#definition.add_optional("ignore_filters", "filter_list", "ignore these filters from the observed SEDs")

# -----------------------------------------------------------------

definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------
