#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.tools.wavelengths import all_regimes

# -----------------------------------------------------------------

filter_types = ["narrow", "broad"]

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_positional_optional("types", "string_list", "show only narrow or broad band filters", choices=filter_types, default=filter_types)
definition.add_positional_optional("matching", "string", "only show filters matching with this string")
definition.add_optional("regime", "string", "wavelength regime for which to list filters", choices=all_regimes())
definition.add_flag("show", "show", True)
definition.add_flag("short", "short: only show the filter names", letter="s")
definition.add_flag("aliases", "show aliases", letter="a")
definition.add_flag("categorize", "categorize per instrument/observatory/system", True, letter="c")
definition.add_flag("plot", "plot the filter transmission curves", letter="p")
definition.add_flag("plot_fwhm", "plot the FWHM ranges of the filters")

# -----------------------------------------------------------------
