#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.config.plot import definition as plot_definition
from pts.core.basics.plot import plotting_libraries

# -----------------------------------------------------------------

default_library = "matplotlib"

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_positional_optional("seds", "filepath_list", "SED files to be plotted")

# Add plotting options
definition.import_section("plot", "plotting options", plot_definition)

# The unit in which to plot
definition.add_optional("wavelength_unit", "length_unit", "unit of wavelength", "micron", convert_default=True)
definition.add_optional("unit", "photometric_unit", "photometric unit", "Jy", convert_default=True)

# The plotting library to use
definition.add_optional("library", "string", "plotting library", default_library, plotting_libraries)

# -----------------------------------------------------------------
