#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.core.basics.plot import mpl, plotting_libraries, pdf, plotting_formats

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Different plotting features
definition.add_flag("instruments", "plot a comparison between the SEDs of the different instruments", True)
definition.add_flag("contributions", "plot the various contributions to the SEDS", True)

# The output directory
definition.add_optional("output", "directory_path", "output directory", fs.cwd())

# The unit in which to plot
definition.add_optional("wavelength_unit", "length_unit", "unit of wavelength", "micron", convert_default=True)
definition.add_optional("unit", "photometric_unit", "photometric unit", "Jy", convert_default=True)

# -----------------------------------------------------------------

# The plotting format
definition.add_optional("format", "string", "plotting format", pdf, plotting_formats)

# The plotting library to use
definition.add_optional("library", "string", "plotting library", mpl, plotting_libraries)

# -----------------------------------------------------------------

# Reference SEDS
definition.add_optional("reference_seds", "filepath_list", "paths of reference SEDs")

# Ignore these filters
definition.add_optional("ignore_filters", "filter_list", "ignore these filters from the observed SEDs")

# -----------------------------------------------------------------
