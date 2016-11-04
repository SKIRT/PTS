#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("ski", "file_path", "path to the ski file")

# Add optional
definition.add_optional("input", "directory_path", "path to the input directory")

# Add optional
definition.add_optional("nwavelengths", "integer", "the number of wavelengths (useful for when a file wavelength grid is used and the file is not present)")
definition.add_optional("ncells", "integer", "number of dust cells (useful only when the ski file includes a tree dust grid)")

# Add flag
definition.add_flag("probe", "probe the number of cells (or memory usage) by launching a dry-run of the simulation, if necessary (tree dust grid and ncells not specified)", True)

# Add optional
definition.add_optional("nprocesses", "integer_list", "number of processes", [1, 2, 4, 8, 12, 20])

# Add flag
#definition.add_flag("data_parallel", "with data parallelization mode")

# Flags
definition.add_flag("show", "show the results", True)
definition.add_flag("plot", "plot the expected memory scaling", False)

# -----------------------------------------------------------------
