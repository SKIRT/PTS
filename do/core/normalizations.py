#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.normalizations Show the normalizations of stellar components in a ski file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.skifile import SkiFile
from pts.core.simulation.skifile import show_normalizations

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Ski path
definition.add_required("ski", "file_path", "path of the ski file")

# Optional
definition.add_optional("flux", "photometric_unit", "also show flux in a particular unit")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("normalizations", definition, "Show the normalizations of stellar components in a ski file")

# -----------------------------------------------------------------

# Load the ski file
ski = SkiFile(config.ski)

# -----------------------------------------------------------------

show_normalizations(ski, flux_unit=config.flux)

# -----------------------------------------------------------------
