#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.relaunch_simulations_screen Relaunch simulations in a screen from a local script file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.screen import ScreenScript

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "script file path")

# Optional arguments
definition.add_optional("parallelization", "parallelization", "new parallelization scheme")

# Read the command line arguments
config = parse_arguments("relaunch_simulations_screen", definition, description="Relaunch simulations in a screen from a local script file")

# -----------------------------------------------------------------

# Load screen script
screen = ScreenScript.from_file(config.filename)

# -----------------------------------------------------------------
