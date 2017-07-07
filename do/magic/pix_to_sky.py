#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.pix_to_sky Convert a pixel coordinate to a sky coordinate for a specific WCS.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

definition = ConfigurationDefinition()

definition.add_required("coordinate", "pixelcoordinate", "the pixel coordinate")
definition.add_required("wcs_path", "file_path", "the path to the file holding the WCS info")

# Get the configuration
config = parse_arguments("pix_to_sky", definition)

# -----------------------------------------------------------------

# Print the sky coordinate
print(config.coordinate.to_sky(CoordinateSystem.from_file(config.wcs_path)))

# -----------------------------------------------------------------
