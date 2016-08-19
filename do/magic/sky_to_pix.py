#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.sky_to_pix Convert a sky coordinate to a pixel coordinate for a specific WCS.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

definition = ConfigurationDefinition()

definition.add_required("coordinate", "skycoordinate", "the sky coordinate")
definition.add_required("wcs_path", "file_path", "the path to the file holding the WCS info")

# Get the configuration
setter = ArgumentConfigurationSetter("pix_to_sky")
config = setter.run(definition)

# -----------------------------------------------------------------

# Print the pixel coordinate
print(config.coordinate.to_pixel(CoordinateSystem.from_file(config.wcs_path)))

# -----------------------------------------------------------------
