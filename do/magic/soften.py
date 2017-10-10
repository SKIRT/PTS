#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.soften Soften edges of an image.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.rgba import RGBAImage
from pts.magic.region.ellipse import PixelEllipseRegion
from pts.magic.basics.stretch import PixelStretch
from pts.magic.basics.coordinate import PixelCoordinate
from pts.core.basics.range import RealRange

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_required("file_path", "file_path", "name of the input image file")

# Parse the command line arguments
config = parse_arguments("soften", definition)

# -----------------------------------------------------------------

image = RGBAImage.from_file(config.filepath)

# -----------------------------------------------------------------

radius = PixelStretch(500.,300.)
center = PixelCoordinate(750., 1200.)
ellipse = PixelEllipseRegion(center, radius)

# -----------------------------------------------------------------

factor_range = RealRange(0.4, 1.2)
image.soften_edges(ellipse, factor_range)

# -----------------------------------------------------------------

image.show()

# -----------------------------------------------------------------
