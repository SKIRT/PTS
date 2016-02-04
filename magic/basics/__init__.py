#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.magic.basics TO DO
#
# This package ...
#

# -----------------------------------------------------------------

# Import classes to make them available at the level of this subpackage
from .vector import Position, Extent
from .layers import Layers
from .mask import Mask
from .region import Region
from .trackrecord import TrackRecord
from .geometry import Ellipse, Rectangle
from .catalogcoverage import CatalogCoverage
from .coordinatesystem import CoordinateSystem
