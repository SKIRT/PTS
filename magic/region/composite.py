#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.composite Contains the CompositeRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .region import Region

# -----------------------------------------------------------------

class CompositeRegion(Region):

    """
    This class ...
    """

# -----------------------------------------------------------------

class PixelCompositeRegion(CompositeRegion):

# -----------------------------------------------------------------

class SkyCompositeRegion(CompositeRegion):

# -----------------------------------------------------------------

class PhysicalCompositeRegion(CompositeRegion):

# -----------------------------------------------------------------
