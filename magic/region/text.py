#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.text Contains the TextRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion
from ..basics.mask import Mask

# -----------------------------------------------------------------

class TextRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set the attributes
        self.center = center
        self.text = text

        # Call the constructor of the base class
        super(TextRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelTextRegion(TextRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be pixel coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

# -----------------------------------------------------------------

class SkyTextRegion(TextRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param ra:
        :param dec:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be sky coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

# -----------------------------------------------------------------

class PhysicalTextRegion(TextRegion, PhysicalCoordinate):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param center:
        :param text:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be physical coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

# -----------------------------------------------------------------
