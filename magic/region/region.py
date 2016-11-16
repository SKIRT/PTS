#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.region Contains the (abstract) Region class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class Region(object):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set the 'exclude' flag
        self.include = kwargs.pop("include", True)

        # Set the appearance info
        self.appearance = kwargs.pop("appearance", dict())

        # Set the meta information
        self.meta = kwargs.pop("meta", dict())

# -----------------------------------------------------------------

class PixelRegion(Region):

    """
    This class ...
    """

    pass

# -----------------------------------------------------------------

class SkyRegion(Region):

    """
    This class ...
    """

    pass

# -----------------------------------------------------------------

class PhysicalRegion(Region):

    """
    This class ...
    """

    pass

# -----------------------------------------------------------------
