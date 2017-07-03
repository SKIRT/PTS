#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.dimensions Contains the Dimensions class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class Dimensions(object):

    """
    This class ...
    """

    def __init__(self, nx, ny):

        """
        The constructor ...
        """

        # Check if integer
        if not isinstance(nx, int): raise ValueError("nx must be integer")
        if not isinstance(ny, int): raise ValueError("ny must be integer")

        # Set the attributes
        self.nx = nx
        self.ny = ny

# -----------------------------------------------------------------
