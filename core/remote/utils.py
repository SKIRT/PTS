#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.utils Contains remote utilities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class HostDownException(Exception):

    """
    This exception should be raised when connection to a host is not possible because it is temporarily down
    """

# -----------------------------------------------------------------

class DetachedCalculation(Exception):
        
   """
   """

# -----------------------------------------------------------------
