#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.welcome Contains a function that prints a welcome message.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection

# -----------------------------------------------------------------

def welcome(width=50, prefix="    "):

    """
    This function ...
    :param width:
    :param prefix:
    :return:
    """

    fmt.print_empty()
    fmt.print_filled("#", prefix=prefix, length=width)
    fmt.print_border("#", prefix=prefix, length=width)
    fmt.print_centered_around_border("Welcome to PTS/Modeling, " + introspection.username(), "#", prefix=prefix, length=width)
    fmt.print_centered_around_border("Author: Sam Verstocken", "#", prefix=prefix, length=width)
    fmt.print_centered_around_border("PTS version: " + introspection.pts_version(), "#", prefix=prefix, length=width)
    fmt.print_centered_around_border("Last update: " + introspection.pts_update_date(), "#", prefix=prefix, length=width)
    fmt.print_border("#", prefix=prefix, length=width)
    fmt.print_centered_around_border("© Astronomical Observatory, Ghent University", "#", prefix=prefix, length=width)
    fmt.print_border("#", prefix=prefix, length=width)
    fmt.print_filled("#", prefix=prefix, length=width)
    fmt.print_empty()

# -----------------------------------------------------------------
