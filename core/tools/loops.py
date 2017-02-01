#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.loops Functions related to loops.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

def repeat(target, ntimes):

    """
    This function ...
    :param target:
    :param ntimes:
    :return:
    """

    for _ in range(ntimes): target()

# -----------------------------------------------------------------
