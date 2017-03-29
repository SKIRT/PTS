#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.utils Contains utilities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class UserIntervention(Exception):

    """
    This exception should be called when user intervention / inspection is required and therefore
    the execution should be halted or aborted.
    """

# -----------------------------------------------------------------
