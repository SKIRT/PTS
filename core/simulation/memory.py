#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.memory Contains the MemoryRequirement class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class MemoryRequirement(object):

    """
    This class ...
    """

    def __init__(self, serial, parallel):

        """
        This function ...
        :param serial:
        :param parallel:
        """

        self.serial = serial
        self.parallel = parallel

# -----------------------------------------------------------------
