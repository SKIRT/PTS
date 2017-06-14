#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.recurrence Contains the Recurrence class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .generation import Generation
from .evaluate import get_parameter_values_from_genome

# -----------------------------------------------------------------

class Recurrence(object):

    """
    This class ...
    """

    def __init__(self, index, individual, generation, original, score):

        """
        This function ...
        :param index:
        :param individual:
        :param generation:
        :param original:
        :param score:
        """

        self.index = index
        self.individual = individual
        self.generation = generation
        self.original = original
        self.score = score

# -----------------------------------------------------------------
