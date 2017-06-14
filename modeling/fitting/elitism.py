#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.elitism Contains the Elitism class.

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

class Elitism(object):

    """
    This class ...
    """

    def __init__(self, index, mother, father, initial_sister, initial_brother, sister, brother, crossover, sister_origins, brother_origins):

        """
        This function ...
        """

        self.index = index
        self.mother = mother
        self.father = father
        self.initial_sister = initial_sister
        self.initial_brother = initial_brother
        self.sister = sister
        self.brother = brother
        self.crossover = crossover
        self.sister_origins = sister_origins
        self.brother_origins = brother_origins

# -----------------------------------------------------------------
