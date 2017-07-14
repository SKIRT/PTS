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

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class Elitism(object):

    """
    This class ...
    """

    def __init__(self, index, replaced, replacement, replaced_score, replacement_score, name=None, replacement_name=None):

        """
        This function ...
        """

        self.index = index
        self.replaced = replaced
        self.replacement = replacement
        self.replaced_score = replaced_score
        self.replacement_score = replacement_score
        self.name = name
        self.replacement_name = replacement_name

# -----------------------------------------------------------------
