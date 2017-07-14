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

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class Recurrence(object):

    """
    This class ...
    """

    def __init__(self, index, individual, generation, original, score, parameters, original_parameters, name=None, original_name=None):

        """
        This function ...
        :param index:
        :param individual:
        :param generation:
        :param original:
        :param score:
        :param parameters:
        :param original_parameters:
        """

        self.index = index
        self.individual = individual
        self.generation = generation
        self.original = original
        self.score = score
        self.parameters = parameters
        self.original_parameters = original_parameters
        self.name = name
        self.original_name = original_name

    # -----------------------------------------------------------------

    @property
    def genome_size(self):

        """
        This function ...
        :return:
        """

        return len(self.individual)

    # -----------------------------------------------------------------

    @lazyproperty
    def differences(self):

        """
        This function ...
        :return:
        """

        return [self.individual.genes[i] != self.original.genes[i] for i in range(self.genome_size)]

    # -----------------------------------------------------------------

    @property
    def ndifferences(self):

        """
        This function ...
        :return:
        """

        return sum(self.differences)

# -----------------------------------------------------------------
