#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.reproduction Contains the ReproductionEvent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class ReproductionEvent(object):

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

    @lazyproperty
    def sister_mutations(self):

        """
        This function ...
        :return:
        """

        return [self.initial_sister.genes[i] != self.sister.genes[i] for i in range(len(self.initial_sister))]

    # -----------------------------------------------------------------

    @lazyproperty
    def brother_mutations(self):

        """
        This function ...
        :return:
        """

        return [self.initial_brother.genes[i] != self.brother.genes[i] for i in range(len(self.initial_brother))]

    # -----------------------------------------------------------------

    @lazyproperty
    def nmutations_sister(self):

        """
        This function ...
        :return: 
        """

        return sum(self.sister_mutations)

    # -----------------------------------------------------------------

    @lazyproperty
    def nmutations_brother(self):

        """
        This function ...
        :return:
        """

        return sum(self.brother_mutations)

    # -----------------------------------------------------------------

    @lazyproperty
    def relative_nmutations_sister(self):

        """
        This function ...
        :return:
        """

        return float(self.nmutations_sister) / len(self.sister)

    # -----------------------------------------------------------------

    @lazyproperty
    def relative_nmutations_brother(self):

        """
        This function ...
        :return:
        """

        return float(self.nmutations_brother) / len(self.brother)

# -----------------------------------------------------------------
