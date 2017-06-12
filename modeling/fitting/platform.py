#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.platform Contains the GenerationPlatform class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from .reproduction import ReproductionEvent

# -----------------------------------------------------------------

class GenerationPlatform(object):
    
    """
    This class...
    """

    def __init__(self, generation):

        """
        The constructor ...
        :param generation:
        :return:
        """

        self.generation = generation

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, generation_path):

        """
        This function ...
        :param generation_path:
        :return:
        """

        # Load generation
        generation = Generation.from_path(generation_path)

        # Create and return platform
        return cls(generation)

    # -----------------------------------------------------------------

    def make_genome(self, genes):

        """
        Tihs function ...
        :param genes:
        :return:
        """

        return self.generation.genome_class(genes=genes)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.generation.fitting_run

    # -----------------------------------------------------------------

    def genome_to_parameters(self, genome):

        """
        This function ...
        :param genome:
        :return:
        """

        # genome, fitting_run, minima, maxima, nbits, parameter_scales, gray=False
        return get_parameter_values_from_genome(genome, self.fitting_run, self.generation.parameter_minima_scalar, self.generation.parameter_maxima_scalar, self.generation.nbits_list, self.generation.parameter_scales, self.generation.gray_code)

    # -----------------------------------------------------------------

    def genes_to_parameters(self, genes):

        """
        This function ...
        :param genes:
        :return:
        """

        genome = self.make_genome(genes)
        return self.genome_to_parameters(genome)

    # -----------------------------------------------------------------

    def crossover(self, mother, father, details=None):

        """
        This function ...
        :param mother:
        :param father:
        :param details:
        :return:
        """

        sister, brother = self.generation.crossover_function(None, mom=mother, dad=father, details=details, count=2)
        return sister, brother

    # -----------------------------------------------------------------

    def crossover_origins(self, size, details):

        """
        This function ...
        :param size:
        :param details:
        :return:
        """

        sister_origins, brother_origins = self.generation.crossover_origins_function(size, details)
        return sister_origins, brother_origins

    # -----------------------------------------------------------------

    @property
    def nreproductions(self):

        """
        This function ...
        :return:
        """

        return len(self.generation.crossover_table)

    # -----------------------------------------------------------------

    @property
    def ncrossovers(self):

        """
        This function ...
        :return:
        """

        return self.generation.crossover_table.ncrossovers

    # -----------------------------------------------------------------

    @property
    def reproductions(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nreproductions): yield self.get_reproduction(index)

    # -----------------------------------------------------------------

    def get_reproduction(self, index):

        """
        This fucntion ...
        :param index:
        :return:
        """

        mother_name, father_name = self.generation.crossover_table.get_parents(index)
        sister_name, brother_name = self.generation.crossover_table.get_children(index)

        # Load parent and new populations
        parents = self.generation.parents
        newborns = self.generation.newborns

        mother = self.make_genome(parents[mother_name])
        father = self.make_genome(parents[father_name])

        sister = self.make_genome(newborns[sister_name])
        brother = self.make_genome(newborns[brother_name])

        crossover = False

        # Crossover happened
        if self.generation.crossover_table.is_crossover(index):

            crossover = True

            # Get crossover method and details
            #method = self.generation.crossover_table.get_crossover_method(index)
            details = self.generation.crossover_table.get_crossover_details(index)

            # Create crossover genomes
            initial_sister, initial_brother = self.crossover(mother, father, details)

            # Get the origins
            sister_origins, brother_origins = self.crossover_origins(len(mother), details)

        # Just cloned
        else: initial_sister, initial_brother, sister_origins, brother_origins = mother, father, None, None

        # Return
        return ReproductionEvent(index, mother, father, initial_sister, initial_brother, sister, brother, crossover, sister_origins, brother_origins)

# -----------------------------------------------------------------
