#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.initial Contains the InitialModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....evolve.engine import GAEngine, RawScoreCriteria
from ....evolve.genomes.list1d import G1DList
from ....evolve import mutators
from ....evolve import initializators
from ....evolve import constants
from .generator import ModelGenerator
from ....core.tools.random import save_state

# -----------------------------------------------------------------

class InitialModelGenerator(ModelGenerator):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(InitialModelGenerator, self).__init__()

        # The genetic algorithm engine
        self.engine = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(InitialModelGenerator, self).setup()

        # Create the first genome
        genome = G1DList(self.nparameters)

        # Set genome options
        genome.setParams(minima=self.parameter_minima, maxima=self.parameter_maxima, bestrawscore=0.00, rounddecimal=2)
        genome.initializator.set(initializators.HeterogeneousListInitializerReal)
        # genome.mutator.set(mutators.HeterogeneousListMutatorRealRange)
        genome.mutator.set(mutators.HeterogeneousListMutatorRealGaussian)

        # Create the genetic algorithm engine
        self.engine = GAEngine(genome)

        # Set options for the engine
        self.engine.terminationCriteria.set(RawScoreCriteria)
        self.engine.setMinimax(constants.minimaxType["minimize"])
        self.engine.setGenerations(5)
        self.engine.setCrossoverRate(0.5)
        self.engine.setPopulationSize(100)
        self.engine.setMutationRate(0.5)

        # Initialize the genetic algorithm
        self.engine.initialize()

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the initial population of models ...")

        # Get the initial population
        population = self.engine.get_population()

        # Loop over the individuals of the population
        for individual in population:

            # Loop over all the genes (parameters)
            for i in range(len(individual)):

                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.parameter_labels_order[i]].append(value)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the genetic algorithm engine
        self.write_engine()

        # Write the state of the random number generator
        self.write_prng()

    # -----------------------------------------------------------------

    def write_engine(self):

        """
        This function ...
        :return:
        """

        # Save the genetic algorithm
        self.engine.saveto(path)

    # -----------------------------------------------------------------

    def write_prng(self):

        """
        This function ...
        :return:
        """

        # Save the state of the random generator
        save_state(random_path)

# -----------------------------------------------------------------
