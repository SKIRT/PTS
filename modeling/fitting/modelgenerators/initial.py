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

        # The dictionary with the list of the model parameters
        self.parameters = dict()

        # The genetic algorithm engine
        self.engine = None

    # -----------------------------------------------------------------

    def run(self, ranges):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(ranges)

        # Generate the model parameters
        self.generate()

    # -----------------------------------------------------------------

    def setup(self, ranges):

        """
        This function ...
        :return:
        """

        # Set minima and maxima for the different genes (model parameters)
        minima = [ranges["FUV young"][0], ranges["FUV ionizing"][0], ranges["Dust mass"][0]]
        maxima = [ranges["FUV young"][1], ranges["FUV ionizing"][1], ranges["Dust mass"][1]]

        # Create the first genome
        genome = G1DList(3)

        # genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
        # genome.initializator.set(initializators.G1DListInitializatorReal)
        # genome.mutator.set(mutators.G1DListMutatorRealGaussian)

        genome.setParams(minima=minima, maxima=maxima, bestrawscore=0.00, rounddecimal=2)
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

        # Get the initial population
        population = self.engine.get_population()

        # Loop over the individuals of the population
        for individual in population:

            young_luminosity = individual[0]
            ionizing_luminosity = individual[1]
            dust_mass = individual[2]

            # Add the parameter values to the dictionary
            self.parameters["FUV young"].append(young_luminosity)
            self.parameters["FUV ionizing"].append(ionizing_luminosity)
            self.parameters["Dust mass"].append(dust_mass)

# -----------------------------------------------------------------
