#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.geneticexplorer Contains the GeneticParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .explorer import ParameterExplorer
from ...core.tools.logging import log
from ...evolve.engine import GAEngine, RawScoreCriteria
from ...evolve.genomes.list1d import G1DList
from ...evolve import mutators
from ...evolve import initializators
from ...evolve import constants

# -----------------------------------------------------------------

class GeneticParameterExplorer(ParameterExplorer):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(GeneticParameterExplorer, self).__init__(config)

        # The genetic algorithm engine
        self.engine = None

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parameter ranges ...")

        # Get the current values in the ski file prepared by InputInitializer
        young_luminosity_guess, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
        ionizing_luminosity_guess, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
        dust_mass_guess = self.ski.get_dust_component_mass(0)

        # Inform the user
        log.info("Determining every possible combination of parameter values ...")

        # Determine the ranges
        #fuv_young_range = self.young_luminosity_range(young_luminosity_guess)
        #fuv_ionizing_range = self.ionizing_luminosity_range(ionizing_luminosity_guess)
        #dust_mass_range = self.dust_mass_range(dust_mass_guess)

        # Create the first genome
        genome = G1DList(3)
        genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
        genome.initializator.set(initializators.G1DListInitializatorReal)
        genome.mutator.set(mutators.G1DListMutatorRealGaussian)

        # Create the genetic algorithm engine
        self.engine = GAEngine(genome)

        self.engine.terminationCriteria.set(RawScoreCriteria)
        self.engine.setMinimax(constants.minimaxType["minimize"])
        self.engine.setGenerations(5)
        self.engine.setCrossoverRate(0.5)
        self.engine.setPopulationSize(100)
        self.engine.setMutationRate(0.5)

        # Initialize the genetic algorithm
        self.engine.initialize()

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

# -----------------------------------------------------------------
