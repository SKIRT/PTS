#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize Contains the Optimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable
from .genomes.list1d import G1DList
from ..core.tools.logging import log
from .initializators import G1DListInitializatorReal, G1DListInitializatorInteger
from .mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary, G1DListMutatorRealGaussian, G1DListMutatorRealRange
from .engine import GeneticEngine, RawScoreCriteria
from . import constants
from ..core.basics.range import RealRange, IntegerRange
from ..core.tools import formatting as fmt

# -----------------------------------------------------------------

class Optimizer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(Optimizer, self).__init__(config)

        # The evaluating function
        self.evaluator = None

        # The initializator function
        self.initializator = None

        # The starting genome
        self.genome = None

        # The genetic algorithm engine
        self.engine = None

        # The best individual
        self.best = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Initialize
        self.initialize()

        # 3. Evolve
        self.evolve()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Optimizer, self).setup(**kwargs)

        # Get evaluator
        if "evaluator" in kwargs: self.evaluator = kwargs.pop("evaluator")

        # Get initializator
        if "initializator" in kwargs: self.initializator = kwargs.pop("initializator")

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing ...")

        # 1. Initialize genome
        self.initialize_genome()

        # 2. Initialize engine
        self.initialize_engine()

    # -----------------------------------------------------------------

    def initialize_genome(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the starting genome ...")

        # Create the genome
        self.genome = G1DList(self.config.nparameters)

        # Set properties
        self.genome.setParams(rangemin=self.config.parameter_range.min, rangemax=self.config.parameter_range.max)
        if self.config.best_raw_score is not None: self.genome.setParams(bestrawscore=self.config.best_raw_score)
        if self.config.round_decimal is not None: self.genome.setParams(rounddecimal=self.config.round_decimal)

        # Set initializator
        if self.initializator is not None: self.genome.initializator.set(self.initializator)
        else:
            if isinstance(self.config.parameter_range, IntegerRange): self.genome.initializator.set(G1DListInitializatorInteger)
            elif isinstance(self.config.parameter_range, RealRange): self.genome.initializator.set(G1DListInitializatorReal)
            else: raise ValueError("Invalid parameter range")

        # Set mutator
        if isinstance(self.config.parameter_range, IntegerRange):
            if self.config.mutation_method == "range": self.genome.mutator.set(G1DListMutatorIntegerRange)
            elif self.config.mutation_method == "gaussian": self.genome.mutator.set(G1DListMutatorIntegerGaussian)
            elif self.config.mutation_method == "binary": self.genome.mutator.set(G1DListMutatorIntegerBinary)
            else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not recognized")
        elif isinstance(self.config.parameter_range, RealRange):
            if self.config.mutation_method == "range": self.genome.mutator.set(G1DListMutatorRealRange)
            elif self.config.mutation_method == "gaussian": self.genome.mutator.set(G1DListMutatorRealGaussian)
            else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not valid for genome of real values")
        else: raise ValueError("Invalid parameter range")

        # Set the evaluator
        self.genome.evaluator.set(self.evaluator)

    # -----------------------------------------------------------------

    def initialize_engine(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the genetic engine ...")

        # Create the genetic engine
        self.engine = GeneticEngine(self.genome)

        # Set options
        self.engine.terminationCriteria.set(RawScoreCriteria)
        self.engine.setMinimax(constants.minimaxType[self.config.min_or_max])
        self.engine.setGenerations(self.config.ngenerations)
        if self.config.crossover_rate is not None: self.engine.setCrossoverRate(self.config.crossover_rate)
        self.engine.setPopulationSize(self.config.nindividuals)
        self.engine.setMutationRate(self.config.mutation_rate)

        # Set elitism options
        self.engine.setElitism(self.config.elitism)
        self.engine.setElitismReplacement(self.config.nelite_individuals)

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolving ...")

        # Let evolve
        self.engine.evolve(freq_stats=self.config.stats_freq, progress_bar=self.config.progress_bar)

        # Get the best individual
        self.best = self.engine.bestIndividual()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the best individual
        show_best(self.best)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------

def show_best(best):

    """
    This function ...
    :return:
    """

    # Get properties of best individual
    score = best.getRawScore()
    fitness = best.getFitnessScore()
    parameters = best.internalParams
    genome = best.genomeList

    print("")
    print(fmt.green + fmt.bold + "Best individual:" + fmt.reset)
    print("")

    print(fmt.yellow + "  Genome:" + fmt.reset)
    print("")
    for gene in genome: print("    " + str(gene))
    print("")

    print(fmt.yellow + "  Score: " + fmt.reset + str(score))

    print("")

    print(fmt.yellow + "  Fitness: " + fmt.reset + str(fitness))

    print("")

# -----------------------------------------------------------------
