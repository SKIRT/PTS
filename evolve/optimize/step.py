#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.step Contains the StepOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..genomes.list1d import G1DList
from ..genomes.list2d import G2DList
from ..genomes.binarystring1d import G1DBinaryString
from ..genomes.binarystring2d import G2DBinaryString
from ...core.tools.logging import log
from ..core.initializators import G1DListInitializatorReal, G1DListInitializatorInteger
from ..core.mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary, G1DListMutatorRealGaussian, G1DListMutatorRealRange
from ..core.engine import GeneticEngine, RawScoreCriteria
from ..core import constants
from ...core.basics.range import RealRange, IntegerRange
from ...core.tools import filesystem as fs
from ...core.tools.random import save_state, load_state
from ...core.basics.configuration import Configuration
from ..core.dbadapters import DBFileCSV, DBSQLite

# -----------------------------------------------------------------

class StepOptimizer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(StepOptimizer, self).__init__(config)

        # The genetic algorithm engine
        self.engine = None

        # The best individual
        self.best = None

        # The current population
        self.population = None

        # The scores of the past generation
        self.scores = None

        # The database
        self.database = None

        # The path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, dir_path):

        """
        This function ...
        :param dir_path:
        :return:
        """

        # Determine the path to the engine
        engine_path = fs.join(dir_path, "engine.pickle")

        # Load the engine
        engine = GeneticEngine.from_file(engine_path)

        # Determine the path to the random generator
        prng_path = fs.join(dir_path, "prng.pickle")

        # Load the random state
        load_state(prng_path)

        # Determine the path to the configuration
        config_path = fs.join(dir_path, "optimizer.cfg")

        # Load the configuration
        config = Configuration.from_file(config_path)

        # Create the optimizer instance
        optimizer = cls(config)

        # Set the engine
        optimizer.engine = engine

        # Set the database
        database_path = fs.join(dir_path, "database.csv")
        optimizer.database = DBFileCSV(filename=database_path)

        #database_path = fs.join(dir_path, "database.db")
        #optimizer.database = DBSQLite(dbname=database_path)
        # Opening is done down below

        # Set the path
        optimizer.path = dir_path

        # Return the optimizer
        return optimizer

    # -----------------------------------------------------------------

    @property
    def initialized(self):

        """
        This function ...
        :return:
        """

        return self.engine is not None

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
        if not self.initialized: self.initialize(**kwargs)

        # 3. Evolve
        else: self.evolve()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

        # Return the current population
        return self.population

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(StepOptimizer, self).setup(**kwargs)

        # Get the engine
        if "engine" in kwargs: self.engine = kwargs.pop("engine")

        # Get the scores
        if "scores" in kwargs: self.scores = kwargs.pop("scores")

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing ...")

        # 1. Initialize databse
        self.initialize_database()

        # 1. Initialize genome
        self.initialize_genome(**kwargs)

        # 2. Initialize engine
        self.initialize_engine(**kwargs)

        # 3. Initialize the evolution (the initial generation)
        self.initialize_evolution()

    # -----------------------------------------------------------------

    def initialize_database(self):

        """
        This function ...
        :return:
        """

        filepath = fs.join(fs.cwd(), "database.csv")

        # Create database adapter
        self.database = DBFileCSV(filename=filepath, identify="run_01", frequency=1, reset=True)

        #filepath = fs.join(fs.cwd(), "database.db")

        #self.database = DBSQLite(identify="ex1", resetDB=True)
        #self.database = DBSQLite(dbname=filepath, identify="ex1", resetDB=True, commit_freq=1, frequency=1)

        #print(self.engine.dbAdapter)

    # -----------------------------------------------------------------

    def initialize_genome(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the starting genome ...")

        # Get evaluator
        #evaluator = kwargs.pop("evaluator")

        # Get the genome
        genome = kwargs.pop("genome").clone() if "genome" in kwargs else None

        # Get initializator
        initializator = kwargs.pop("initializator") if "initializator" in kwargs else None

        # Get mutator
        mutator = kwargs.pop("mutator") if "mutator" in kwargs else None

        # Get crossover function
        crossover = kwargs.pop("crossover") if "crossover" in kwargs else None

        # Create the genome if necessary
        if genome is None:

            # List genome
            if self.config.genome_type == "list":

                # Create 1D genome
                if self.config.genome_dimension == 1: genome = G1DList(self.config.nparameters)
                elif self.config.genome_dimension == 2: genome = G2DList(self.config.nparameters, self.config.nparameters2)
                else: raise ValueError("Dimensions > 2 are not supported")

            # Binary string genome
            elif self.config.genome_type == "binary_string":

                # 1D or 2D
                if self.config.genome_dimension == 1: genome = G1DBinaryString(self.config.nparameters)
                elif self.config.genome_dimension == 2: genome = G2DBinaryString(self.config.nparameters, self.config.nparameters2)
                else: raise ValueError("Dimensions > 2 are not supported")

            # Invalid option
            else: raise ValueError("Genome type must be 'list' or 'binary_string")

        # Set basic properties
        if self.config.parameter_range is not None: genome.setParams(rangemin=self.config.parameter_range.min, rangemax=self.config.parameter_range.max)
        if self.config.best_raw_score is not None: genome.setParams(bestrawscore=self.config.best_raw_score)
        if self.config.round_decimal is not None: genome.setParams(rounddecimal=self.config.round_decimal)

        # Set initializator
        if initializator is not None: genome.initializator.set(initializator)
        else:
            if isinstance(self.config.parameter_range, IntegerRange): genome.initializator.set(G1DListInitializatorInteger)
            elif isinstance(self.config.parameter_range, RealRange): genome.initializator.set(G1DListInitializatorReal)
            else: raise ValueError("Invalid parameter range")

        # Set mutator
        if mutator is not None: genome.mutator.set(mutator)
        else:
            if isinstance(self.config.parameter_range, IntegerRange):
                if self.config.mutation_method == "range": genome.mutator.set(G1DListMutatorIntegerRange)
                elif self.config.mutation_method == "gaussian": genome.mutator.set(G1DListMutatorIntegerGaussian)
                elif self.config.mutation_method == "binary": genome.mutator.set(G1DListMutatorIntegerBinary)
                else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not recognized")
            elif isinstance(self.config.parameter_range, RealRange):
                if self.config.mutation_method == "range": genome.mutator.set(G1DListMutatorRealRange)
                elif self.config.mutation_method == "gaussian": genome.mutator.set(G1DListMutatorRealGaussian)
                else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not valid for genome of real values")
            else: raise ValueError("Invalid parameter range")

        # Set crossover
        if crossover is not None: genome.crossover.set(crossover)

        # Set the evaluator
        #genome.evaluator.set(evaluator)

        # Return the genome
        return genome

    # -----------------------------------------------------------------

    def initialize_engine(self, genome, **kwargs):

        """
        This function ...
        :param genome:
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Initializing the genetic engine ...")

        # Create the genetic engine
        self.engine = GeneticEngine(genome)

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

        # Get the callback function
        callback = kwargs.pop("callback") if "callback" in kwargs else None

        # Get the database adapter function
        #adapter = kwargs.pop("adapter") if "adapter" in kwargs else None

        # Get kwargs for the different functions
        evaluator_kwargs = kwargs.pop("evaluator_kwargs") if "evaluator_kwargs" in kwargs else None
        initializator_kwargs = kwargs.pop("initializator_kwargs") if "initializator_kwargs" in kwargs else None
        mutator_kwargs = kwargs.pop("mutator_kwargs") if "mutator_kwargs" in kwargs else None
        crossover_kwargs = kwargs.pop("crossover_kwargs") if "crossover_kwargs" in kwargs else None
        callback_kwargs = kwargs.pop("callback_kwargs") if "callback_kwargs" in kwargs else None

        # Set the callback function
        if callback is not None: self.engine.stepCallback.set(callback)

        # Set the database adapter
        #if adapter is not None:
        #    adapter.open(self.engine)
        #    self.engine.setDBAdapter(adapter)

        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Set kwargs
        if evaluator_kwargs is not None: self.engine.set_kwargs("evaluator", evaluator_kwargs)
        if initializator_kwargs is not None: self.engine.set_kwargs("initializator", initializator_kwargs)
        if mutator_kwargs is not None: self.engine.set_kwargs("mutator", mutator_kwargs)
        if crossover_kwargs is not None: self.engine.set_kwargs("crossover", crossover_kwargs)
        if callback_kwargs is not None: self.engine.set_kwargs("callback", callback_kwargs)

    # -----------------------------------------------------------------

    def initialize_evolution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initialize the evolution ...")

        # Initialize the genetic algorithm
        self.engine.initialize()

        # Get the initial population
        self.population = self.engine.get_population()

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolve ...")

        # Set the database adapter again
        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Make sure we don't stop unexpectedly, always increment the number of generations we want
        # (it doesn't matter what the value is)
        self.engine.setGenerations(self.engine.getGenerations() + 1)

        # Set the scores from the previous generation
        self.set_scores()

        # Generate the new models
        self.generate_new_models()

    # -----------------------------------------------------------------

    def set_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting scores from previous generation ...")

        # Check if the scores are set
        if self.scores is None: raise ValueError("The scores are not set")

        # Set the scores
        self.engine.set_scores(self.scores)

    # -----------------------------------------------------------------

    def generate_new_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new population ...")

        # Generate the new population
        self.engine.generate_new_population()

        # Get the new population
        self.population = self.engine.new_population

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # INform the user
        log.info("Writing ...")

        # Write the new generation
        self.write_population()

    # -----------------------------------------------------------------

    def write_population(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the population genomes ...")

        # Loop over the individuals
        for individual in self.population:

            pass

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check if path is present
        if self.path is None: raise ValueError("Path is not set")

        # Save to the current path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, dir_path):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving to '" + dir_path + "' ...")

        # Save the engine
        self.save_engine(dir_path)

        # Save the state of the prng
        self.save_state(dir_path)

        # Save the configuration
        self.save_config(dir_path)

        # Save the database
        self.save_database(dir_path)

    # -----------------------------------------------------------------

    def save_engine(self, dir_path):

        """
        This function ...
        :return:
        """

        # Determine the path to the engine
        engine_path = fs.join(dir_path, "engine.pickle")

        # Save the engine
        self.engine.saveto(engine_path)

    # -----------------------------------------------------------------

    def save_state(self, dir_path):

        """
        This function ...
        :return:
        """

        # Determine the path to the prng
        prng_path = fs.join(dir_path, "prng.pickle")

        # Save the prng state
        save_state(prng_path)

    # -----------------------------------------------------------------

    def save_config(self, dir_path):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = fs.join(dir_path, "optimizer.cfg")

        # Save the configuration
        self.config.saveto(path)

    # -----------------------------------------------------------------

    def save_database(self, dir_path):

        """
        This function ...
        :param dir_path:
        :return:
        """

        #pass
        self.database.close()

# -----------------------------------------------------------------
