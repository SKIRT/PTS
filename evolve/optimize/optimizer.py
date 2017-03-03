#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.optimizer Contains the Optimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

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
from ...core.tools import formatting as fmt
from ...core.tools import stringify
from ..core.adapters import DBFileCSV, DBSQLite
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class Optimizer(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(Optimizer, self).__init__(config, interactive)

        # The database adapter
        self.database = None

        # The statistics table adapter
        self.statistics = None

        # The intial genome
        self.initial_genome = None

        # The genetic algorithm engine
        self.engine = None

        # The best individual
        self.best = None

        # The generations plotter
        self.generations_plotter = None

        # The plot path
        self.plot_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Optimizer, self).setup(**kwargs)

        # Set the plot path
        self.plot_path = fs.create_directory_in(self.config.path, "plot")

    # -----------------------------------------------------------------

    def initialize_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the statistics table ...")

        # Determine the file path
        if self.config.writing.statistics_path is not None: filepath = fs.absolute_or_in(self.config.writing.statistics_path, self.output_path)
        else: filepath = self.output_path_file("statistics.csv")

        # Check the file
        if fs.is_file(filepath):
            log.debug("The statistics file is already present. Adding a new run to the file.")
            reset = False
        else:
            log.debug("Creating a new statistics file ...")
            reset = True

        # Create the adapter
        self.statistics = DBFileCSV(filename=filepath, frequency=self.config.statistics_frequency, reset=reset, identify=self.config.run_id)

    # -----------------------------------------------------------------

    def initialize_database(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the database ...")

        # Determine the file path
        if self.config.writing.database_path is not None: filepath = fs.absolute_or_in(self.config.writing.database_path, self.output_path)
        else: filepath = self.output_path_file("database.db")

        # Check the file
        if fs.is_file(filepath):
            log.debug("The database file is already present. Adding a new run to the database.")
            reset = False
        else:
            log.debug("Creating a new database file ...")
            reset = True

        # Create the database adapter
        self.database = DBSQLite(dbname=filepath, identify=self.config.run_id, resetDB=reset, commit_freq=self.config.database_frequency, frequency=self.config.database_frequency)

    # -----------------------------------------------------------------

    def initialize_genome(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Initializing the starting genome ...")

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

        # Set the initial genome
        self.initial_genome = genome

    # -----------------------------------------------------------------

    def initialize_engine(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Initializing the genetic engine ...")

        # Create the genetic engine
        self.engine = GeneticEngine(self.initial_genome)

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
        #callback = kwargs.pop("callback") if "callback" in kwargs else None

        # Get kwargs for the different functions
        evaluator_kwargs = kwargs.pop("evaluator_kwargs") if "evaluator_kwargs" in kwargs else None
        initializator_kwargs = kwargs.pop("initializator_kwargs") if "initializator_kwargs" in kwargs else None
        mutator_kwargs = kwargs.pop("mutator_kwargs") if "mutator_kwargs" in kwargs else None
        crossover_kwargs = kwargs.pop("crossover_kwargs") if "crossover_kwargs" in kwargs else None
        #callback_kwargs = kwargs.pop("callback_kwargs") if "callback_kwargs" in kwargs else None

        # Set the callback function
        #if callback is not None: self.engine.stepCallback.set(callback)

        # Set generations plotter
        self.generations_plotter = kwargs.pop("generations_plotter", None)
        self.generations_plotter.output_path = self.plot_path
        if self.generations_plotter is not None: self.engine.stepCallback.set(self.generations_plotter.add_generation)

        # Set the database adapter
        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Set the adapter for the statistics table
        self.statistics.open(self.engine)
        self.engine.setDBAdapter(self.statistics)

        # Set kwargs
        if evaluator_kwargs is not None: self.engine.set_kwargs("evaluator", evaluator_kwargs)
        if initializator_kwargs is not None: self.engine.set_kwargs("initializator", initializator_kwargs)
        if mutator_kwargs is not None: self.engine.set_kwargs("mutator", mutator_kwargs)
        if crossover_kwargs is not None: self.engine.set_kwargs("crossover", crossover_kwargs)
        #if callback_kwargs is not None: self.engine.set_kwargs("callback", callback_kwargs)

    # -----------------------------------------------------------------

    @abstractmethod
    def evolve(self):

        """
        This fucntion ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def show(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def write(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def plot(self):

        """
        This function ...
        :return:
        """

        pass

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
    print(stringify.stringify_list_fancy([gene for gene in genome], 100, ", ", "    ")[1])
    print("")

    print(fmt.yellow + "  Score: " + fmt.reset + str(score))

    print("")

    print(fmt.yellow + "  Fitness: " + fmt.reset + str(fitness))

    print("")

# -----------------------------------------------------------------
