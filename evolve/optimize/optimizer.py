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

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..genomes.list1d import G1DList
from ..genomes.list2d import G2DList
from ..genomes.binarystring1d import G1DBinaryString
from ..genomes.binarystring2d import G2DBinaryString
from ...core.tools.logging import log
from ..core.initializators import G1DListInitializatorReal, G1DListInitializatorInteger, HeterogeneousListInitializerReal, HeterogeneousListInitializerInteger
from ..core.crossovers import G1DListCrossoverSinglePoint, G1DListCrossoverTwoPoint, G1DListCrossoverUniform, G1DListCrossoverOX
from ..core.crossovers import G1DListCrossoverEdge, G1DListCrossoverCutCrossfill, G1DListCrossoverRealSBX
from ..core.crossovers import G2DListCrossoverUniform, G2DListCrossoverSingleVPoint, G2DListCrossoverSingleHPoint
from ..core.mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary
from ..core.mutators import G1DListMutatorRealGaussian, G1DListMutatorRealRange
from ..core.mutators import HeterogeneousListMutatorRealRange, HeterogeneousListMutatorRealGaussian
from ..core.mutators import HeterogeneousListMutatorIntegerRange, HeterogeneousListMutatorIntegerGaussian
from ..core.engine import GeneticEngine, RawScoreCriteria
from ..core import constants
from ...core.basics.range import RealRange, IntegerRange
from ...core.tools import formatting as fmt
from ...core.tools import stringify
from ..core.adapters import DBFileCSV, DBSQLite
from ...core.tools import filesystem as fs
from ...core.tools import types

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

        # The parameter minima and maxima (for heterogeneous genomes)
        self.parameter_minima = None
        self.parameter_maxima = None

        # The parameter centers and sigmas (for heterogeneous genomes with Gaussian mutator and/or initializer)
        self.parameter_centers = None
        self.parameter_sigmas = None

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_base_type(self):

        """
        This function ...
        :return:
        """

        # Parameter range is defined, not hetereogeneous
        if self.config.parameter_range is not None:

            if isinstance(self.config.parameter_range, IntegerRange): return "integer"
            elif isinstance(self.config.parameter_range, RealRange): return "real"
            else: raise ValueError("Invalid parameter range")

        # Parameter minima are defined, heterogeneous
        elif self.parameter_minima is not None:

            if types.is_integer_type(self.parameter_minima[0]): return "integer"
            elif types.is_real_type(self.parameter_minima[0]): return "real"
            else: raise ValueError("Invalid parameter minima")

        # Parameter maxima are defined, heterogeneous
        elif self.parameter_maxima is not None:

            if types.is_integer_type(self.parameter_maxima[0]): return "integer"
            elif types.is_real_type(self.parameter_maxima[0]): return "real"
            else: raise ValueError("Invalid parameter maxima")

        # Parameter type is defined in the configuration
        elif self.config.parameter_type is not None: return self.config.parameter_type

        # No clue
        else: raise ValueError("Parameter type cannot be determined or is invalid")

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

        # Set the parameter minima and maxima
        if "minima" in kwargs: self.parameter_minima = kwargs.pop("minima")
        if "maxima" in kwargs: self.parameter_maxima = kwargs.pop("maxima")

        # Set the parameter centers and sigmas
        if "centers" in kwargs: self.parameter_centers = kwargs.pop("centers")
        if "sigmas" in kwargs: self.parameter_sigmas = kwargs.pop("sigmas")

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
        if "genome" in kwargs: self.initial_genome = kwargs.pop("genome").clone()

        # Create the genome if necessary
        if self.initial_genome is None: self.create_genome()

        # Set genome properties
        self.set_genome_properties()

        # Set initializator
        self.set_genome_initializator(kwargs)

        # Set mutator
        self.set_genome_mutator(kwargs)

        # Set crossover
        self.set_genome_crossover(kwargs)

    # -----------------------------------------------------------------

    def create_genome(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the genome ...")

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

        # Set the genome
        self.initial_genome = genome

    # -----------------------------------------------------------------

    def set_genome_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting genome properties ...")

        # Debugging
        log.debug("Setting basic properties ...")

        # Set basic properties
        if self.config.parameter_range is not None: self.initial_genome.setParams(rangemin=self.config.parameter_range.min, rangemax=self.config.parameter_range.max)
        if self.config.best_raw_score is not None: self.initial_genome.setParams(bestrawscore=self.config.best_raw_score)
        if self.config.round_decimal is not None: self.initial_genome.setParams(rounddecimal=self.config.round_decimal)

        # Debugging
        log.debug("Setting parameter minima and maxima ...")

        # Set minima or maxima for heterogeneous genome lists
        if self.parameter_minima is not None: self.initial_genome.setParams(minima=self.parameter_minima)
        if self.parameter_maxima is not None: self.initial_genome.setParams(maxima=self.parameter_maxima)

        # Debugging
        log.debug("Setting parameter centers and sigmas ...")

        # Set parameter centers and sigmas
        if self.parameter_centers is not None: self.initial_genome.setParams(centers=self.parameter_centers)
        if self.parameter_sigmas is not None: self.initial_genome.setParams(sigmas=self.parameter_sigmas)

    # -----------------------------------------------------------------

    def set_genome_initializator(self, kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the initializator ...")

        # Get initializator
        initializator = kwargs.pop("initializator") if "initializator" in kwargs else None

        # Set initializator
        if initializator is not None: self.initial_genome.initializator.set(initializator)
        else: self.initial_genome.initializator.set(self.get_initializator())

    # -----------------------------------------------------------------

    def set_genome_mutator(self, kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the mutator ...")

        # Get mutator
        mutator = kwargs.pop("mutator") if "mutator" in kwargs else None

        # Set mutator
        if mutator is not None: self.initial_genome.mutator.set(mutator)
        else: self.initial_genome.mutator.set(self.get_mutator())

    # -----------------------------------------------------------------

    def set_genome_crossover(self, kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the crossover ...")

        # Get crossover function
        crossover = kwargs.pop("crossover") if "crossover" in kwargs else None

        # Set crossover
        if crossover is not None: self.initial_genome.crossover.set(crossover)
        else: self.initial_genome.crossover.set(self.get_crossover())

    # -----------------------------------------------------------------

    @property
    def genome_dimension(self):

        """
        This function ...
        :return:
        """

        if self.initial_genome is not None: return self.initial_genome.dimension
        else: return self.config.genome_dimension

    # -----------------------------------------------------------------

    @property
    def genome_size(self):

        """
        This function ...
        :return:
        """

        # Value for 1D, tuple for 2D

        if self.initial_genome is not None: return self.initial_genome.getSize()
        elif self.config.nparameters2 is None: return self.config.nparameters
        else: return (self.config.nparameters, self.config.nparameters2)

    # -----------------------------------------------------------------

    def get_initializator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the initializator type ...")

        # Integer type
        if self.parameter_base_type == "integer":

            if self.config.heterogeneous: return HeterogeneousListInitializerInteger
            else: return G1DListInitializatorInteger

        # Real type
        elif self.parameter_base_type == "real":

            if self.config.heterogeneous: return HeterogeneousListInitializerReal
            else: return G1DListInitializatorReal

        # Invalid
        else: raise ValueError("Invalid parameter type")

    # -----------------------------------------------------------------

    def get_mutator(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Getting the mutator type ...")

        # Integer type
        if self.parameter_base_type == "integer":

            # Range-based mutator
            if self.config.mutation_method == "range":

                # Choose class
                if self.config.heterogeneous: return HeterogeneousListMutatorIntegerRange
                else: return G1DListMutatorIntegerRange

            # Gaussian mutator
            elif self.config.mutation_method == "gaussian":

                # Choose class
                if self.config.heterogeneous: return HeterogeneousListMutatorIntegerGaussian
                else: return G1DListMutatorIntegerGaussian

            # Binary mutator
            elif self.config.mutation_method == "binary":

                # Choose class
                if self.config.heterogeneous: raise ValueError("Cannot use binary mutation on heterogeneous genomes")
                else: return G1DListMutatorIntegerBinary

            # Invalid
            else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not recognized")

        # Real type
        elif self.parameter_base_type == "real":

            # Range-based mutator
            if self.config.mutation_method == "range":

                # Choose class
                if self.config.heterogeneous: return HeterogeneousListMutatorRealRange
                else: return G1DListMutatorRealRange

            # Gaussian mutator
            elif self.config.mutation_method == "gaussian":

                # Choose class
                if self.config.heterogeneous: return HeterogeneousListMutatorRealGaussian
                else: return G1DListMutatorRealGaussian

            # Invalid
            else: raise ValueError("Mutation method '" + self.config.mutation_method + "' not valid for genome of real values")

        # Invalid
        else: raise ValueError("Invalid parameter type")

    # -----------------------------------------------------------------

    def get_crossover(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the crossover type ...")

        # 1D genome
        if self.config.genome_dimension == 1:

            # Check whether the genome size is larger than one. In that case, always use the uniform crossover
            if self.genome_size == 1:

                # Give a warning that we are not using the user specificied option
                log.warning("Uniform crossover will be used because the genome size is only one")
                return G1DListCrossoverUniform

            # Check whether the genome size is not zero
            elif self.genome_size == 0: raise ValueError("The genome size cannot be zero")

            # Single-point crossover
            if self.config.crossover_method == "single_point": return G1DListCrossoverSinglePoint

            # Dual-point crossover
            elif self.config.crossover_method == "two_point": return G1DListCrossoverTwoPoint

            # Uniform
            elif self.config.crossover_method == "uniform": return G1DListCrossoverUniform

            # OX
            elif self.config.crossover_method == "OX": return G1DListCrossoverOX

            # Edge
            elif self.config.crossover_method == "edge": return G1DListCrossoverEdge

            # Cut crossfill
            elif self.config.crossover_method == "cut_crossfill": return G1DListCrossoverCutCrossfill

            # Real SBX
            elif self.config.crossover_method == "real_SBX": return G1DListCrossoverRealSBX

            # Invalid
            else: raise ValueError("Invalid crossover method for one-dimensional genomes")

        # 2D genome
        elif self.config.genome_dimension == 2:

            # Uniform
            if self.config.crossover_method == "uniform": return G2DListCrossoverUniform

            # Vertical single-point
            elif self.config.crossover_method == "single_vertical_point": return G2DListCrossoverSingleVPoint

            # Horizontal single-point
            elif self.config.crossover_method == "single_horizontal_point": return G2DListCrossoverSingleHPoint

            # Invalid
            else: raise ValueError("Invalid crossover method for two-dimensional genomes")

        # Not supported number of dimensions
        else: raise ValueError("Dimensions > 2 are not supported")

    # -----------------------------------------------------------------

    def initialize_engine(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Initializing the genetic engine ...")

        # 1. Create the genetic engine
        self.engine = GeneticEngine(self.initial_genome)

        # 2. Set options
        self.set_engine_options()

        # 3. Set plotter
        self.set_engine_plotter(kwargs)

        # 4. Set database and statistics
        self.set_engine_database_and_statistics()

        # 5. Set kwargs
        self.set_engine_kwargs(kwargs)

    # -----------------------------------------------------------------

    def set_engine_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting genetic engine options ...")

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

    def set_engine_plotter(self, kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Setting plotter for the genetic engine ...")

        # Set generations plotter
        self.generations_plotter = kwargs.pop("generations_plotter", None)
        if self.generations_plotter is not None:
            self.generations_plotter.output_path = self.plot_path
            self.engine.stepCallback.set(self.generations_plotter.add_generation)

    # -----------------------------------------------------------------

    def set_engine_database_and_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the database and statistics adapters ...")

        # Set the database adapter
        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Set the adapter for the statistics table
        self.statistics.open(self.engine)
        self.engine.setDBAdapter(self.statistics)

    # -----------------------------------------------------------------

    def set_engine_kwargs(self, kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Setting the kwargs for the engine functions ...")

        # Get kwargs for the different functions
        evaluator_kwargs = kwargs.pop("evaluator_kwargs") if "evaluator_kwargs" in kwargs else None
        initializator_kwargs = kwargs.pop("initializator_kwargs") if "initializator_kwargs" in kwargs else None
        mutator_kwargs = kwargs.pop("mutator_kwargs") if "mutator_kwargs" in kwargs else None
        crossover_kwargs = kwargs.pop("crossover_kwargs") if "crossover_kwargs" in kwargs else None

        # Set kwargs
        if evaluator_kwargs is not None: self.engine.set_kwargs("evaluator", evaluator_kwargs)
        if initializator_kwargs is not None: self.engine.set_kwargs("initializator", initializator_kwargs)
        if mutator_kwargs is not None: self.engine.set_kwargs("mutator", mutator_kwargs)
        if crossover_kwargs is not None: self.engine.set_kwargs("crossover", crossover_kwargs)

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
