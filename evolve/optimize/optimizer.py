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
import numpy as np
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..core.engine import GeneticEngine, RawScoreCriteria
from ...core.basics.range import RealRange, IntegerRange
from ...core.tools import formatting as fmt
from ...core.tools import stringify
from ..core.adapters import DBFileCSV, DBSQLite, PopulationsFile
from ...core.tools import filesystem as fs
from ...core.tools import types
from ...core.tools import sequences
from ...core.tools.serialization import write_dict
from ..core.population import Population, NamedPopulation
from ...core.tools import numbers
from ...core.tools.stringify import tostr
from .components import get_genome_type, is_1d_genome, is_2d_genome, get_mutator, create_genome, get_crossover_method, get_crossover, get_selector, get_scaling, get_initializator
from ...core.tools.random import setup_prng
from .parameters import get_binary_genome_from_scaled_parameters, scale_parameter_sets
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class Optimizer(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Optimizer, self).__init__(*args, **kwargs)

        # The database adapter
        self.database = None

        # The statistics table adapter
        self.statistics = None

        # The populations file
        self.populations = None

        # The intial genome
        self.initial_genome = None

        # The initial population
        self.initial_population = None

        # The genetic algorithm engine
        self.engine = None

        # The best individual
        self.best = None

        # The generations plotter
        self.generations_plotter = None

        # The plot path
        self.plot_path = None

        # The parameter minima and maxima (for heterogeneous 1D genomes)
        self.parameter_minima = None
        self.parameter_maxima = None

        # The parameter centers and sigmas (for heterogeneous 1D genomes with Gaussian mutator and/or initializer)
        self.parameter_centers = None
        self.parameter_sigmas = None

        # The parameter range (for homogeneous 1D list genomes)
        self.parameter_range = None

        # Parameters for the initial generation
        self.initial_parameters = None

        # Number of digits for the parameters
        self.ndigits = None

        # Number of binary digits for the parameters
        self.nbits = None

        # The scales for the different parameters
        self.scales = None

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_base_type(self):

        """
        This function ...
        :return:
        """

        # Parameter range is defined, homogeneous
        if self.parameter_range is not None:

            if isinstance(self.parameter_range, IntegerRange): return "integer"
            elif isinstance(self.parameter_range, RealRange): return "real"
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

        # Parameter centers and sigmas are defined:
        # NO: DOESN'T SAY ANYTHING: CENTER AND SIGMA CAN BE REAL BUT WE STILL WANT TO GENERATE INTEGER GENES !!
        # (e.g. HeterogeneousListMutatorIntegerGaussian)
        #elif self.parameter_centers is not None and self.parameter_sigmas is not None: pass

        # Parameter type is defined in the configuration
        elif self.config.parameter_type is not None:

            # Check that the parameter type is set
            if self.config.parameter_type is None: raise ValueError("Parameter type (real, integer) must be set when parameter_centers and parameter_sigmas are defined")
            return self.config.parameter_type

        # No clue
        else: raise ValueError("Parameter type cannot be determined or is invalid")

    # -----------------------------------------------------------------

    @lazyproperty
    def is_integer_parameter(self):

        """
        This function ...
        :return: 
        """

        return self.parameter_base_type == "integer"

    # -----------------------------------------------------------------

    @lazyproperty
    def is_real_parameter(self):

        """
        THis function ...
        :return: 
        """

        return self.parameter_base_type == "real"

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

        # Set the parameter range (for uniform genome lists)
        if "parameter_range" in kwargs: self.parameter_range = kwargs.pop("parameter_range")
        else: self.parameter_range = self.config.parameter_range

        # Prepare the parameter ranges
        self.prepare_parameter_ranges()

        # Set initial parameters
        if "initial_parameters" in kwargs: self.initial_parameters = kwargs.pop("initial_parameters")

        # Check whether the number of initial parameter sets is equal to the number of individuals per generation
        if self.initial_parameters is not None:
            if len(self.initial_parameters) != self.config.nindividuals: raise ValueError("The number of initial parameter sets should be equal to the configured number of individuals per generation")

        # Set ndigits
        if "ndigits" in kwargs: self.ndigits = kwargs.pop("ndigits")

        # Set nbits
        if "nbits" in kwargs: self.nbits = kwargs.pop("nbits")

        # Check whether NBITS is given if binary_string is genome representation
        if self.binary_string_genome and self.nbits is None: raise ValueError("Number of bits must be specified for binary genomes")

        # Debugging
        if self.nbits is not None: log.debug("Number of bits for each parameter: " + " ".join(str(n) for n in self.nbits))

        # Get the scales
        if "scales" in kwargs: self.scales = kwargs.pop("scales")

    # -----------------------------------------------------------------

    def prepare_parameter_ranges(self):

        """
        This function ...
        :return: 
        """

        # If minima and maxima are defined but not parameter_range, but the genome is not necessarily a heterogeneous genome,
        # we need to set parameter_range as equal to both parameter ranges IF POSSIBLE ?
        if self.parameter_range is None and self.parameter_minima is not None and self.parameter_maxima is not None:

            # Check if they're the same
            if sequences.all_equal(self.parameter_minima) and sequences.all_equal(self.parameter_maxima):

                # Get the minimal and maximal value
                min_value = self.parameter_minima[0]
                max_value = self.parameter_maxima[0]

                # Create the parameter range
                if types.is_integer_type(min_value):
                    if not types.is_integer_type(max_value): raise ValueError("Minimum and maximum values should be of the same type")
                    self.parameter_range = IntegerRange(min_value, max_value)
                elif types.is_real_type(min_value):
                    if not types.is_real_type(max_value): raise ValueError("Minimum and maximum values should be of the same type")
                    self.parameter_range = RealRange(min_value, max_value)
                else: raise ValueError("Unknown type for maxima and minima")

                # Not the same
            elif not self.config.heterogeneous: raise ValueError("Parameter range must be defined for non-heterogeneous genomes")

        # Debugging
        if self.parameter_range is not None: log.debug("The parameter range is " + tostr(self.parameter_range, fancy=True))

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_scales(self):

        """
        This function ...
        :return: 
        """

        if self.scales is None: return [self.config.default_scale] * self.config.nparameters # DOESN'T WORK FOR 2D GENOMES OF COURSE
        else: return self.scales

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_minima_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.parameter_minima is None: return None

        minima = []

        # Loop over the scales and minima
        for scale, minimum in zip(self.parameter_scales, self.parameter_minima):

            # Convert to sclae
            value = numbers.to_scale(minimum, scale)

            # Add the minimum value (in lin or log space)
            minima.append(value)

        # Return the minima
        return minima

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_maxima_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.parameter_maxima is None: return None

        maxima = []

        # Loop over the scales and minima
        for scale, maximum in zip(self.parameter_scales, self.parameter_maxima):

            # Convert to scale
            value = numbers.to_scale(maximum, scale)

            # Add the maximum value (in lin or log space)
            maxima.append(value)

        # Return the maxima
        return maxima

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_centers_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.parameter_centers is None: return None

        centers = []

        # Loop over the scales and centers
        for scale, center in zip(self.parameter_scales, self.parameter_centers):

            # Convert to scale
            value = numbers.to_scale(center, scale)

            # Add the center value
            centers.append(value)

        # Return the centers
        return centers

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_sigmas_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.parameter_sigmas is None: return None

        sigmas = []

        # Loop over the scales and sigmas
        for scale, sigma in zip(self.parameter_scales, self.parameter_sigmas):

            # Convert to scale
            value = numbers.to_scale(sigma, scale)

            # Add the sigma value
            sigmas.append(value)

        # Return the sigmas
        return sigmas

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_range_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.parameter_range is None: return None

        if not sequences.all_equal(self.parameter_scales): return None # raise ValueError("Cannot return a scaled parameter ranges because the scales are different for the different parameters")

        else:

            scale = self.parameter_scales[0]

            # Convert and return
            if scale == "linear": return self.parameter_range
            elif scale == "logarithmic": return RealRange(np.log10(self.parameter_range.min), np.log10(self.parameter_range.max))
            else: raise ValueError("Invalid scale " + scale)

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_parameters_scaled(self):

        """
        This function ...
        :return: 
        """

        if self.initial_parameters is None: return None

        # Scale the initial parameter sets according to the parameter scales
        return scale_parameter_sets(self.initial_parameters, self.parameter_scales)

    # -----------------------------------------------------------------

    def initialize_prng(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the random number generator ...")

        # Set the seed of the PRNG
        setup_prng(self.config.seed)

    # -----------------------------------------------------------------

    def initialize_adapters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Initializing the adapters ...")

        # 1. Initialize the statistics table
        self.initialize_statistics()

        # 2. Initialize databse
        self.initialize_database()

        # 3. Initilaize the populations table
        self.initialize_populations()

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
        self.statistics = DBFileCSV(filename=filepath, frequency=self.config.statistics_frequency, reset=reset,
                                    identify=self.config.run_id, name=self.config.statistics_name)

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
        self.database = DBSQLite(dbname=filepath, identify=self.config.run_id, resetDB=reset,
                                 commit_freq=self.config.database_commit_frequency, frequency=self.config.database_frequency,
                                 resetIdentify=False, name=self.config.database_name)

    # -----------------------------------------------------------------

    def initialize_populations(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Initializing the populations ...")

        # Determine the file path
        if self.config.writing.populations_path is not None: filepath = fs.absolute_or_in(self.config.writing.populations_path, self.output_path)
        else: filepath = self.output_path_file("populations.dat")

        # Check the file
        if fs.is_file(filepath):
            log.debug("The populations file is already present. Adding a new run to the file.")
            reset = False
        else:
            log.debug("Creating a new populations file ...")
            reset = True

        # Create the populations file adapter
        self.populations = PopulationsFile(filepath=filepath, identify=self.config.run_id, reset=reset,
                                           frequency=self.config.populations_frequency, name=self.config.populations_name)

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

    @property
    def genome_type(self):

        """
        This function ...
        :return: 
        """

        # Return configured option
        if self.initial_genome is None: return self.config.genome_type

        # Return the genome type
        return get_genome_type(self.initial_genome)

    # -----------------------------------------------------------------

    @property
    def list_genome(self):

        """
        This function ...
        :return: 
        """

        return self.genome_type == "list"

    # -----------------------------------------------------------------

    @property
    def binary_string_genome(self):

        """
        This function ...
        :return: 
        """

        return self.genome_type == "binary_string"

    # -----------------------------------------------------------------

    @property
    def internal_population(self):

        """
        This function ...
        :return:
        """

        if self.engine is None: return None
        else: return self.engine.population

    # -----------------------------------------------------------------

    @property
    def crossover_method(self):

        """
        This function ...
        :return:
        """

        if self.initial_genome is not None: return get_crossover_method(self.initial_genome.crossover)
        else: return self.config.crossover_method

    # -----------------------------------------------------------------

    def create_genome(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the genome ...")

        # Number of digits has to be specified
        if self.binary_string_genome and self.ndigits is None: raise ValueError("Number of digits (in base-10) has to be specified for binary string conversion")

        # Create the genome
        self.initial_genome = create_genome(self.genome_dimension, self.genome_type, self.config.nparameters, self.config.nparameters2, nbits=self.total_nbits)

    # -----------------------------------------------------------------

    @property
    def total_nbits(self):

        """
        This function ...
        :return: 
        """

        total_nbits = sum(self.nbits)
        return total_nbits

    # -----------------------------------------------------------------

    @lazyproperty
    def bit_slices(self):

        """
        This function ...
        :return: 
        """

        slices = []

        tempsum = 0
        for index in range(len(self.nbits)):
            nbits = self.nbits[index]
            slices.append(slice(tempsum, tempsum+nbits))
            tempsum += nbits

        return slices

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
        if self.parameter_range_scaled is not None: self.initial_genome.setParams(rangemin=self.parameter_range_scaled.min, rangemax=self.parameter_range_scaled.max)
        if self.config.best_raw_score is not None: self.initial_genome.setParams(bestrawscore=self.config.best_raw_score)
        if self.config.round_decimal is not None: self.initial_genome.setParams(rounddecimal=self.config.round_decimal)

        # Debugging
        log.debug("Setting parameter minima and maxima ...")

        # Set minima or maxima for heterogeneous genome lists
        if self.parameter_minima_scaled is not None: self.initial_genome.setParams(minima=self.parameter_minima_scaled)
        if self.parameter_maxima_scaled is not None: self.initial_genome.setParams(maxima=self.parameter_maxima_scaled)

        # Debugging
        log.debug("Setting parameter centers and sigmas ...")

        # Set parameter centers and sigmas
        if self.parameter_centers_scaled is not None: self.initial_genome.setParams(centers=self.parameter_centers_scaled)
        if self.parameter_sigmas_scaled is not None: self.initial_genome.setParams(sigmas=self.parameter_sigmas_scaled)

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
    def is_1d_genome(self):

        """
        This function ...
        :return:
        """

        return is_1d_genome(self.initial_genome)

    # -----------------------------------------------------------------

    @property
    def is_2d_genome(self):

        """
        This function ...
        :return:
        """

        return is_2d_genome(self.initial_genome)

    # -----------------------------------------------------------------

    @property
    def genome_size(self):

        """
        This function ...
        :return:
        """

        # Value for 1D, tuple for 2D

        # Get directly from initial genome
        if self.initial_genome is not None: return self.initial_genome.getSize()

        # 1D
        elif self.config.nparameters2 is None:

            # 1D
            if self.list_genome: return self.config.nparameters
            elif self.binary_string_genome: return self.total_nbits
            else: raise ValueError("Invalid genome type")

        # 2D
        else:

            # 2D
            if self.list_genome: return self.config.nparameters, self.config.nparameters2
            elif self.binary_string_genome: raise NotImplementedError("Not implemented")
            else: raise ValueError("Invalid genome type")

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return: 
        """

        # Is this the best implementation?
        if self.list_genome: return len(self.initial_genome)
        elif self.binary_string_genome: return len(self.nbits)
        else: raise ValueError("Genome type not recognized")

    # -----------------------------------------------------------------

    def get_initializator(self):

        """
        This function ...
        :return:
        """

        # Get
        return get_initializator(self.genome_type, parameter_type=self.parameter_base_type, hetero=self.config.heterogeneous)

    # -----------------------------------------------------------------

    def get_mutator(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Getting the mutator type ...")

        # Get the mutator class
        if self.list_genome: return get_mutator(self.genome_type, self.genome_dimension, self.config.mutation_method, parameter_type=self.parameter_base_type, hetero=self.config.heterogeneous)
        elif self.binary_string_genome: return get_mutator(self.genome_type, self.genome_dimension, self.config.binary_mutation_method)

    # -----------------------------------------------------------------

    def get_crossover(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the crossover type ...")

        # Get the crossover
        return get_crossover(self.genome_type, self.genome_dimension, self.config.crossover_method, genome_size=self.genome_size)

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
        self.create_engine()

        # 2. Set options
        self.set_engine_options()

        # 3. Set the selector
        self.set_engine_selector()

        # 4. Set plotter
        self.set_engine_plotter(kwargs)

        # 5. Set database and statistics
        self.set_engine_adapters()

        # 6. Set kwargs
        self.set_engine_kwargs(kwargs)

    # -----------------------------------------------------------------

    @property
    def has_initial_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.initial_parameters is not None

    # -----------------------------------------------------------------

    def create_engine(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the genetic engine ...")

        #  Create the engine, passing the initial genome, and the 'named_individuals' flag
        if self.has_initial_parameters: self.engine = GeneticEngine(self.create_initial_population())
        else: self.engine = GeneticEngine(self.initial_genome, named_individuals=self.config.named_individuals)

    # -----------------------------------------------------------------

    def create_initial_population(self):

        """
        This function doesn't set any attributes but is used to create an initial population
        :return: 
        """

        # Inform the user
        log.info("Creating the initial population ...")

        # Initialize the population
        if self.config.named_individuals: population = NamedPopulation(size=self.config.nindividuals)
        else: population = Population(size=self.config.nindividuals)

        # Loop over the parameter sets
        for parameters in self.initial_parameters_scaled:

            # Make a new genome by cloning the initial genome
            genome = self.initial_genome.clone()

            # Clear the genes
            genome.clear()

            # Set the genes
            if self.list_genome: genome.set_genes(parameters)
            elif self.binary_string_genome:

                # PARAMETERS IS ALREADY SCALED!!!
                genes = get_binary_genome_from_scaled_parameters(parameters, self.parameter_minima_scaled, self.parameter_maxima_scaled, self.nbits, gray=self.config.gray_code)

                # Set the binary string genome
                genome.set_genes(genes)

            # Invalid
            else: raise ValueError("Unrecognized genome type")

            # Add the genome (individual) to the population
            population.append(genome)

        #print(population.internalPop.names)

        # Return the population
        return population

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
        self.engine.setMinimax(self.config.min_or_max)
        self.engine.setGenerations(self.config.ngenerations)
        if self.config.crossover_rate is not None: self.engine.setCrossoverRate(self.config.crossover_rate)
        self.engine.setPopulationSize(self.config.nindividuals)
        self.engine.setMutationRate(self.config.mutation_rate)

        # Set elitism options
        self.engine.setElitism(self.config.elitism)
        self.engine.setElitismReplacement(self.config.nelite_individuals)

    # -----------------------------------------------------------------

    def set_engine_selector(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the genetic engine selector ...")

        # Get the selector
        selector = get_selector(self.config.selector_method)

        # Set the selector
        self.engine.selector.set(selector)

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

    def set_engine_adapters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the engine adapters ...")

        # Set the database adapter
        self.database.open(self.engine)
        self.engine.add_database_adapter(self.database)

        # Set the adapter for the statistics table
        self.statistics.open(self.engine)
        self.engine.add_database_adapter(self.statistics)

        # Set the populations adapter
        self.populations.open(self.engine)
        self.engine.add_database_adapter(self.populations)

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

    def initialize_population(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Initializing the population ...")

        # Get the internal population
        self.initial_population = self.engine.get_population()

        # Set population minimax type
        self.set_population_minimax()

        # Set scaling method
        self.set_population_scaling()

    # -----------------------------------------------------------------

    def set_population_minimax(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the population minimize/maximize flag ...")

        # Set
        self.initial_population.setMinimax(self.config.min_or_max)

    # -----------------------------------------------------------------

    def set_population_scaling(self):

        """
        This function ...
        :return: 
        """
        
        # Inform the user
        log.info("Setting the population scaling method ...")

        # Get the scaling function
        scaling = get_scaling(self.config.scaling_method)

        # Set scaling
        self.initial_population.scaleMethod.set(scaling)

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

    def write_database(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the database ...")

        # Commit all changes and close the database
        self.database.commit_and_close()

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the statistics ...")

        # Commit all changes and close
        self.statistics.commit_and_close()

    # -----------------------------------------------------------------

    def write_populations(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the populations ...")

        # Commit all changes and close
        self.populations.commit_and_close()

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the parameters of the genetic engine and population ...")

        # Write engine parameters
        self.write_engine_parameters()

        # Write genome parameters
        if self.initial_genome is not None: self.write_genome_parameters()

        # Write population parameters
        if self.initial_population is not None: self.write_population_parameters()

    # -----------------------------------------------------------------

    def write_engine_parameters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the parameters of the engine ...")

        # Get the parameters
        parameters = self.engine.get_parameters()

        # Function slots
        selectors = self.engine.selector.names
        callback = self.engine.stepCallback.names
        termination = self.engine.terminationCriteria.names

        # Add to parameters
        parameters["selectors"] = selectors
        parameters["callback"] = callback
        parameters["termination"] = termination

        # Other properties
        minimax = self.engine.minimax
        generations = self.engine.nGenerations
        crossover_rate = self.engine.pCrossover
        #population_size = self.engine.popSize
        mutation_rate = self.engine.pMutation
        elitism = self.engine.elitism
        nelite_individuals = self.engine.nElitismReplacement

        # Add to parameters
        parameters["minimax"] = minimax
        parameters["generations"] = generations
        parameters["crossover_rate"] = crossover_rate
        #parameters["population_size"] = population_size
        parameters["mutation_rate"] = mutation_rate
        parameters["elitism"] = elitism
        parameters["nelite_individuals"] = nelite_individuals

        # Determine the path
        path = self.output_path_file("engine.param")

        # Debugging
        log.debug("Writing parameter dictionary to '" + path + "' ...")

        # Write as dict
        write_dict(parameters, path)

    # -----------------------------------------------------------------

    def write_genome_parameters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the parameters of the initial genome ...")

        # Get the parameters
        parameters = self.initial_genome.get_parameters()

        # Function slots
        initializators = self.initial_genome.initializator.names
        mutators = self.initial_genome.mutator.names
        crossover = self.initial_genome.crossover.names

        # Add to parameters
        parameters["initializators"] = initializators
        parameters["mutators"] = mutators
        parameters["crossover"] = crossover

        # Other properties
        dimension = self.initial_genome.dimension
        size = self.initial_genome.getSize()

        # Add to parameters
        parameters["dimension"] = dimension
        parameters["size"] = size

        # Determine the path
        path = self.output_path_file("genome.param")

        # Debugging
        log.debug("Writing parameter dictionary to '" + path + "' ...")

        # Write as dict
        write_dict(parameters, path)

    # -----------------------------------------------------------------

    def write_population_parameters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the parameters of the initial population ...")

        # Get the parameters
        parameters = self.initial_population.get_parameters()

        # Function slots
        scaling_functions = self.initial_population.scaleMethod.names
        parameters["scaling"] = scaling_functions

        # Determine the path
        path = self.output_path_file("population.param")

        # Debugging
        log.debug("Writing parameter dictionary to '" + path + "' ...")

        # Write as dict
        write_dict(parameters, path)

    # -----------------------------------------------------------------

    def write_best(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the best individual ...")

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
    print(stringify.stringify_list_fancy([gene for gene in genome], width=100, delimiter=", ", lines_prefix="    ")[1])
    print("")

    print(fmt.yellow + "  Score: " + fmt.reset + str(score))

    print("")

    print(fmt.yellow + "  Fitness: " + fmt.reset + str(fitness))

    print("")

# -----------------------------------------------------------------

def get_best_parameter_values_per_generation(database_path, populations_path, run_id):

    """
    This function ...
    :param database_path:
    :param populations_path:
    :param run_id:
    :return:
    """

    import pandas as pd
    from ..analyse.database import get_best_individual_key_for_generation
    from .stepwise import load_populations

    key = []
    for g_n in np.arange(11):
        key.append(get_best_individual_key_for_generation(database_path, run_id, g_n, "min"))

    pop = load_populations(populations_path)

    vred = pop.values()
    va_t = []
    ke_t = []

    va = pd.DataFrame(columns=[['Name', 'Binary']])
    for ind, k in enumerate(key):
        for i in vred:
            ke_t.append(i[ind].keys()[0])
            va_t.append(i[ind].values()[0])

    va['Name'] = ke_t
    va['Binary'] = va_t

    #va['Values'] = None

    va.index.name = 'Generation'

    #
    return va

# -----------------------------------------------------------------
