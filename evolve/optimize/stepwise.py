#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.stepwise Contains the StepWiseOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..core.engine import GeneticEngine
from ...core.tools import filesystem as fs
from ...core.tools.random import save_state, load_state
from ...core.basics.configuration import Configuration
from ..core.adapters import DBFileCSV, DBSQLite, PopulationsFile
from .optimizer import Optimizer
from ..core.population import NamedPopulation
from .tables import ElitismTable, CrossoverTable, ScoresTable, RecurrenceTable
from ..analyse.database import load_database, get_score_for_individual
from .parameters import get_parameters_from_genome, get_binary_genome_from_parameters
from .parameters import equal_genomes
from ..core import constants
from ...core.basics.map import Map
from ...core.basics.configurable import write_input
from ...core.tools import types
from ...core.tools.stringify import tostr
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

newborns_filename = "newborns.dat"
parents_filename = "parents.dat"

# -----------------------------------------------------------------

class StepWiseOptimizer(Optimizer):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(StepWiseOptimizer, self).__init__(*args, **kwargs)

        # The scores of the past generation
        self.scores = None

        # The individual names in the same order as the scores
        self.scores_names = None

        # A check for the parameters of the models that are scored
        self.scores_check = None

        # The elitism table
        self.elitism_table = None

        # The crossover table
        self.crossover_table = None

        # The scores table
        self.scores_table = None

        # The previous population
        self.previous_population = None

        # The previous recurrence table
        self.previous_recurrence = None

        # Recurrence table
        self.recurrence_table = None

        # The input that was passed to the run function
        self.input = None

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, dir_path, run_id, statistics_name=None, database_name=None, populations_name=None,
                       frequency=None, commit_frequency=None):

        """
        This function ...
        :param dir_path:
        :param run_id:
        :param statistics_name:
        :param database_name:
        :param populations_name:
        :param frequency:
        :param commit_frequency:
        :return:
        """

        # Debugging
        log.debug("Loading optimizer for run '" + run_id + "' from '" + dir_path + "' ...")

        # Determine the path to the engine
        engine_path = fs.join(dir_path, "engine.pickle")

        # Determine the path to the random generator
        prng_path = fs.join(dir_path, "prng.pickle")

        # Determine the path to the configuration
        config_path = fs.join(dir_path, "optimizer.cfg")

        # Set the path to the statistics file (opening is done during initialization or evolution)
        statistics_path = fs.join(dir_path, "statistics.csv")

        # Set the path to the database (opening is done during initialization or evolution)
        database_path = fs.join(dir_path, "database.db")

        # Set the path to the populations file (opening is done during initialization or evolution)
        populations_path = fs.join(dir_path, "populations.dat")

        # Create the optimizer instance and return it
        return cls.from_paths(dir_path, engine_path, prng_path, config_path, statistics_path, database_path,
                              populations_path, run_id, statistics_name=statistics_name, database_name=database_name,
                              populations_name=populations_name, frequency=frequency, commit_frequency=commit_frequency)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, output_path, engine_path, prng_path, config_path, statistics_path, database_path,
                   populations_path, run_id, statistics_name=None, database_name=None, populations_name=None,
                   frequency=None, commit_frequency=None):

        """
        This function ...
        :param output_path:
        :param engine_path:
        :param prng_path:
        :param config_path:
        :param statistics_path:
        :param database_path:
        :param populations_path:
        :param run_id:
        :param statistics_name:
        :param database_name:
        :param populations_name:
        :param frequency:
        :param commit_frequency:
        :return:
        """

        # Debugging
        log.debug("")
        log.debug("Loading optimizer for run '" + run_id + "' from:")
        log.debug("")
        log.debug(" - Engine: " + engine_path)
        log.debug(" - PRNG: " + prng_path)
        log.debug(" - Configuration: " + config_path)
        log.debug(" - Statistics: " + statistics_path)
        log.debug(" - Database: " + database_path)
        log.debug(" - Populations: " + populations_path)
        log.debug(" - Output: " + output_path)
        log.debug("")

        # Set optional values to default values
        if statistics_name is None: statistics_name = constants.CDefCSVName
        if database_name is None: database_name = constants.CDefSQLiteName
        if populations_name is None: populations_name = constants.CDefPopulationsName
        if frequency is None: frequency = 1
        if commit_frequency is None: commit_frequency = constants.CDefSQLiteStatsCommitFreq

        # Load the engine
        engine = GeneticEngine.from_file(engine_path)

        # Load the random state
        load_state(prng_path)

        # Load the configuration
        config = Configuration.from_file(config_path)

        # Create the optimizer instance
        optimizer = cls(config)

        # Set the engine
        optimizer.engine = engine

        # Load the statistics (opening is done during initialization or evolution)
        if not fs.is_file(statistics_path): raise IOError("The statistics file could not be found at '" + statistics_path + "'")
        # Check whether there is data?
        #if not fs.contains_lines(statistics_path): raise IOError("The statistics file is empty")
        log.debug("Loading the statistics from '" + statistics_path + "' ...")
        optimizer.statistics = DBFileCSV(filename=statistics_path, reset=False, identify=run_id, name=statistics_name, frequency=frequency)

        # Load the database (opening is done during initialization or evolution)
        if not fs.is_file(database_path): raise IOError("The database could not be found at '" + database_path + "'")
        log.debug("Loading the database from '" + database_path + "' ...")
        optimizer.database = DBSQLite(dbname=database_path, resetDB=False, identify=run_id, resetIdentify=False, name=database_name, frequency=frequency, commit_freq=commit_frequency)

        # Load the populations file (opening is done during initialization or evoluation)
        if not fs.is_file(populations_path): raise IOError("The populations file could not be found at '" + populations_path + "'")
        # Check whether there is data? -> NO, BECAUSE DATA IS ONLY ADDED AT THE END OF THE SECOND RUN (WHEN INITIAL HAS BEEN SCORED, AND GENERATION0 GENERATED)
        #if not fs.contains_lines(populations_path): raise IOError("The populations file is empty")
        log.debug("Loading the populations file from '" + populations_path + "' ...")
        optimizer.populations = PopulationsFile(filepath=populations_path, reset=False, identify=run_id, name=populations_name, frequency=frequency)

        # Set the path
        optimizer.config.output = output_path

        # Set the individual paths
        optimizer.config.writing.engine_path = engine_path
        optimizer.config.writing.prng_path = prng_path
        optimizer.config.writing.config_path = config_path
        optimizer.config.writing.statistics_path = statistics_path
        optimizer.config.writing.database_path = database_path
        optimizer.config.writing.populations_path = populations_path

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Initialize
        if not self.initialized: self.initialize(**kwargs)

        # 3. Finish
        elif self.config.finish: self.finish()

        # 3. Evolve
        else: self.evolve()

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

        # Get the input
        self.input = copy.copy(kwargs)

        # Call the setup function of the base class
        super(StepWiseOptimizer, self).setup(**kwargs)

        # Get the engine
        if "engine" in kwargs: self.engine = kwargs.pop("engine")

        # Get the scores
        if "scores" in kwargs: self.scores = kwargs.pop("scores")

        # Get the scores names
        if "scores_names" in kwargs: self.scores_names = kwargs.pop("scores_names")

        # Get the scores check
        if "scores_check" in kwargs: self.scores_check = kwargs.pop("scores_check")

        # Get the output path
        if "output" in kwargs: self.config.output = kwargs.pop("output")

        # Get the previous population
        if "previous_population" in kwargs: self.previous_population = kwargs.pop("previous_population")

        # Get the previous recurrence data
        if "previous_recurrence" in kwargs: self.previous_recurrence = kwargs.pop("previous_recurrence")

        # Set best to None: only when finish_evoluation is run, best should be set
        self.best = None

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing ...")

        # 1. Initialize the random number generator
        self.initialize_prng()

        # 2. Initialize the adapters
        self.initialize_adapters()

        # 3. Initialize genome
        self.initialize_genome(**kwargs)

        # 4. Initialize engine
        self.initialize_engine(**kwargs)

        # 5. Initialize the initial population
        self.initialize_population()

        # 6. Initialize the evolution (the initial generation)
        self.initialize_evolution()

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

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finishing ...")

        # Set the database adapters again
        self.set_engine_adapters()

        # Set the scores from the previous generation
        self.set_scores()

        # Get the best individual
        self.best = self.engine.finish_evolution()

    # -----------------------------------------------------------------

    @property
    def check_recurrence(self):

        """
        This function ...
        :return: 
        """

        return self.config.check_recurrence and not self.engine.is_initial_generation

        # if the engine is at its first generation, we just did initial, so no data available in the populations.dat
        # -> NOT TRUE, the GENERATION INDEX IS ONLY INCREMENTED AT THE SCORING PHASE

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolve ...")

        # Set the database and statistics adapters again
        self.set_engine_adapters()

        # Make sure we don't stop unexpectedly, always increment the number of generations we want
        # (it doesn't matter what the value is)
        self.engine.setGenerations(self.engine.getGenerations() + 1)

        # Set the scores from the previous generation
        self.set_scores()

        # Generate new population
        self.generate_new_population()

        # Check recurrence
        if self.check_recurrence: self.set_recurrence()

    # -----------------------------------------------------------------

    @property
    def all_scores(self):

        """
        This function ...
        :return: 
        """

        if self.previous_recurrence is None: return self.scores
        if self.previous_recurrence.nrecurrences == 0: return self.scores

        scores = []

        passed_scores = iter(self.scores)

        # Loop over the names in the previous population, because the order IS IMPORTANT
        for key in self.previous_population:

            if key in self.previous_recurrence: score = self.previous_recurrence.get_score_for_individual(key)
            else: score = passed_scores.next()

            # Add the score
            scores.append(score)

        # Check the length of the scores list
        if len(scores) != len(self.previous_population): raise RuntimeError("Something went wrong")

        # Return the scores
        return scores

    # -----------------------------------------------------------------

    @property
    def all_scores_names(self):

        """
        This function ...
        :return:
        """

        # If no names are passed, return None
        if self.scores_names is None: return None

        # No recurrency
        if self.previous_recurrence is None: return self.scores_names
        if self.previous_recurrence.nrecurrences == 0: return self.scores_names

        # Initialize a list for the names
        names = []

        passed_names = iter(self.scores_names)

        for key in self.previous_population:

            if key in self.previous_recurrence: name = key
            else: name = passed_names.next()

            # Add the name
            names.append(name)

        # Check the length of the list
        if len(names) != len(self.previous_population): raise RuntimeError("Something went wrong")

        # Return the names
        return names

    # -----------------------------------------------------------------

    @property
    def all_scores_check(self):

        """
        This function ...
        :return: 
        """

        # If no check is passed, return None
        if self.scores_check is None: return None

        # No recurrency
        if self.previous_recurrence is None: return self.scores_check
        if self.previous_recurrence.nrecurrences == 0: return self.scores_check

        # Initialize a list for the checks
        checks = []

        passed_checks = iter(self.scores_check)

        # Loop over the names in the previous population, because the order IS IMPORTANT
        for key in self.previous_population:

            # Recurrent
            if key in self.previous_recurrence: parameters = get_parameters_from_genome(self.previous_population[key], self.parameter_minima_scaled, self.parameter_maxima_scaled, self.nbits, self.parameter_scales, gray=self.config.gray_code)

            # Not recurrent
            else: parameters = passed_checks.next()

            # Add the parameters
            checks.append(parameters)

        # Check the length of the checks list
        if len(checks) != len(self.previous_population): raise RuntimeError("Something went wrong")

        # Return the checks
        return checks

    # -----------------------------------------------------------------

    @property
    def all_scores_check_converted(self):

        """
        This function ...
        :return: 
        """

        # If list genome: do nothing
        if self.list_genome: return self.all_scores_check
        elif self.binary_string_genome:

            # Initialize list for the checks
            checks = []

            scores_check = self.all_scores_check
            if scores_check is None: return None

            # Check whether parameter minima and maxima are defined
            if self.parameter_minima is None or self.parameter_maxima is None: raise ValueError("Parameter minima and maxima should be defined")

            # Loop over the (unscaled) parameter sets
            for parameters in scores_check:

                # Generate and add the binary string genome from the (unscaled) parameters
                binary_string = get_binary_genome_from_parameters(parameters, self.parameter_minima_scaled, self.parameter_maxima_scaled, self.nbits, self.parameter_scales, gray=self.config.gray_code)
                checks.append(binary_string)

            # Return the converted checks
            return checks

        # Invalid
        else: raise ValueError("Unrecognized genome type")

    # -----------------------------------------------------------------

    @lazyproperty
    def binary_parameters(self):

        """
        This function ...
        :return: 
        """

        # Create binary parameters map
        binary_parameters = Map()
        #binary_parameters.minima = self.parameter_minima_scaled
        #binary_parameters.maxima = self.parameter_maxima_scaled
        binary_parameters.minima = self.parameter_minima
        binary_parameters.maxima = self.parameter_maxima
        binary_parameters.ndigits = self.ndigits # list
        binary_parameters.nbits = self.nbits # list
        binary_parameters.gray = self.config.gray_code
        binary_parameters.scales = self.parameter_scales

        # Return
        return binary_parameters

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

        # Check whether the number of bits is defined
        if self.nbits is None: raise ValueError("The number of bits for each parameter should be defined")

        # Get scores and checks
        scores = self.all_scores
        names = self.all_scores_names
        checks = self.all_scores_check_converted

        # Debugging
        log.debug("All scores (with recurrent models) are: " + " ".join(str(score) for score in scores))
        log.debug("All names (with recurrent models) are: " + " ".join(name for name in names))
        log.debug("All checks (with recurrent models) are: " + " ".join(str(check) for check in checks))

        # Set the scores
        elitism_data = self.engine.set_scores(scores, names, checks, rtol=self.config.check_rtol, atol=self.config.check_atol, binary_parameters=self.binary_parameters)

        # Create elitism table
        if elitism_data is not None: self.elitism_table = ElitismTable.from_data(elitism_data)
        else: log.warning("No elitism has been performed in this generation")

        # Create the scores table
        self.scores_table = ScoresTable.from_data(self.population_scores, min_or_max=self.min_or_max)

    # -----------------------------------------------------------------

    def generate_new_population(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new population ...")

        # Generate the new population
        crossover_data = self.engine.generate_new_population()

        # Create the crossover table
        if crossover_data is not None: self.crossover_table = CrossoverTable.from_data(crossover_data, self.genome_type, self.crossover_method)
        else: raise RuntimeError("Could not get the crossover data")

    # -----------------------------------------------------------------

    @property
    def newborns(self):

        """
        This function ...
        :return:
        """

        return self.engine.new_population

    # -----------------------------------------------------------------

    @property
    def parents(self):

        """
        This function ...
        :return:
        """

        # NEW
        if self.config.finish: return None
        else: return self.internal_population

    # -----------------------------------------------------------------

    @property
    def population(self):

        """
        This function ...
        :return:
        """

        if self.newborns is not None: return self.newborns
        else: return self.parents

    # -----------------------------------------------------------------

    @property
    def population_scores(self):

        """
        This function return sthe scores of the internal population, as a dictionary
        :return:
        """

        return self.engine.scores

    # -----------------------------------------------------------------

    def set_recurrence(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Looking for recurrence of old individuals within the new population ...")

        # Initialize the recurrence table
        self.recurrence_table = RecurrenceTable()

        # Get run ID
        run_id = self.populations.identify

        # Load the populations data
        populations = load_populations(self.populations.filepath)
        populations_run = populations[run_id]

        # Index of the new generation
        generation = self.engine.currentGeneration + 1

        # Load the database
        database = load_database(self.database.dbName)

        # Loop over the individual names (of the newborns)
        for name in self.individual_names:

            # Get the individual
            individual = self.population[name]

            # Check recurrence
            generation_index, key, parameters, original_parameters = find_recurrent_individual(populations_run, individual, generation, rtol=self.config.recurrence_rtol, atol=self.config.recurrence_atol, binary_parameters=self.binary_parameters, return_comparison=True)

            # If not found, skip
            if generation_index is None: continue

            # Determine the generation name
            #generation_name = "Generation" + str(generation_index-1)

            # Debugging
            log.debug("Individual '" + name + "' is recurrent: individual '" + str(key) + "' from generation " + str(generation_index-1))

            # Otherwise, look for the (raw) score in the database
            score = get_score_for_individual(database, run_id, generation_index, key)

            # Debugging
            log.debug("The score of this individual was " + str(score))

            # Convert parameters to strings
            parameters_string = tostr(parameters)
            original_parameters_string = tostr(original_parameters)

            # Add entry to the recurrence table
            self.recurrence_table.add_entry(name, generation_index, key, score, parameters_string, original_parameters_string)

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

        # Inform the user
        log.info("Writing ...")

        # Write the input
        self.write_input()

        # Save the engine
        self.write_engine()

        # Save the state of the prng
        self.write_random_state()

        # Save the configuration
        self.write_config()

        # Save the database
        self.write_database()

        # Save the statistics
        self.write_statistics()

        # Write the populations
        self.write_populations()

        # Write the new generation
        if self.newborns is not None: self.write_newborns()

        # Write the internal generation
        if self.parents is not None: self.write_parents()

        # Write the best individual
        self.write_best()

        # Write the parameters
        self.write_parameters()

        # Write the crossover data
        if self.crossover_table is not None: self.write_crossover()

        # Write the scores
        if self.scores_table is not None: self.write_scores()

        # Write the elitism data
        if self.elitism_table is not None: self.write_elitism()

        # Write the recurrency data
        if self.recurrence_table is not None: self.write_recurrence()

    # -----------------------------------------------------------------

    def write_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the input ...")

        # Determine the path
        if self.config.writing.input_path is not None: input_path = fs.create_absolute_or_in(self.config.writing.input_path, self.output_path)
        else: input_path = self.output_path_directory("input")

        # Save the input
        write_input(self.input, input_path)

    # -----------------------------------------------------------------

    def write_engine(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the engine ...")

        # Determine the path to the engine
        if self.config.writing.engine_path is not None:

            # Single path
            if types.is_string_type(self.config.writing.engine_path):

                # Get the path
                engine_path = fs.absolute_or_in(self.config.writing.engine_path, self.output_path)

                # Save the engine
                self.engine.saveto(engine_path)

            # Multiple paths
            elif types.is_string_sequence(self.config.writing.engine_path):

                # Loop over the paths
                for path in self.config.writing.engine_path:

                    # Get path
                    engine_path = fs.absolute_or_in(path, self.output_path)

                    # Save the engine
                    self.engine.saveto(engine_path)

            # Invalid
            else: raise ValueError("Invalid type for engine_path: must be string or string sequence")

        # Path not specified
        else:

            # Get the path
            engine_path = self.output_path_file("engine.pickle")

            # Save the engine
            self.engine.saveto(engine_path)

    # -----------------------------------------------------------------

    def write_random_state(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the state of the random number generator ...")

        # Determine the path to the prng
        if self.config.writing.prng_path is not None:

            # Single path
            if types.is_string_type(self.config.writing.prng_path):

                # Determine the path
                prng_path = fs.absolute_or_in(self.config.writing.prng_path, self.output_path)

                # Save the state
                save_state(prng_path)

            # Multiple paths
            elif types.is_string_sequence(self.config.writing.prng_path):

                for path in self.config.writing.prng_path:

                    # Determine the path
                    prng_path = fs.absolute_or_in(path, self.output_path)

                    # Save the state
                    save_state(prng_path)

            # Invalid
            else: raise ValueError("Invalid type for prng_path: must be string or string sequence")

        # No path specified
        else:

            # Set the path
            prng_path = self.output_path_file("prng.pickle")

            # Save the prng state
            save_state(prng_path)

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the configuration file ...")

        # Determine the path
        if self.config.writing.config_path is not None:

            # Single path
            if types.is_string_type(self.config.writing.config_path):

                # Determine path
                config_path = fs.absolute_or_in(self.config.writing.config_path, self.output_path)

                # Save the configuration
                self.config.saveto(config_path)

            # Multiple paths
            elif types.is_string_sequence(self.config.writing.config_path):

                # Loop over the paths
                for path in self.config.writing.config_path:

                    # Get path
                    config_path = fs.absolute_or_in(path, self.output_path)

                    # Save the configuration
                    self.config.saveto(config_path)

            # Invalid
            else: raise ValueError("Invalid type for config_path: must be string or string sequence")

        # No path specified
        else:

            # Set the path
            config_path = self.output_path_file("optimizer.cfg")

            # Save the configuration
            self.config.saveto(config_path)

    # -----------------------------------------------------------------

    def get_individual(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return self.population[name]

    # -----------------------------------------------------------------

    def get_parameters(self, name):

        """
        This fucntion ...
        :param name: 
        :return: 
        """

        # Get the genome of the individual
        genome = self.population[name]

        # Get the real parameters, unscaled
        parameters = get_parameters_from_genome(genome, self.parameter_minima, self.parameter_maxima, self.nbits, self.parameter_scales, gray=self.config.gray_code)

        # Round ?? MAYBE NOT -> causes error, after this, the scores_check contains this rounded value,
        # after which it is scaled again to log scale, then represented in binary, then converted back to log and rounded again
        # to check --> TOO MANY ERRORS!
        #if self.ndigits is not None: parameters = round_parameters(parameters, self.ndigits)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    @property
    def is_named_population(self):

        """
        This function ...
        :return: 
        """

        return isinstance(self.internal_population, NamedPopulation)

    # -----------------------------------------------------------------

    @property
    def min_or_max(self):

        """
        This function ...
        :return:
        """

        return self.internal_population.minimax

    # -----------------------------------------------------------------

    @property
    def individual_names(self):

        """
        This function ...
        :return: 
        """

        if not self.is_named_population: raise ValueError("The population is not a named population")
        return self.population.names

    # -----------------------------------------------------------------

    @property
    def individual_keys(self):

        """
        This function ...
        :return: 
        """

        return self.population.keys

    # -----------------------------------------------------------------

    @property
    def new_individual_names(self):

        """
        This function ...
        :return: 
        """

        if not self.is_named_population: raise ValueError("The population is not a named population")

        if self.recurrence_table is None: return self.individual_names

        names = []
        for name in self.individual_names:

            if name in self.recurrence_table: continue
            names.append(name)

        return names

    # -----------------------------------------------------------------

    @property
    def new_individual_keys(self):

        """
        This function ...
        :return: 
        """

        if self.recurrence_table is None: return self.individual_keys

        keys = []
        for key in self.individual_keys:

            if key in self.recurrence_table: continue
            keys.append(key)

        return keys

    # -----------------------------------------------------------------

    def write_newborns(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the newborn genomes ...")

        # Determine the path
        if self.config.writing.newborns_path is not None: path = fs.absolute_or_in(self.config.writing.newborns_path, self.output_path)
        else: path = self.output_path_file(newborns_filename)

        # Write
        write_population(self.newborns, path)

    # -----------------------------------------------------------------

    def write_parents(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parent genomes ...")

        # Determine the path
        if self.config.writing.parents_path is not None: path = fs.absolute_or_in(self.config.writing.parents_path, self.output_path)
        else: path = self.output_path_file(parents_filename)

        # Write
        write_population(self.parents, path)

    # -----------------------------------------------------------------

    def write_crossover(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the crossover table ...")

        # Determine the path
        if self.config.writing.crossover_table_path is not None: path = fs.absolute_or_in(self.config.writing.crossover_table_path, self.output_path)
        else: path = self.output_path_file("crossover.dat")

        # Save the table
        self.crossover_table.saveto(path)

    # -----------------------------------------------------------------

    def write_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the scores table ...")

        # Determine the path
        if self.config.writing.scores_table_path is not None: path = fs.absolute_or_in(self.config.writing.scores_table_path, self.output_path)
        else: path = self.output_path_file("scores.dat")

        # Save the table
        self.scores_table.saveto(path)

    # -----------------------------------------------------------------

    def write_elitism(self):

        """
        This fucntion ...
        :return: 
        """

        # Inform the user
        log.info("Writing the elitism table ...")

        # Determine the path
        if self.config.writing.elitism_table_path is not None: path = fs.absolute_or_in(self.config.writing.elitism_table_path, self.output_path)
        else: path = self.output_path_file("elitism.dat")

        # Save the table
        self.elitism_table.saveto(path)

    # -----------------------------------------------------------------

    def write_recurrence(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the recurrence data ...")

        # Determine the path
        if self.config.writing.recurrent_path is not None: path = fs.absolute_or_in(self.config.writing.recurrent_path, self.output_path)
        else: path = self.output_path_file("recurrence.dat")

        # Save the dictionary
        self.recurrence_table.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------

def write_population(population, path):

    """
    This function ...
    :param population: 
    :param path: 
    :return: 
    """

    # Loop over the individuals
    entries = []
    for key in population.keys:

        # Get the individual
        individual = population[key]

        # Add entry to the list
        entries.append((key, individual.genomeList)) # .genomeList works for all G1D genomes (list, binary string, ...)

    # Loop over the entries, write to file
    with open(path, 'w') as population_file:

        # Loop over the entries
        for key, genes in entries: print(key + " " + str(genes), file=population_file)

# -----------------------------------------------------------------

def load_population(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    # ORDERED IS IMPORTANT!!! KEYS ARE USED TO GET LIST OF INDIVIDUAL NAMES FOR CREATING ALL_SCORES
    population = OrderedDict()

    # loop over the lines in the file
    for line in fs.read_lines(path):

        parts = line.split(" ")
        key = parts[0]
        rest = line.split(key + " ")[1]
        genome = eval(rest)

        population[key] = genome

    # Return the population
    return population

# -----------------------------------------------------------------

def load_populations(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    populations = OrderedDict()

    # Loop over the lines in the file
    for line in fs.read_lines(path):

        #print("LINE", line)

        parts = line.split(" ")

        run_name = parts[0]
        generation = int(parts[1])
        key = parts[2]

        # Get the genome
        rest = line.split(key + " ")[1]
        genome = eval(rest)

        # Create empty list for each run initially
        if run_name not in populations: populations[run_name] = []

        # Check whether there is place in the list for this population yet (e.g. generation = 2, len must be at least 3)
        if len(populations[run_name]) < generation + 1: populations[run_name].append(OrderedDict())

        # Add the genome
        populations[run_name][generation][key] = genome

    # Return the populations data
    return populations

# -----------------------------------------------------------------

def find_recurrent_individual(populations, individual, current_generation, rtol=1e-5, atol=1e-8, binary_parameters=None, return_comparison=False):

    """
    This function ...
    :param populations:
    :param individual:
    :param current_generation:
    :param rtol:
    :param atol:
    :param binary_parameters:
    :param return_comparison:
    :return: 
    """

    # Look for a match with previous generations (populations) of the same run
    # Loop over the previous generations
    generation_indices = range(current_generation)
    for generation_index in reversed(generation_indices): # look from the last to the first generations

        # Loop over the individuals in this generation
        for key in populations[generation_index]:

            # Get the genome
            genome = populations[generation_index][key]

            # If comparison values have to be returned
            if return_comparison:

                # If equal, return
                equal, parameters, original_parameters = equal_genomes(individual.genomeList, genome, rtol=rtol, atol=atol, binary_parameters=binary_parameters, return_comparison=True)
                if equal: return generation_index, key, parameters, original_parameters

            # Just return generation index and individual key
            else:

                # If the genomes are equal, return the generation index and the individual key
                if equal_genomes(individual.genomeList, genome, rtol=rtol, atol=atol, binary_parameters=binary_parameters): return generation_index, key

    # Nothing found
    if return_comparison: return None, None, None, None
    else: return None, None

# -----------------------------------------------------------------
