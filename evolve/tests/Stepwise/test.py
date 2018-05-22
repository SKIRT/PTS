#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.range import RealRange
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log
from pts.core.tools.loops import repeat_check
from pts.evolve.optimize.stepwise import StepWiseOptimizer
from pts.core.tools import types
from pts.evolve.optimize.tables import ScoresTable
from pts.modeling.fitting.tables import GenerationsTable, ParametersTable
from pts.modeling.fitting.explorer import GenerationInfo
from pts.core.tools import stringify
from pts.modeling.fitting.tables import BestParametersTable
from pts.core.tools import time, tables
from pts.evolve.analyse.database import get_scores, get_scores_named_individuals, load_database
from pts.evolve.analyse.statistics import get_best_score_for_generation, load_statistics
from pts.core.tools import sequences
from pts.modeling.fitting.tables import IndividualsTable
from pts.core.tools.stringify import tostr
from pts.do.commandline import Command
from pts.evolve.optimize.tables import ElitismTable
from pts.evolve.analyse.database import get_best_individual_key_all_generations

# -----------------------------------------------------------------

description = "finding the maximum of the function defined by Charbonneau (1995) using the StepWiseOptimizer"

# -----------------------------------------------------------------

# Charbonneau parameter 'n'
default_n = 9

# -----------------------------------------------------------------

# Define properties
nparameters = 2
parameter_range = RealRange(0., 1.)
#best_raw_score = float('inf')
best_raw_score = 100
min_or_max = "maximize"

# -----------------------------------------------------------------

run_name = "run_0"

# -----------------------------------------------------------------

# labels and ranges
free_parameter_labels = ["x", "y"]
free_parameter_ranges = {"x": parameter_range, "y": parameter_range}

# Parameter units
parameter_units = dict()

# -----------------------------------------------------------------

# Real best values
real_best_values = {"x": 0.5, "y": 0.5}

# -----------------------------------------------------------------

class StepWiseTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(StepWiseTest, self).__init__(*args, **kwargs)

        # Paths
        self.main_engine_path = None
        self.main_prng_path = None
        self.optimizer_config_path = None
        self.generations_path = None
        self.generations_table_path = None
        self.plot_path = None

        # The database path
        self.database_path = None

        # The statistics path
        self.statistics_path = None

        # The optimizer
        self.optimizer = None

        # Initial flag
        self.initial = None

        # The parameter ranges
        self.ranges = dict()

        # The scores and check
        self.scores = None
        self.scores_check = None

        # The generation info
        self.generation_info = None

        # The generations table
        self.generations_table = None

        # Set the path to the best parameters table
        self.best_parameters_table_path = None

        # The best parameters table
        self.best_parameters_table = None

        ## Per generation

        # The individual names
        self.individual_names = None

        # The scores table
        self.scores_table = None

        # The individuals table
        self.individuals_table = None

        # The parameters table
        self.parameters_table = None

        # Paths
        self.generation_path = None
        self.individuals_table_path = None
        self.parameters_table_path = None
        self.scores_table_path = None

        # The dictionary with the list of the model parameters
        #self.parameters = defaultdict(list)
        self.parameters = defaultdict(dict)

        ## END

        # The final population
        self.population = None

        # The best individual
        self.best = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Optimize
        self.optimize()

        # 3. Get best parameter values
        self.get_best_parameter_values()

        # 3. Writing
        self.write()

        # 4. Plotting
        if self.config.plot: self.plot()

        # 5. Test
        self.test()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(StepWiseTest, self).setup(**kwargs)

        # Set paths
        self.main_engine_path = fs.join(self.path, "engine.pickle")
        self.main_prng_path = fs.join(self.path, "prng.pickle")
        self.optimizer_config_path = fs.join(self.path, "optimizer.cfg")
        self.generations_path = fs.create_directory_in(self.path, "generations")
        self.generations_table_path = fs.join(self.path, "generations.dat")
        self.plot_path = fs.create_directory_in(self.path, "plot")

        # Set the path to the database
        self.database_path = fs.join(self.path, "database.db")

        # Set the path to the statistics file
        self.statistics_path = fs.join(self.path, "statistics.csv")

        # Set the path to the best parameters table
        self.best_parameters_table_path = fs.join(self.path, "best_parameters.dat")

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Optimizing ...")

        # 1. Start: launch the initial generation
        self.start()

        # 2. Advance: launch generations 0 -> (n-1)
        repeat_check(self.advance, self.config.ngenerations)

        # 3. Finish
        self.finish()

    # -----------------------------------------------------------------

    def start(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Starting with a randomly created generation ...")

        # Create the generations table
        self.create_generations_table()

        # Create the best parameters table
        self.create_best_parameters_table()

        # Explore
        self.explore()

    # -----------------------------------------------------------------

    def create_generations_table(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the generations table ...")

        # Create the table
        self.generations_table = GenerationsTable(parameters=free_parameter_labels, units=parameter_units)

    # -----------------------------------------------------------------

    def create_best_parameters_table(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the best parameters table ...")

        # Create the table
        self.best_parameters_table = BestParametersTable(parameters=free_parameter_labels, units=parameter_units)

    # -----------------------------------------------------------------

    @property
    def individual_names_dict(self):

        """
        This function ...
        :return: 
        """

        individual_names = dict()
        for generation_name in self.generations_table.generation_names:
            individual_names[generation_name] = self.parameters_table.simulation_names
        return individual_names

    # -----------------------------------------------------------------

    def advance(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Advancing the fitting with a new generation ...")

        # Check whether there is a generation preceeding this one
        if self.generations_table.last_generation_name is None: raise RuntimeError("Preceeding generation cannot be found")

        # Debugging
        log.debug("Previous generation: " + self.generations_table.last_generation_name)

        # Debugging
        if self.generations_table.has_finished: log.debug("There are finished generations: " + tostr(self.generations_table.finished_generations))
        if has_unevaluated_generations(self.generations_table, self.best_parameters_table): log.debug("There are unevaluated generations: " + tostr(get_unevaluated_generations(self.generations_table, self.best_parameters_table)))

        # If some generations have finished, fit the SED
        if self.generations_table.has_finished and has_unevaluated_generations(self.generations_table, self.best_parameters_table): self.score()

        # If all generations have finished, explore new generation of models
        if self.generations_table.all_finished: self.explore()

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Evaluating last generation ...")

        # Check the current number of generations
        current_ngenerations = self.generations_table.ngenerations
        if current_ngenerations <= 1: raise RuntimeError("Need at least one generation after the initial generation to finish the fitting")
        if current_ngenerations == 0: raise RuntimeError("There are no generations")

        # Check if there are unfinished generations
        has_unfinished = has_unfinished_generations(self.generations_table)
        if has_unfinished: log.warning("There are unfinished generations, but evaluating finished simulations anyway ...")

        # Check if there are unevaluated generations
        if not has_unevaluated_generations(self.generations_table, self.best_parameters_table): log.success("All generations have already been evaluated")

        # Do the SED fitting step
        self.score()

        # Set the scores to the optimizer
        self.finish_optimizer()

        # Success
        if not has_unfinished: log.success("Succesfully evaluated all generations")

    # -----------------------------------------------------------------

    def finish_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Finishing the optimizer ...")

        # 1. Load the optimizer
        self.load_optimizer()

        # 2. Set settings
        self.set_optimizer_settings()

        # 3. Set finish flag
        self.optimizer.config.finish = True

        # 4. Run the optimizer
        self.run_optimizer(get_parameters=False)

        # 5. Set the population
        self.population = self.optimizer.population

        # 6. Get the best individual
        self.best = self.optimizer.best

        # Debugging
        log.debug("Best individual:")
        if log.is_debug: print(self.best)

    # -----------------------------------------------------------------

    def explore(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Exploring the parameter space ...")

        # 1. Reset the generation
        self.reset_generation()

        # 2. Create the generation info
        self.create_generation_info()

        # 3. Create the generation directory
        self.create_generation_directory()

        # 4. Generate the parameters
        self.generate()

        # 5. Add the generation to the table (after generate because ranges have to be set)
        self.add_generation()

        # Fill individuals and parameters tables
        self.fill_individuals_and_parameters_tables()

        # 6. Evaluate
        self.evaluate()

        # 7. Write exploration output
        self.write_exploration()

    # -----------------------------------------------------------------

    @property
    def generation_names(self):

        """
        This function ...
        :return: 
        """

        return self.generations_table.generation_names

    # -----------------------------------------------------------------

    def reset_generation(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Resetting the generation ...")

        # Set the generation info to None
        self.generation = None

        # Set the individual names list to None
        self.individual_names = None

        # Set paths to None
        self.generation_path = None
        self.individuals_table_path = None
        self.parameters_table_path = None
        self.scores_table_path = None

        # Parameters table
        self.parameters_table = None

        # Individuals table
        self.individuals_table = None

        # Scores table
        self.scores_table = None

        # The dictionary with the list of the model parameters
        self.parameters = defaultdict(dict)

    # -----------------------------------------------------------------

    def get_generation_name(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        return str("Generation" + str(index))

    # -----------------------------------------------------------------

    def create_generation_info(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the generation info ...")

        # Not the initial generation
        if "initial" in self.generation_names:

            # Set index and name
            generation_index = self.last_genetic_generation_index + 1
            generation_name = self.get_generation_name(generation_index)

        # Initial generation
        else:

            generation_index = None
            generation_name = "initial"

        # Create the generation info object
        self.generation_info = GenerationInfo()

        # Set the generation info
        self.generation_info.name = generation_name
        self.generation_info.index = generation_index
        self.generation_info.method = "genetic"
        self.generation_info.wavelength_grid_level = None
        #self.generation_info.nsimulations = self.config.nindividuals # NO: DO IT AFTER THE MODELS AR GENERATED TO GET THE ACTUAL NUMBER (recurrence)
        self.generation_info.npackages = None
        self.generation_info.selfabsorption = None
        self.generation_info.transient_heating = None

    # -----------------------------------------------------------------

    def add_generation(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Adding the generation to the table ...")

        # Add an entry to the generations table
        self.generations_table.add_entry(self.generation, self.ranges)

    # -----------------------------------------------------------------

    def path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.generations_path, generation_name)

    # -----------------------------------------------------------------

    def individuals_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.path_for_generation(generation_name), "individuals.dat")

    # -----------------------------------------------------------------

    def parameters_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.path_for_generation(generation_name), "parameters.dat")

    # -----------------------------------------------------------------

    def scores_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.path_for_generation(generation_name), "scores.dat")

    # -----------------------------------------------------------------

    def elitism_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.path_for_generation(generation_name), "elitism.dat")

    # -----------------------------------------------------------------

    def create_generation_directory(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the generation directory ...")

        # Create
        self.generation_path = self.path_for_generation(self.generation_info.name)
        fs.create_directory(self.generation_path)

        # Set paths
        self.individuals_table_path = self.individuals_table_path_for_generation(self.generation_info.name)
        self.parameters_table_path = self.parameters_table_path_for_generation(self.generation_info.name)
        self.scores_table_path = self.scores_table_path_for_generation(self.generation_info.name)

        # Initialize the individuals table
        self.individuals_table = IndividualsTable()
        self.individuals_table.saveto(self.individuals_table_path)

        # Initialize the parameters table
        self.parameters_table = ParametersTable(parameters=free_parameter_labels, units=parameter_units)
        self.parameters_table.saveto(self.parameters_table_path)

        # Initialize the scores table
        self.scores_table = ScoresTable(min_or_max="max")
        self.scores_table.saveto(self.scores_table_path)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the parameter values ...")

        # 1. Reset
        self.reset_optimizer()

        # 2. Setup
        self.set_optimizer()

        # 3. Set parameters
        self.set_parameters()

        # 4. Run the optimizer
        self.run_optimizer()

    # -----------------------------------------------------------------

    def reset_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Resetting the optimizer ...")

        # Set to None
        self.optimizer = None

    # -----------------------------------------------------------------

    def set_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the optimizer ...")

        # Load or create the optimizer
        if fs.is_file(self.main_engine_path): self.load_optimizer()
        else: self.create_optimizer()

        # Set settings
        self.set_optimizer_settings()

    # -----------------------------------------------------------------

    def load_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the optimizer ...")

        # Load the optimizer from files
        self.optimizer = StepWiseOptimizer.from_paths(self.path, self.main_engine_path, self.main_prng_path, self.optimizer_config_path, self.statistics_path, self.database_path, run_name)

        # Set initial flag
        self.initial = False

    # -----------------------------------------------------------------

    def create_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the optimizer ...")

        # Create a new optimizer and set paths
        self.optimizer = StepWiseOptimizer()
        self.optimizer.config.output = self.path
        self.optimizer.config.writing.engine_path = self.main_engine_path
        self.optimizer.config.writing.prng_path = self.main_prng_path
        self.optimizer.config.writing.config_path = self.optimizer_config_path
        self.optimizer.config.writing.statistics_path = self.statistics_path
        self.optimizer.config.writing.database_path = self.database_path
        self.optimizer.config.run_id = run_name

        # Set initial flag
        self.initial = True

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting optimizer settings ...")

        ## In order of optimizer configuration

        # Parameters
        self.optimizer.config.nparameters = nparameters

        # Number of generations and the number of individuals per generation
        self.optimizer.config.ngenerations = self.config.ngenerations
        self.optimizer.config.nindividuals = self.config.nindividuals

        # User
        self.optimizer.config.mutation_rate = self.config.genetic.mutation_rate
        self.optimizer.config.crossover_rate = self.config.genetic.crossover_rate

        # Fixed
        self.optimizer.config.stats_freq = 1
        self.optimizer.config.best_raw_score = 0.

        # User
        self.optimizer.config.round_decimal = self.config.genetic.round_decimal
        self.optimizer.config.mutation_method = self.config.genetic.mutation_method

        # User, scaling
        self.optimizer.config.scaling_method = self.config.genetic.scaling_method

        # User, selector
        self.optimizer.config.selector_method = self.config.genetic.selector_method

        # Fixed
        self.optimizer.config.min_or_max = min_or_max
        self.optimizer.config.database_frequency = 1
        self.optimizer.config.statistics_frequency = 1

        # Fixed
        #self.optimizer.config.output = self.fitting_run.path

        # Fixed
        self.optimizer.config.elitism = self.config.genetic.elitism

        # Fixed
        self.optimizer.config.nelite_individuals = self.config.genetic.nelite_individuals

        # Set heterogeneous flag
        #self.optimizer.config.heterogeneous = True
        self.optimizer.config.heterogeneous = False

        # Set named individuals
        self.optimizer.config.named_individuals = True

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Loop over the free parameters
        for label in free_parameter_labels:

            # Get range
            parameter_range = free_parameter_ranges[label]

            # Add the parameter
            self.add_parameter(label, parameter_range)

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in free_parameter_labels: minima.append(self.ranges[label].min)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in free_parameter_labels: maxima.append(self.ranges[label].max)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    @property
    def parameter_minima_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in free_parameter_labels:

            min_value = self.ranges[label].min

            # Convert if necessary
            if label in parameter_units and parameter_units[label] is not None:
                unit = parameter_units[label]
                min_value = min_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(min_value)
            min_value = float(min_value)

            # Add to list
            minima.append(min_value)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in free_parameter_labels:

            max_value = self.ranges[label].max

            # Convert if necessary
            if label in parameter_units and parameter_units[label] is not None:
                unit = parameter_units[label]
                max_value = max_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(max_value)
            max_value = float(max_value)

            # Add to list
            maxima.append(max_value)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    def add_parameter(self, label, parameter_range):

        """
        This function ...
        :param label:
        :param parameter_range:
        :return:
        """

        self.ranges[label] = parameter_range

    # -----------------------------------------------------------------

    @property
    def nmodels(self):
        """
        This function ...
        :return:
        """

        return len(self.parameters[free_parameter_labels[0]])

    # -----------------------------------------------------------------

    def run_optimizer(self, get_parameters=True):

        """
        This function ...
        :param get_parameters:
        :return: 
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores
        if not self.initial: self.set_scores()

        evaluator = None
        evaluator_kwargs = None

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check, minima=self.parameter_minima_scalar,
                           maxima=self.parameter_maxima_scalar, evaluator=evaluator, evaluator_kwargs=evaluator_kwargs)

        # Get the parameter values of the new models
        if get_parameters: self.get_model_parameters()

        # Set the number of individuals for this generation
        self.generation_info.nsimulations = self.nmodels

    # -----------------------------------------------------------------

    @property
    def last_genetic_generation_name(self):

        """
        This function ...
        :return: 
        """

        highest_index = -1
        name = None

        # Find the name of the generation with the highest index
        for i in range(len(self.generations_table)):
            if not self.generations_table["Generation index"].mask[i]:
                index = self.generations_table["Generation index"][i]
                if index > highest_index:
                    highest_index = index
                    name = self.generations_table["Generation name"][i]

        # Return the name of the generation with the highest index
        return name

    # -----------------------------------------------------------------

    @property
    def last_genetic_or_initial_generation_name(self):

        """
        This function ...
        :return: 
        """

        name = self.last_genetic_generation_name

        # Check whether the initial generation exists
        if name is None and "initial" in self.generations_table["Generation name"]: name = "initial"

        # Return the name
        return name

    # -----------------------------------------------------------------

    @property
    def last_genetic_generation_index(self):

        """
        This function ...
        :return:
        """

        highest_index = -1

        # Find the highest index
        for i in range(len(self.generations_table)):
            if not self.generations_table["Generation index"].mask[i]:
                index = self.generations_table["Generation index"][i]
                if index  > highest_index: highest_index = index

        # Return the highest generation index
        return highest_index

    # -----------------------------------------------------------------

    def individuals_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return IndividualsTable.from_file(self.individuals_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def parameters_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return ParametersTable.from_file(self.parameters_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def scores_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return ScoresTable.from_file(self.scores_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def elitism_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return ElitismTable.from_file(self.elitism_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def set_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting scores from previous generation ...")

        # Load the parameters table from the previous generation
        parameters_table = self.parameters_table_for_generation(self.last_genetic_or_initial_generation_name)

        # Load the chi squared table from the previous generation
        scores_table = self.scores_table_for_generation(self.last_genetic_or_initial_generation_name)

        # List of chi squared values in the same order as the parameters table
        score_values = []

        # Check whether the chi-squared and parameter tables match
        for i in range(len(parameters_table)):
            individual_name = parameters_table["Simulation name"][i]
            score = scores_table.score_for(individual_name)
            score_values.append(score)

        # Check individual values with parameter table of the last generation
        check = []
        for label in free_parameter_labels:
            values = parameters_table[label]
            check.append(values)

        # Set the scores
        self.scores = score_values
        self.scores_check = check

    # -----------------------------------------------------------------

    def get_model_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model parameters ...")

        # THIS WAS THE IMPLEMENTATION BEFORE RECURRENCE CHECKING:

        # Loop over the individuals of the population
        #for individual in self.optimizer.population:

        # Set the individual names
        #self.individual_names = self.optimizer.population.names

        # Loop over the individual names
        #for name in self.individual_names:

            # Get the individual
            #individual = self.optimizer.population[name]

            # Loop over all the genes (parameters)
            #for i in range(len(individual)):

                # Get the parameter value
                #value = individual[i]

                # Add the parameter value to the dictionary
                #self.parameters[free_parameter_labels[i]].append(value)

                #self.parameters[free_parameter_labels[i]][name] = value

        # NEW IMPLEMENTATION: AS IN GENETICMODELGENERATOR

        # Loop over the individual names
        for name in self.individual_names:

            # Get the parameters
            for index, parameter_value in enumerate(self.optimizer.get_parameters(name)):

                # Add the parameter value to the dictionary
                self.parameters[free_parameter_labels[index]][name] = parameter_value

    # -----------------------------------------------------------------

    def fill_individuals_and_parameters_tables(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Filling the individuals and parameters tables ...")

        counter = 0

        # Loop over the individual names
        for name in self.individual_names:

            # Get the parameter values
            parameter_values = get_parameter_values_for_named_individual(self.parameters, name)

            # Generate the 'simulation' name
            simulation_name = "individual" + str(counter)

            # Debugging
            log.debug("Adding an entry to the individuals table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            log.debug(" - Individual_name: " + name)
            log.debug("")

            # Add entry
            self.individuals_table.add_entry(simulation_name, name)

            # Debugging
            log.debug("Adding entry to the parameters table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            for label in parameter_values: log.debug(" - " + label + ": " + stringify.stringify_not_list(parameter_values[label])[1])
            log.debug("")

            # Add an entry to the parameters table
            self.parameters_table.add_entry(simulation_name, parameter_values)

            # Increment counter
            counter += 1

    # -----------------------------------------------------------------

    def evaluate(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Evaluating the generation ...")

        # Loop over the individuals
        for name in self.individual_names:

            # Determine the simulation name
            simulation_name = self.individuals_table.get_simulation_name(name)

            # Find the index in the table for this generation
            index = tables.find_index(self.generations_table, self.generation_info.name, "Generation name")

            # Get the number of simulations for this generation
            nsimulations = self.generations_table["Number of simulations"][index]

            # Inform the user
            log.info("Calculating the value ...")

            # Get the parameter values
            parameter_values = self.parameters_table.parameter_values_for_simulation(simulation_name)

            # Calculate the value
            value = eval_func_xy(parameter_values["x"], parameter_values["y"])

            # Inform the user
            log.info("Adding the score for the current model to the score table ...")

            # Debugging
            log.debug("The score for this individual is " + str(value))

            # Add entry
            self.scores_table.add_entry(simulation_name, value)

            # Get the number of entries in the chi squared table
            nfinished_simulations = len(self.scores_table)

            # If this is the last simulation
            if nsimulations == nfinished_simulations + 1: self.generations_table.set_finishing_time(self.generation_info.name, time.timestamp())

    # -----------------------------------------------------------------

    def write_exploration(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing after exploration step ...")

        # 1. Write individuals
        self.write_individuals()

        # 2. Write parameters
        self.write_parameters()

        # 3. Write scores
        self.write_scores()

    # -----------------------------------------------------------------

    def write_individuals(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the individuals table ...")

        # Save the individuals table
        self.individuals_table.saveto(self.individuals_table_path)

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model parameters table ...")

        # Save the parameters table
        self.parameters_table.saveto(self.parameters_table_path)

    # -----------------------------------------------------------------

    def write_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the scores table ...")

        # Save the chi squared table
        self.scores_table.saveto(self.scores_table_path)

    # -----------------------------------------------------------------

    @property
    def finished_generations(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.finished_generations

    # -----------------------------------------------------------------

    def score(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Scoring engine ...")

        # Check if there are finished generations
        if len(self.finished_generations) == 0: raise RuntimeError("There are no finished generations")

        # 2. Get the parameters of the best models for each generation
        self.get_best_parameters()

    # -----------------------------------------------------------------

    def is_finished_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.generations_table.is_finished(generation_name)

    # -----------------------------------------------------------------

    def best_parameter_values_for_generation(self, generation_name, return_chi_squared=False, only_finished=True):

        """
        This function ...
        :param generation_name:
        :param return_chi_squared:
        :param only_finished:
        :return:
        """

        # Check if the generation is finished (if this is required by the caller)
        if only_finished:
            if not self.is_finished_generation(generation_name): raise RuntimeError("The generation '" + generation_name + "' is not yet finished")

        # Open the chi squared table
        #chi_squared_table = self.chi_squared_table_for_generation(generation_name)
        scores_table = self.scores_table_for_generation(generation_name)

        # Get the name of the simulation with the lowest chi squared value
        #best_simulation_name = chi_squared_table.best_simulation_name
        best_individual_name = scores_table.best_individual_name

        # Open the parameters table for this generation
        parameters_table = self.parameters_table_for_generation(generation_name)

        # Return the parameters of the best simulation
        if return_chi_squared:
            chi_squared = scores_table.score_for(best_individual_name)
            return parameters_table.parameter_values_for_simulation(best_individual_name), chi_squared
        else: return parameters_table.parameter_values_for_simulation(best_individual_name)

    # -----------------------------------------------------------------

    def get_best_parameters(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Getting the parameter values of the best model for the finished generations (if not already done) ...")

        # Loop over the finished generations
        for generation_name in self.finished_generations:

            # Check if the generation is already in the best parameters table
            if generation_name in self.best_parameters_table.generation_names: continue

            # Debugging
            log.info("Getting the best parameter values for generation '" + generation_name + "' ...")

            # Otherwise, add the best parameter values
            values, chi_squared = self.best_parameter_values_for_generation(generation_name, return_chi_squared=True)

            # Add an entry to the best parameters table file
            self.best_parameters_table.add_entry(generation_name, values, chi_squared)

    # -----------------------------------------------------------------

    def get_best_parameter_values(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the best parameter values ...")

        # Get the best parameter values
        self.best_parameter_values = self.modeler.modeler.fitter.fitting_run.best_parameter_values

        # Debugging
        log.debug("The best parameter values are:")
        log.debug("")
        for parameter_name in self.best_parameter_values: log.debug(" - " + parameter_name + ": " + tostr(self.best_parameter_values[parameter_name], scientific=True, fancy=True, ndigits=parameter_ndigits[parameter_name]))
        log.debug("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the generations table
        self.write_generations_table()

    # -----------------------------------------------------------------

    def write_generations_table(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the generations table ...")

        # Save
        self.generations_table.saveto(self.generations_table_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the function
        self.plot_function()

        # Plot the scores
        self.plot_scores()

        # Plot the heatmap
        self.plot_heatmap()

    # -----------------------------------------------------------------

    def plot_function(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting the function ...")

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #jet = plt.get_cmap('jet')
        cmap = plt.get_cmap("viridis")

        # Determine range and number of samples
        nxsamples = 200
        nysamples = 200
        min_value = parameter_range.min
        max_value = parameter_range.max

        x = np.linspace(min_value, max_value, nxsamples)
        y = np.linspace(min_value, max_value, nysamples)

        # Calculate data
        X, Y = np.meshgrid(x, y)
        Z = eval_func_xy(X, Y)
        max_z = np.max(Z)

        # Plot the surface
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmap, linewidth=0)
        ax.set_zlim3d(0, max_z)

        # Show the
        #plt.show()

        path = fs.join(self.plot_path, "charbonneau.pdf")

        plt.savefig(path)
        plt.close()

    # -----------------------------------------------------------------

    @property
    def statistics(self):

        """
        THis function ...
        :return: 
        """

        return load_statistics(self.statistics_path)

    # -----------------------------------------------------------------

    @property
    def database(self):

        """
        This function ...
        :return: 
        """

        return load_database(self.database_path)

    # -----------------------------------------------------------------

    def plot_scores(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting the scores ...")

        # Settings
        settings = dict()
        settings["fitness"] = False
        settings["output"] = self.plot_path

        # Input
        input_dict = dict()
        input_dict["database"] = self.database

        # Plot
        command = Command("plot_scores", "plot the scores", settings, input_dict)
        plotter = self.run_command(command)

    # -----------------------------------------------------------------

    def plot_heatmap(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting heatmap ...")

        # Settings
        settings = dict()
        settings["fitness"] = False
        settings["output"] = self.plot_path

        # Input
        input_dict = dict()
        input_dict["database"] = self.database

        # Plot
        command = Command("plot_heat_map", "plot heatmap", settings, input_dict)
        plotter = self.run_command(command)

    # -----------------------------------------------------------------

    def test(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing ...")

        # Check the best value
        self.check_best()

        # Check the database
        self.check_database()

        # Check the statistics
        self.check_statistics()

    # -----------------------------------------------------------------

    def check_best(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the best individual ...")

        # Get the dictionary of the best parameter values from the optimizer
        best_parameters = get_parameter_values_from_individual(self.best)

        # Get the dictionary of the best parameter values from the tables
        best_parameters_table = self.get_best_parameters_from_tables()

        # Loop over the free parameter labels
        print("")
        for label in free_parameter_labels:

            print(label + ":")
            print("")

            real = real_best_values[label]
            best = best_parameters[label]
            abs_diff = abs(best_parameters[label] - real_best_values[label])
            rel_diff = abs_diff / real
            best_tables = best_parameters_table[label]

            # Show
            print(" - real value: " + tostr(real))
            print(" - best individual value: " + tostr(best))
            print(" - absolute difference: " + tostr(abs_diff))
            print(" - relative difference: " + tostr(rel_diff) + " (" + tostr(rel_diff*100) + "%)")
            print(" - best value from tables: " + tostr(best_tables))

            print("")

    # -----------------------------------------------------------------

    def get_best_parameters_from_tables(self):

        """
        Check whether we find the same best individual from the scores table
        :return: 
        """

        # Get the scores table of the last generation
        scores_table = self.scores_table_for_generation(self.last_genetic_generation_name)

        # Get the name of the 'simulation' with the best score
        simulation_name = scores_table.best_individual_name

        # Get the corresponding parameters from the parameters table of the last generation
        parameters_table = self.parameters_table_for_generation(self.last_genetic_generation_name)

        # Get the parameters
        best_parameters = parameters_table.parameter_values_for_simulation(simulation_name)

        # Return the parameters
        return best_parameters

    # -----------------------------------------------------------------

    def check_database(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the database ...")

        # Loop over the generations
        for index in range(self.config.ngenerations):

            # Determine the generation name
            generation_name = self.get_generation_name(index)

            # Load the scores table
            scores_table = self.scores_table_for_generation(generation_name)

            # Load the individuals table
            individuals_table = self.individuals_table_for_generation(generation_name)

            # Get the scores from the database
            database_index = index + 1
            scores_database = get_scores_named_individuals(self.database_path, run_name, database_index)

            # Keep track of the mismatches
            mismatches = []

            # Loop over the indnividual names
            for name in scores_database:

                # Get the score from the database
                score_database = scores_database[name]

                # Get the simulation name
                simulation_name = individuals_table.get_simulation_name(name)

                # Get the score
                score = scores_table.score_for(simulation_name)

                # Check if equal
                equal = np.isclose(score, score_database)
                if not equal: mismatches.append((score, score_database))

                # Report
                if len(mismatches) == 0: log.success(generation_name + ": OK")
                else:
                    log.error(generation_name + ": " + str(len(mismatches)) + " mismatch(es):")
                    for score, score_database in mismatches: log.error("   " + str(score) + " , " + str(score_database))

    # -----------------------------------------------------------------

    def check_database_unnamed(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the database ...")

        # Loop over the generations
        for index in range(self.config.ngenerations):

            # Determine the generation name
            generation_name = self.get_generation_name(index)

            # Load the parameters table
            parameters = self.parameters_table_for_generation(generation_name)

            # Load the scores table
            scores_table = self.scores_table_for_generation(generation_name)

            # Sort scores
            scores_table.sort_as(parameters.simulation_names)

            # Get individual names and scores
            individual_names = scores_table.individual_names
            scores = scores_table.scores

            # Get the scores from the database
            scores_database = get_scores(self.database_path, run_name, index)

            # Check if the lists of scores contain the same elements
            if sequences.contains_same_elements(scores, scores_database):

                # Success
                log.success(generation_name + ": elements OK")

                # Check order
                if scores == scores_database: log.success(generation_name + ": order OK")
                else: log.error(generation_name + ": order not OK")

            # The scores in both lists are not the same
            else:

                log.error(generation_name + ":")

                # First check
                first_check = sequences.elements_not_in_other(scores, scores_database)
                log.error("Scores not present in scores table from database (" + str(len(first_check)) + "):")
                print(stringify.stringify_list_fancy(first_check)[1])

                # Second check
                second_check = sequences.elements_not_in_other(scores_database, scores)
                log.error("Scores not present in database from scores table (" + str(len(second_check)) + "):")
                print(stringify.stringify_list_fancy(second_check)[1])

    # -----------------------------------------------------------------

    def get_best_score_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Get the scores table
        scores_table = self.scores_table_for_generation(generation_name)

        # Return the best score
        return scores_table.best_score

    # -----------------------------------------------------------------

    def check_statistics(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Checking the statistics file ...")

        # Loop over the generations
        print("")
        for index in range(self.config.ngenerations):

            # Determine the generation name
            generation_name = self.get_generation_name(index)

            # Get the best score from the statistics
            statistics_index = index + 1
            score_statistics = get_best_score_for_generation(self.statistics_path, run_name, statistics_index, minmax="max")

            # Get the best score from the scores table
            score_table = self.get_best_score_for_generation(generation_name)

            # Get the best score by re-evaluating the best individual
            score_best = eval_func_xy(self.best[0], self.best[1])

            # Rel diff
            rel_diff = abs(score_table - score_statistics) / score_statistics

            print(generation_name + ":")
            print("")

            if index == self.config.ngenerations - 1: print(" - Score of best individual: " + tostr(score_best))
            print(" - Best score from statistics: " + tostr(score_statistics))
            print(" - Best score from scores table: " + tostr(score_table))
            print(" - Relative difference: " + tostr(rel_diff) + " (" + tostr(rel_diff * 100) + "%)")

            print("")

# -----------------------------------------------------------------

def charbonneau(x, y, n):

    """
    This function ...
    :param x: 
    :param y: 
    :param n: 
    :return: 
    """

    # Calculate z and return
    z = (16 * x * (1. - x) * y * (1. - y) * np.sin(n * np.pi * x) * np.sin(n * np.pi * y)) ** 2
    return z

# -----------------------------------------------------------------

def eval_func_xy(x, y):

    """
    This function ...
    :param x:
    :param y:
    :param 9:
    :return:
    """

    return charbonneau(x, y, default_n)

# -----------------------------------------------------------------

def eval_func(chromosome, **kwargs):

    """
    The evaluation function
    """

    # Get x and y
    x = chromosome[0]
    y = chromosome[1]

    # Return the value of the function
    return eval_func_xy(x, y)

# -----------------------------------------------------------------

def evolve_callback(ga_engine, **kwargs):

    """
    This function ...
    :param ga_engine:
    :param kwargs:
    :return:
    """

    # Get the coords
    #coords = kwargs.pop("coordinates")

    #if ga_engine.currentGeneration % 10 == 0:
    #    best = ga_engine.bestIndividual()
    #    write_tour_to_img(coords, best, "tsp_result_%d.png" % (ga_engine.currentGeneration,))
    #return False


    #return coords, cm

    pass

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

def setup(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
# OPTIMIZE
# -----------------------------------------------------------------

# Settings
#settings_optimize = dict()
#settings_optimize["output"] = None
#settings_optimize["nparameters"] = nparameters
#settings_optimize["nindividuals"] = nindividuals
#settings_optimize["parameter_range"] = parameter_range
#settings_optimize["best_raw_score"] = best_raw_score
##settings_optimize["round_decimal"] = round_decimal
#settings_optimize["ngenerations"] = ngenerations
#settings_optimize["mutation_rate"] = mutation_rate
#settings_optimize["crossover_rate"] = crossover_rate
#settings_optimize["stats_freq"] = stats_freq
##settings_optimize["mutation_method"] = mutation_method
#settings_optimize["min_or_max"] = min_or_max

# Input
#input_optimize = dict()
##input_optimize["genome"] = genome
#input_optimize["evaluator"] = eval_func
##input_optimize["initializator"] = G1DListTSPInitializator
#input_optimize["mutator"] = G1DListMutatorSwap
#input_optimize["crossover"] = G1DListCrossoverOX
#input_optimize["callback"] = evolve_callback
##input_optimize["adapter"] = sqlite_adapter

# Create dictionary for extra arguments to the evalutor function
##input_optimize["evaluator_kwargs"] = {"distances": cm}
##input_optimize["callback_kwargs"] = {"coordinates": coords}

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    return

    # Remove 'callback' from the settings dictionary: we don't want to plot now
    ##if "callback" in input_optimize: del input_optimize["callback"]

    # Solve the problem with the original Pyevolve implementation
    #best = reference.call(settings_optimize, input_optimize)

    # Show the best individual
    #show_best(best)

# -----------------------------------------------------------------

def get_generation_names(generations_table):

    """
    This function ...
    :param generations_table:
    :return:
    """

    # Return the generation names
    return generations_table.generation_names

# -----------------------------------------------------------------

def get_unevaluated_generations(generations_table, best_parameters_table):

    """
    This function ...
    :param generations_table: 
    :param best_parameters_table:
    :return: 
    """

    generation_names = []

    # Loop over the generations
    for generation_name in get_generation_names(generations_table):
        if not best_parameters_table.has_generation(generation_name): generation_names.append(generation_name)

    # Return the generation names
    return generation_names

# -----------------------------------------------------------------

def has_unevaluated_generations(generations_table, best_parameters_table):

    """
    This function ...
    :param generations_table:
    :param best_parameters_table:
    :return:
    """

    # Loop over the generations
    for generation_name in get_generation_names(generations_table):

        # If at least one generation is not evaluated, return True
        if not best_parameters_table.has_generation(generation_name): return True

    # No generation was encountered that was not completely evaluated
    return False

# -----------------------------------------------------------------

def has_unfinished_generations(generations_table):

    """
    This function ...
    :param generations_table:
    :return:
    """

    # Open the generations table
    table = generations_table
    return table.has_unfinished

# -----------------------------------------------------------------

def get_parameter_values_for_named_individual(parameters, name):

    """
    This function ...
    :param parameters: 
    :param name: 
    :return: 
    """

    # Set the parameter values as a dictionary for this individual model
    parameter_values = dict()
    for label in free_parameter_labels:

        # Get the value for this model from the generator
        value = parameters[label][name]

        # Get unit (if any)
        unit = get_parameter_unit(label)

        # Set value with unit
        if unit is not None: parameter_values[label] = value * unit

        # Set dimensionless value
        else: parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_parameter_values_for_individual(parameters, index):

    """
    This function ...
    :param parameters:
    :param index:
    :return:
    """

    # Set the parameter values as a dictionary for this individual model
    parameter_values = dict()
    for label in free_parameter_labels:

        # Get the value for this model from the generator and get the unit defined for this parameter
        value = parameters[label][index]

        # Get unit (if any)
        unit = get_parameter_unit(label)

        # Set value with unit
        if unit is not None: parameter_values[label] = value * unit

        # Set dimensionless value
        else: parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_parameter_values_from_individual(genome):

    """
    This function ...
    :param genome: 
    :return: 
    """

    # Create a dictionary for the parameter values
    parameter_values = dict()

    # Loop over all the genes (parameters)
    for index in range(len(genome)):

        # Get the parameter value
        value = genome[index]

        # Get the label of the parameter
        label = free_parameter_labels[index]

        # Get the unit for the parameter
        unit = get_parameter_unit(label)

        # Get the unit for the parameter
        if unit is not None: parameter_values[label] = value * unit

        # Scalar parameter
        else: parameter_values[label] = value

    # Return the parameter values
    return parameter_values

# -----------------------------------------------------------------

def get_parameter_unit(label):

    """
    This function ...
    :param label:
    :return:
    """

    if label in parameter_units:
        if parameter_units[label] is None: return None
        elif parameter_units[label] == "": return None
        else: return parameter_units[label]
    else: return None

# -----------------------------------------------------------------

def get_best_parameter_values(database_path, run_name, populations):

    """
    This function ...
    :param database_path:
    :param run_name:
    :param populations:
    :return: 
    """

    # NEW: THIS FUNCTION WAS CREATED BECAUSE RECURRENCE WAS IMPLEMENTED: THIS MEANS THAT OUR OWN TABLES
    # (THOSE WHO CONTAIN ONLY MODELS THAT HAVE TO BE SIMULATED AND HAVE NOT OCCURED AND SCORED BEFORE)

    #from .component import get_statistics_path, get_populations
    #from .component import get_database_path

    # Get path
    #statistics_path = get_statistics_path(self.modeling_path)
    #database_path = get_database_path(self.modeling_path)

    # Get generation and individual
    generation_index, individual_key = get_best_individual_key_all_generations(database_path, run_name, minmax=min_or_max)

    # Look in the populations data for the parameters, for this fitting run
    #populations = get_populations(self.modeling_path)[run_name]

    # Get the parameter values for the generation and individual
    individuals_generation = populations[generation_index]

    # Get genome
    genome = individuals_generation[individual_key]

    # NEW: EXPERIMENTAL:
    # BE AWARE: IF THIS IS CHANGED, ALSO CHANGE IN OPTIMIZER -> set_nbits()
    nbits_list = []
    for index in range(len(self.ndigits_list)):
        ndigits = self.ndigits_list[index]
        low = self.parameter_minima_scalar[index]
        high = self.parameter_maxima_scalar[index]
        nbits = numbers.nbits_for_ndigits_experimental(ndigits, low, high)
        nbits_list.append(nbits)

    # Convert
    if self.genetic_settings.gray_code: parameters = gray_binary_string_to_parameters(genome, self.parameter_minima_scalar, self.parameter_maxima_scalar, nbits_list)
    else: parameters = binary_string_to_parameters(genome, self.parameter_minima_scalar, self.parameter_maxima_scalar, nbits_list)

    values = dict()

    for label_index, label in enumerate(free_parameter_labels):

        # Add unit
        if label in parameter_units:
            value = parameters[label_index] * parameter_units[label]
        else: value = parameters[label_index]

        values[label] = value

    # Return the dictionary
    return values

# -----------------------------------------------------------------
