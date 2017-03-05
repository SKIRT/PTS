#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.genetic Contains the GeneticModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....evolve.core.engine import GeneticEngine
from .generator import ModelGenerator
from ....core.tools import filesystem as fs
from ....core.tools.random import save_state, load_state
from ....evolve.optimize.stepwise import StepWiseOptimizer

# -----------------------------------------------------------------

class GeneticModelGenerator(ModelGenerator):

    """
    This function ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).__init__(config, interactive)

        # The scores (only if this is not the initial generation)
        self.scores = None
        self.scores_check = None

        # The optimizer
        self.optimizer = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This funtion ...
        :return:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).setup(**kwargs)

        # Get the fitting run
        self.fitting_run = kwargs.pop("fitting_run")

        # Re-invoke existing optimizer run
        if fs.is_file(self.fitting_run.main_engine_path): self.optimizer = StepWiseOptimizer.from_paths(self.fitting_run.path,
                                                                                                        self.fitting_run.main_engine_path,
                                                                                                        self.fitting_run.main_prng_path,
                                                                                                        self.fitting_run.optimizer_config_path,
                                                                                                        self.statistics_path, self.database_path)

        # New optimizer run
        else:

            # Create a new optimizer and set paths
            self.optimizer = StepWiseOptimizer()
            self.optimizer.config.output = self.fitting_run.path
            self.optimizer.config.writing.engine_path = self.fitting_run.main_engine_path
            self.optimizer.config.writing.prng_path = self.fitting_run.main_prng_path
            self.optimizer.config.writing.config_path = self.fitting_run.optimizer_config_path
            self.optimizer.config.writing.statistics_path = self.statistics_path
            self.optimizer.config.writing.database_path = self.database_path

        # Set settings
        self.set_optimizer_settings()

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This function ...
        :return:
        """

        ## In order of optimizer configuration

        # User
        self.optimizer.config.mutation_rate = self.fitting_run.genetic_settings.mutation_rate
        self.optimizer.config.crossover_rate = self.fitting_run.genetic_settings.crossover_rate

        # Fixed
        self.optimizer.config.stats_freq = 1
        self.optimizer.config.best_raw_score = 0.

        # User
        self.optimizer.config.rounddecimal = self.fitting_run.genetic_settings.rounddecimal
        self.optimizer.config.mutation_method = self.fitting_run.genetic_settings.mutation_method

        # Fixed
        self.optimizer.config.min_or_max = "minimize"
        self.optimizer.config.run_id = self.fitting_run.name
        self.optimizer.config.database_frequency = 1
        self.optimizer.config.statistics_frequency = 1

        # Fixed
        #self.optimizer.config.output = self.fitting_run.path

        # Fixed
        self.optimizer.config.elitism = True

        # Fixed
        self.optimizer.config.nelite_individuals = self.fitting_run.genetic_settings.nelite_individuals

    # -----------------------------------------------------------------

    def setup_from_initial(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(InitialModelGenerator, self).setup(**kwargs)

        # Create the first genome
        genome = G1DList(self.nparameters)

        # Set genome options
        genome.setParams(minima=self.parameter_minima, maxima=self.parameter_maxima, bestrawscore=0.00, rounddecimal=2)
        genome.initializator.set(initializators.HeterogeneousListInitializerReal)
        genome.mutator.set(mutators.HeterogeneousListMutatorRealRange)

        # If gaussian is used
        #genome.mutator.set(mutators.HeterogeneousListMutatorRealGaussian)
        #genome.setParams(centers=centers, sigmas=sigmas) # means is a list with the centers of the gaussian function [for each parameter], sigmas is a list with the standard deviations [for each parameter]

        # Create the genetic algorithm engine
        self.engine = GeneticEngine(genome)

        # Set options for the engine
        self.engine.terminationCriteria.set(RawScoreCriteria)
        self.engine.setMinimax(constants.minimaxType["minimize"])
        self.engine.setGenerations(5) # not important in this case
        self.engine.setCrossoverRate(self.config.crossover_rate)
        self.engine.setPopulationSize(self.config.nmodels)
        self.engine.setMutationRate(self.config.mutation_rate)

        # Initialize the genetic algorithm
        self.engine.initialize()

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores
        self.set_scores()

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check)

        # Get the parameter values of the new models
        self.get_model_parameters()

    # -----------------------------------------------------------------

    def generate_from_initial(self):

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
                self.parameters[self.free_parameter_labels[i]].append(value)

    # -----------------------------------------------------------------

    def set_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting scores from previous generation ...")

        # Load the parameters table from the previous generation
        parameters_table = self.fitting_run.parameters_table_for_generation(self.fitting_run.last_genetic_or_initial_generation_name)

        # Load the chi squared table from the previous generation
        chi_squared_table = self.fitting_run.chi_squared_table_for_generation(self.fitting_run.last_genetic_or_initial_generation_name)

        # List of chi squared values in the same order as the parameters table
        chi_squared_values = []

        # Check whether the chi-squared and parameter tables match
        for i in range(len(parameters_table)):
            simulation_name = parameters_table["Simulation name"][i]
            chi_squared = chi_squared_table.chi_squared_for(simulation_name)
            chi_squared_values.append(chi_squared)

        # Get the scores
        scores = chi_squared_table["Chi squared"]

        # Check individual values with parameter table of the last generation
        check = []
        for label in self.fitting_run.free_parameter_labels:
            values = parameters_table[label]
            check.append(values)

        # Set the scores
        self.scores = scores
        self.scores_check = check

    # -----------------------------------------------------------------

    def get_model_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model parameters ...")

        # Loop over the individuals of the population
        for individual in self.optimizer.population:

            # Loop over all the genes (parameters)
            for i in range(len(individual)):

                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.fitting_run.free_parameter_labels[i]].append(value)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the state of the optimizer
        #self.write_optimizer()

# -----------------------------------------------------------------
