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
from .generator import ModelGenerator
from ....core.tools import filesystem as fs
from ....evolve.optimize.stepwise import StepWiseOptimizer
from ....evolve.optimize.continuous import ContinuousOptimizer
from ..evaluate import evaluate

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

        # Whether or not this is the initial generation
        self.initial = None

        # The evaluator function and its kwargs
        self.evaluator = None
        self.evaluator_kwargs = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).setup(**kwargs)

        # Get the fitting run
        self.fitting_run = kwargs.pop("fitting_run")

        # If the number of generations in one run is more than one
        if self.config.ngenerations > 1:

            # Create a new continuous optimizer
            self.optimizer = ContinuousOptimizer()
            self.optimizer.config.output = self.fitting_run.path
            #self.optimizer.config.writing.engine_path = self.fitting_run.main_engine_path
            #self.optimizer.config.writing.prng_path = self.fitting_run.main_prng_path
            #self.optimizer.config.writing.config_path = self.fitting_run.optimizer_config_path
            #self.optimizer.config.writing.statistics_path = self.statistics_path
            #self.optimizer.config.writing.database_path = self.database_path

            self.optimizer.config.run_id = self.fitting_run.name

            # Set initial flag
            self.initial = True

            # Set the evaluator function
            self.evaluator = evaluate

            # Set the evaluator kwargs
            self.evaluator_kwargs = dict()
            self.evaluator_kwargs["fitting_run"] = self.fitting_run #self.fitting_run.name

        else:

            # Re-invoke existing optimizer run
            if fs.is_file(self.fitting_run.main_engine_path):

                # Load the optimizer from files
                self.optimizer = StepWiseOptimizer.from_paths(self.fitting_run.path,
                                                                self.fitting_run.main_engine_path,
                                                                self.fitting_run.main_prng_path,
                                                                self.fitting_run.optimizer_config_path,
                                                                self.statistics_path, self.database_path, self.fitting_run.name)
                # Set initial flag
                self.initial = False

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
                self.optimizer.config.run_id = self.fitting_run.name

                # Set initial flag
                self.initial = True

        # Set parameters
        self.set_parameters()

        # Set settings
        self.set_optimizer_settings()

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Use: self.parameter_ranges_for_generation(self, generation_name) ????

            # Get range
            parameter_range = self.fitting_run.free_parameter_ranges[label]

            # Add the parameter
            self.add_parameter(label, parameter_range)

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This function ...
        :return:
        """

        ## In order of optimizer configuration

        # Parameters
        self.optimizer.config.nparameters = self.fitting_run.nfree_parameters

        # Number of generations and the number of individuals per generation
        self.optimizer.config.ngenerations = self.config.ngenerations
        self.optimizer.config.nindividuals = self.config.nmodels

        # User
        self.optimizer.config.mutation_rate = self.fitting_run.genetic_settings.mutation_rate
        self.optimizer.config.crossover_rate = self.fitting_run.genetic_settings.crossover_rate

        # Fixed
        self.optimizer.config.stats_freq = 1
        self.optimizer.config.best_raw_score = 0.

        # User
        self.optimizer.config.round_decimal = self.fitting_run.genetic_settings.round_decimal
        self.optimizer.config.mutation_method = self.fitting_run.genetic_settings.mutation_method

        # User, scaling
        self.optimizer.config.scaling_method = self.fitting_run.genetic_settings.scaling_method

        # User, selector
        self.optimizer.config.selector_method = self.fitting_run.genetic_settings.selector_method

        # Fixed
        self.optimizer.config.min_or_max = "minimize"
        #self.optimizer.config.run_id = self.fitting_run.name # THIS IS NOW DONE IN THE SETUP
        self.optimizer.config.database_frequency = 1
        self.optimizer.config.statistics_frequency = 1

        # Fixed
        #self.optimizer.config.output = self.fitting_run.path

        # Fixed
        self.optimizer.config.elitism = True

        # Fixed
        self.optimizer.config.nelite_individuals = self.fitting_run.genetic_settings.nelite_individuals

        # Set heterogeneous flag
        self.optimizer.config.heterogeneous = True

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores
        if not self.initial: self.set_scores()

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check, minima=self.parameter_minima_scalar,
                           maxima=self.parameter_maxima_scalar, evaluator=self.evaluator, evaluator_kwargs=self.evaluator_kwargs)

        # Get the parameter values of the new models
        self.get_model_parameters()

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
