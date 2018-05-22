#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.solve.extremizer Contains the Extremizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..optimize.continuous import ContinuousOptimizer
from ...core.basics.configuration import ConfigurationDefinition
from ..config.optimize import crossover_methods, mutation_methods, scaling_methods, selector_methods, genome_types
from ...core.tools.random import skirt_seed

# -----------------------------------------------------------------

default_genome_type = "list"
default_mutation_rate = 0.03
default_crossover_rate = 0.65
default_crossover_method = "single_point"
default_round_decimal = 3
default_mutation_method = "range"
default_elitism = True
default_nelite_individuals = 1
default_scaling_method = "linear"
default_selector_method = "roulette_wheel"

default_binary_mutation_method = "flip"
binary_mutation_methods = ["flip", "swap"]

# -----------------------------------------------------------------

# Create the genetic definition
genetic_definition = ConfigurationDefinition(write_config=False)

# Random seed
genetic_definition.add_optional("seed", "positive_integer", "random seed to use for initiating the genetic algorithm", skirt_seed)

# Add settings
genetic_definition.add_optional("genome_type", "string", "genome type", default_genome_type, choices=genome_types)
genetic_definition.add_optional("mutation_rate", "real", "mutation rate", default_mutation_rate)
genetic_definition.add_optional("crossover_rate", "real", "crossover rate", default_crossover_rate)
genetic_definition.add_optional("crossover_method", "string", "crossover method", default_crossover_method, choices=crossover_methods)
genetic_definition.add_optional("round_decimal", "integer", "round everything to this decimal place", default_round_decimal)
genetic_definition.add_optional("mutation_method", "string", "mutation method", default_mutation_method, choices=mutation_methods)
genetic_definition.add_optional("binary_mutation_method", "string", "mutation method for binary string genomes", default_binary_mutation_method, choices=binary_mutation_methods)
genetic_definition.add_optional("scaling_method", "string", "scaling method", default_scaling_method, choices=scaling_methods)
genetic_definition.add_optional("selector_method", "string", "selector method", default_selector_method, choices=selector_methods)

# Flags
genetic_definition.add_flag("elitism", "enable elitism", default_elitism)

# Advanced
genetic_definition.add_optional("nelite_individuals", "positive_integer", "number of individuals to take as elite", default_nelite_individuals)
genetic_definition.add_flag("gray_code", "use Gray coding for binary genome representations", True)

# -----------------------------------------------------------------

def find_minimum(function, initial_values, ranges, genetic_settings=None):

    """
    This function ...
    :param function:
    :param initial_values:
    :param ranges:
    :param genetic_settings:
    :return: 
    """

    return find_extremum(function, initial_values, ranges, "min", genetic_settings=genetic_settings)

# -----------------------------------------------------------------

def find_maximum(function, initial_values, ranges, genetic_settings=None):

    """
    This function ...
    :param function: 
    :param initial_values: 
    :param ranges:
    :param genetic_settings:
    :return: 
    """

    return find_extremum(function, initial_values, ranges, "max", genetic_settings=genetic_settings)

# -----------------------------------------------------------------

def find_extremum(function, initial_values, ranges, min_or_max, genetic_settings=None):

    """
    This function ...
    :param function:
    :param initial_values:
    :param ranges:
    :param min_or_max:
    :param genetic_settings:
    :return: 
    """

    # Settings
    settings = dict()
    if genetic_settings is not None: settings["genetic"] = genetic_settings
    settings["min_or_max"] = min_or_max

    # Input
    input_dict = dict()
    input_dict["function"] = function
    input_dict["initial_values"] = initial_values
    input_dict["ranges"] = ranges

    # Create the extremizer
    extremizer = Extremizer(settings)

    # Run
    extremizer.run(**input_dict)

    # Return the best parameter values
    return extremizer.result

# -----------------------------------------------------------------

class Extremizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param config: 
        :param interactive: 
        """

        # Call the constructor of the base class
        super(Extremizer, self).__init__(*args, **kwargs)

        # The function to be minimized or maximized
        self.function = None

        # The initial values and ranges
        self.initial_values = None
        self.ranges = None

        # The optimizer
        self.optimizer = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 2. Optimize
        self.optimize()

        # 3. Plot
        self.plot()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(Extremizer, self).setup(**kwargs)

        # Get the function
        self.function = kwargs.pop("function")

        # Get the initial values
        self.initial_values = kwargs.pop("initial_values")

        # Get the ranges
        self.ranges = kwargs.pop("ranges")

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return: 
        """

        return len(self.initial_values)

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return: 
        """

        return [parameter_range.min for parameter_range in self.ranges]

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return: 
        """

        return [parameter_range.max for parameter_range in self.ranges]

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Running optimization ...")

        # Settings
        settings_optimize = dict()
        settings_optimize["output"] = None
        settings_optimize["nparameters"] = self.nparameters
        settings_optimize["nindividuals"] = self.config.genetic.nindividuals
        #settings_optimize["best_raw_score"] = best_raw_score
        settings_optimize["round_decimal"] = self.config.genetic.round_decimal
        settings_optimize["ngenerations"] = self.config.genetic.ngenerations
        settings_optimize["mutation_rate"] = self.config.genetic.mutation_rate
        settings_optimize["crossover_rate"] = self.config.genetic.crossover_rate
        #settings_optimize["stats_freq"] = stats_freq
        settings_optimize["mutation_method"] = self.config.genetic.mutation_method
        settings_optimize["min_or_max"] = self.config.min_or_max

        # Set random seed
        settings_optimize["seed"] = self.config.genetic.seed

        # Input
        input_optimize = dict()
        input_optimize["evaluator"] = self.function

        # Minima and maxima
        input_optimize["minima"] = self.parameter_minima
        input_optimize["maxima"] = self.parameter_maxima

        # Create the optimizer
        self.optimizer = ContinuousOptimizer(settings_optimize)

        # Run the optimizer
        self.optimizer.run(**input_optimize)

    # -----------------------------------------------------------------

    @property
    def result(self):

        """
        This function ...
        :return: 
        """

        return self.optimizer.best

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
