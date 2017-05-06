#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# CHOICES

genome_dimensions = [1, 2]

crossover_methods = ["single_point", "two_point", "uniform", "OX", "edge", "cut_crossfill", "real_SBX", "single_vertical_point", "single_horizontal_point"]
mutation_methods = ["range", "gaussian", "binary"]
min_or_max = ["minimize", "maximize"]

scaling_methods = ["linear", "sigma_truncation", "power_law", "boltzmann", "exponential", "saturated"]
selector_methods = ["rank", "uniform", "tournament", "tournament_alternative", "roulette_wheel"]

# -----------------------------------------------------------------

# DEFAULT VALUES

default_genome_type = "list"
default_dimension = 1
default_nparameters = 20
default_nindividuals = 80
default_ngenerations = 100
default_mutation_rate = 0.02
default_crossover_rate = 0.9
default_crossover_method = "single_point"
default_stats_freq = 10
default_round_decimal = 2
default_mutation_method = "range"
default_min_or_max = "maximize"
default_database_frequency = 1
default_statistics_frequency = 1

default_elitism = True

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add optional settings
definition.add_optional("genome_type", "string", "type of genome", default_genome_type, choices=["list", "binary_string"])
definition.add_optional("genome_dimension", "positive_integer", "dimension of the genome", default_dimension, choices=genome_dimensions)
definition.add_optional("nparameters", "positive_integer", "number of parameters (genes) of the genome (in the first dimension)", default_nparameters)
definition.add_optional("nparameters2", "positive_integer", "number of parameters of the genome in the second dimension")
definition.add_optional("nindividuals", "integer", "number of individuals in one generation", default_nindividuals)
definition.add_optional("parameter_range", "range", "allowed range of the parameter")
definition.add_optional("best_raw_score", "real", "best raw score")
definition.add_optional("ngenerations", "positive_integer", "number of generations", default_ngenerations)
definition.add_optional("mutation_rate", "real", "mutation rate", default_mutation_rate)
definition.add_optional("crossover_rate", "real", "crossover rate", default_crossover_rate)
definition.add_optional("crossover_method", "string", "crossover method", default_crossover_method, choices=crossover_methods)
definition.add_optional("stats_freq", "integer", "frequency of printing statistics (in number of generations)", default_stats_freq)
definition.add_optional("round_decimal", "integer", "round everything to this decimal place", default_round_decimal)
definition.add_optional("mutation_method", "string", "mutation method", default_mutation_method, choices=mutation_methods)
definition.add_optional("min_or_max", "string", "minimize or maximize", default_min_or_max, choices=min_or_max)
definition.add_optional("run_id", "string", "identifier for this run", default="run0")
definition.add_optional("database_frequency", "positive_integer", "frequency of appending to the database (in the number of generations)", default_database_frequency)
definition.add_optional("statistics_frequency", "positive_integer", "frequency of appending to the statistics table", default_statistics_frequency)
definition.add_optional("scaling_method", "string", "scaling method", "linear", choices=scaling_methods)
definition.add_optional("selector_method", "string", "selector method", "rank", choices=selector_methods)

# Other
definition.add_optional("output", "directory_path", "output directory")

# Flags
definition.add_flag("elitism", "enable elitism", default_elitism)
definition.add_flag("show", "show results", True)
definition.add_flag("write", "write results", True)
definition.add_flag("plot", "plot results", False)
definition.add_flag("finish", "finish the evolution: set the scores of the last generation but don't generate a new population", False)
definition.add_flag("heterogeneous", "genomes use heterogeneous quantities (parameter minima and maxima - or centers and sigmas for gaussian initializers and mutators - must be specified as input to the 'run' function")
definition.add_flag("named_individuals", "use named individuals", False)

# Writing options
definition.add_section("writing", "writing options")
definition.sections["writing"].add_optional("engine_path", "string", "path for the genetic engine")
definition.sections["writing"].add_optional("prng_path", "string", "path for the prng")
definition.sections["writing"].add_optional("config_path", "string", "path for the configuration file")
definition.sections["writing"].add_optional("statistics_path", "string", "path for the statistics file")
definition.sections["writing"].add_optional("database_path", "string", "path for the database")
definition.sections["writing"].add_optional("population_path", "string", "path for the population data file")
definition.sections["writing"].add_optional("elitism_table_path", "string", "path for the elitism table")

# Advanced
definition.add_optional("nelite_individuals", "positive_integer", "number of individuals to take as elite", 1)

# Must only be defined when specifying only centers and sigmas for gaussian initializers and mutators, and no minima or maxima, or parameter_range
definition.add_optional("parameter_type", "string", "parameter type", choices=["real", "integer"])

# Recurrence of the same individuals
definition.add_flag("check_recurrence", "check for the recurrence of the same individuals", False)
definition.add_optional("recurrence_rtol", "positive_real", "relative tolerance for comparing equality of individuals for checking recurrence", 1e-5)
definition.add_optional("recurrence_atol", "positive_real", "absolute tolerance for comparing equality of individuals for checking recurrence", 1e-8)

# -----------------------------------------------------------------
