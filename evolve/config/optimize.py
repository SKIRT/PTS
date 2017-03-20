#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add optional settings
definition.add_optional("genome_type", "string", "type of genome", default="list", choices=["list", "binary_string"])
definition.add_optional("genome_dimension", "positive_integer", "dimension of the genome", default=1, choices=[1, 2])
definition.add_optional("nparameters", "positive_integer", "number of parameters (genes) of the genome (in the first dimension)", 20)
definition.add_optional("nparameters2", "positive_integer", "number of parameters of the genome in the second dimension")
definition.add_optional("nindividuals", "integer", "number of individuals in one generation", 80)
definition.add_optional("parameter_range", "range", "allowed range of the parameter")
definition.add_optional("best_raw_score", "real", "best raw score")
definition.add_optional("ngenerations", "positive_integer", "number of generations", 100)
definition.add_optional("mutation_rate", "real", "mutation rate", 0.02)
definition.add_optional("crossover_rate", "real", "crossover rate", 0.9)
definition.add_optional("stats_freq", "integer", "frequency of statistics (in number of generations)", 50)
definition.add_optional("best_raw_score", "real", "best score for an individual")
definition.add_optional("round_decimal", "integer", "round everything to this decimal place", 2)
definition.add_optional("mutation_method", "string", "mutation method", choices=["range", "gaussian", "binary"])
definition.add_optional("min_or_max", "string", "minimize or maximize", choices=["minimize", "maximize"])
definition.add_optional("run_id", "string", "identifier for this run", default="run0")
definition.add_optional("database_frequency", "positive_integer", "frequency of appending to the database (in the number of generations)", 1)
definition.add_optional("statistics_frequency", "positive_integer", "frequency of appending to the statistics table", 1)

# Other
definition.add_optional("output", "directory_path", "output directory")

# Flags
definition.add_flag("elitism", "enable elitism", True)
definition.add_flag("show", "show results", True)
definition.add_flag("write", "write results", True)
definition.add_flag("plot", "plot results", False)
definition.add_flag("finish", "finish the evolution: set the scores of the last generation but don't generate a new population", False)

definition.add_flag("heterogeneous", "genomes use heterogeneous quantities (parameter minima and maxima - or centers and sigmas for gaussian initializers and mutators - must be specified as input to the 'run' function")

# Writing options
definition.add_section("writing", "writing options")
definition.sections["writing"].add_optional("engine_path", "string", "path for the genetic engine")
definition.sections["writing"].add_optional("prng_path", "string", "path for the prng")
definition.sections["writing"].add_optional("config_path", "string", "path for the configuration file")
definition.sections["writing"].add_optional("statistics_path", "string", "path for the statistics file")
definition.sections["writing"].add_optional("database_path", "string", "path for the database")

# Advanced
definition.add_optional("nelite_individuals", "positive_integer", "number of individuals to take as elite", 1)

# Must only be defined when specifying only centers and sigmas for gaussian initializers and mutators, and no minima or maxima, or parameter_range
definition.add_optional("parameter_type", "string", "parameter type", choices=["real", "integer"])

# -----------------------------------------------------------------
