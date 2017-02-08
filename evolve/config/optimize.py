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
definition.add_optional("nparameters", "integer", "number of parameters (genes) in the genome", 20)
definition.add_optional("nindividuals", "integer", "number of individuals in one generation", 80)
definition.add_optional("parameter_range", "range", "allowed range of the parameter")
definition.add_optional("best_raw_score", "real", "best raw score")
definition.add_optional("ngenerations", "positive_integer", "number of generations", 100)
definition.add_optional("mutation_rate", "real", "mutation rate", 0.02)
definition.add_optional("crossover_rate", "real", "crossover rate", 0.9)
definition.add_optional("stats_freq", "integer", "frequency of statistics (in number of generations)", 50)
definition.add_optional("best_raw_score", "real", "best score for an individual")
definition.add_optional("rounddecimal", "integer", "round everything to this decimal place")
definition.add_optional("mutation_method", "string", "mutation method", choices=["range", "gaussian", "binary"])
definition.add_optional("min_or_max", "string", "minimize or maximize", choices=["minimize", "maximize"])

# Flags
definition.add_flag("elitism", "enable elitism", True)
definition.add_flag("progress_bar", "use a progress bar to show the progress of the evoluation")
definition.add_flag("show", "show results", True)
definition.add_flag("write", "write results", True)

# Advanced
definition.add_optional("nelite_individuals", "positive_integer", "number of individuals to take as elite", 1)

# -----------------------------------------------------------------
