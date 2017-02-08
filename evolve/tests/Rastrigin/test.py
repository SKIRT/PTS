#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import inspect

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.core.basics.range import RealRange
from pts.evolve import reference
from pts.evolve.optimize import show_best

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "optimizing the Rastrigin function, a deceptive function"

# -----------------------------------------------------------------

# Define properties
nparameters = 20
nindividuals = 80
parameter_range = RealRange(-5.2, 5.30)
best_raw_score = 0.0
round_decimal = None
ngenerations = 800
mutation_rate = 0.05
crossover_rate = None
stats_freq = 50
mutation_method = "gaussian" # or range, or binary
min_or_max = "minimize"

# -----------------------------------------------------------------

def rastringin(genome):

    """
    This is the Rastringin function, a deceptive function
    :param genome:
    :return:
    """

    n = len(genome)
    total = 0
    for i in xrange(n):
        total += genome[i] ** 2 - 10 * math.cos(2 * math.pi * genome[i])

    return (10 * n) + total

# -----------------------------------------------------------------

# Initialize list for the commands
commands = []

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
settings_optimize = dict()
settings_optimize["output"] = None
settings_optimize["nparameters"] = nparameters
settings_optimize["nindividuals"] = nindividuals
settings_optimize["parameter_range"] = parameter_range
settings_optimize["best_raw_score"] = best_raw_score
settings_optimize["round_decimal"] = round_decimal
settings_optimize["ngenerations"] = ngenerations
settings_optimize["mutation_rate"] = mutation_rate
settings_optimize["crossover_rate"] = crossover_rate
settings_optimize["stats_freq"] = stats_freq
settings_optimize["mutation_method"] = mutation_method
settings_optimize["min_or_max"] = min_or_max

# Other
settings_optimize["progress_bar"] = True

# Input
input_optimize = dict()
input_optimize["evaluator"] = rastringin

# Construct the command
optimize = Command("optimize", "optimize the Rastrigin function", settings_optimize, input_optimize, cwd=".")

# Add the command
commands.append(optimize)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    # Solve the problem with the original Pyevolve implementation
    best = reference.call(settings_optimize, rastringin)

    # Show the best individual
    show_best(best)

# -----------------------------------------------------------------
