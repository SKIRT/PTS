#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.evolve Evolve the genetic algorithm: produce the next generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools import time
from pts.evolve.optimize.stepwise import StepWiseOptimizer
from pts.evolve.optimize.optimizer import show_best
from pts.evolve.tests.TravelingSalesman.test import settings_optimize, input_optimize, eval_func
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import setup_log

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("ngenerations", "positive_integer", "number of generations (steps)", 10)

# Create the configuration
config = parse_arguments("test_tsp_steps", definition)

# Set logging
log = setup_log("DEBUG")

# -----------------------------------------------------------------

evaluator_kwargs = input_optimize["evaluator_kwargs"]

# -----------------------------------------------------------------

# Determine temporary directory
path = introspection.create_temp_dir(time.unique_name("optimize"))

# Change working directory
fs.change_cwd(path)

# -----------------------------------------------------------------

# Create the step optimizer
optimizer = StepWiseOptimizer(settings_optimize)

# Run the optimizer
population = optimizer.run(**input_optimize)

# Remove the variable, to load it again later
del optimizer

# -----------------------------------------------------------------

best = None

# -----------------------------------------------------------------

# Repeat for a certain number of times
for index in range(config.ngenerations):

    # Inform the user
    log.info("Evaluating and setting scores for the new population ...")

    # Initialize the list of scores
    scores = []

    # Determine the scores
    for genome in population:

        # Determine the score
        score = eval_func(genome, **evaluator_kwargs)

        # Add the score
        scores.append(score)

    # -----------------------------------------------------------------

    # Load the optimizer
    optimizer = StepWiseOptimizer.from_directory(path)

    # Set new input
    new_input = {"scores": scores}

    # -----------------------------------------------------------------

    # Last generation
    if index == config.ngenerations - 1:

        # Inform the user
        log.info("Finishing evolution ...")

        # Set finish flag
        optimizer.config.finish = True

    # Continue with next step
    else: log.info("Advancing evolution ...")

    # Run again
    population = optimizer.run(**new_input)

    # Get the best individual
    best = optimizer.best

    # Remove the variable, to load it again later
    del optimizer

# -----------------------------------------------------------------

# Show the best
show_best(best)

# -----------------------------------------------------------------
