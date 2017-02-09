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
from pts.evolve.optimize.step import StepOptimizer
from pts.evolve.tests.TravelingSalesman.test import settings_optimize, input_optimize, eval_func

# -----------------------------------------------------------------

evaluator_kwargs = input_optimize["evaluator_kwargs"]

# -----------------------------------------------------------------

# Determine temporary directory
path = fs.create_directory_in(introspection.pts_temp_dir, time.unique_name("optimize"))

# Change working directory
fs.change_cwd(path)

# Create the step optimizer
optimizer = StepOptimizer(settings_optimize)

# Run the optimizer
population = optimizer.run(**input_optimize)

# Save the optimizer
optimizer.saveto(path)

# Remove the variable, to load it again later
del optimizer

# -----------------------------------------------------------------

# Initialize the list of scores
scores = []

# Determine the scores
for genome in population:

    # Determine the score
    score = eval_func(genome, **evaluator_kwargs)

    # Add the score
    scores.append(score)

# -----------------------------------------------------------------

# Load the optimizer again
optimizer = StepOptimizer.from_file(path)

# Set new input
new_input = {"scores": scores}

# Run again
population = optimizer.run(**new_input)

# Save the optimizer
optimizer.saveto(path)

# Remove the variable, to load it again later
del optimizer

# -----------------------------------------------------------------



print(path)

# -----------------------------------------------------------------
