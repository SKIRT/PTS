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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.evolve.core.dbadapters import DBFileCSV
from pts.core.tools.logging import setup_log

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("ngenerations", "positive_integer", "number of generations (steps)", 10)

# Create the configuration
setter = ArgumentConfigurationSetter("test_tsp_steps")
config = setter.run(definition)

# Set logging
log = setup_log("DEBUG")

# -----------------------------------------------------------------

evaluator_kwargs = input_optimize["evaluator_kwargs"]

# -----------------------------------------------------------------

# Determine temporary directory
path = fs.create_directory_in(introspection.pts_temp_dir, time.unique_name("optimize"))

# Change working directory
fs.change_cwd(path)

# Determine the path to the database
#database_path = fs.join(path, "database.csv")
#database_path = fs.join(path, "database.db")

# Create database adapter
#adapter = DBFileCSV(filename=database_path, identify="run_01", frequency=1, reset=True)
#adapter = DBSQLite(filename=database_path, identify="test", resetDB=True)

# Set the adapter
#input_optimize["adapter"] = adapter

# -----------------------------------------------------------------

# Create the step optimizer
optimizer = StepOptimizer(settings_optimize)

# Run the optimizer
population = optimizer.run(**input_optimize)

# Save the optimizer
optimizer.saveto(path)

# Remove the variable, to load it again later
del optimizer

# -----------------------------------------------------------------

# Repeat for a certain number of times
for _ in range(config.ngenerations):

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

    #print(repr(optimizer.engine))

    # Set new input
    new_input = {"scores": scores}

    # Run again
    population = optimizer.run(**new_input)

    # Save the optimizer
    optimizer.saveto(path)

    # Remove the variable, to load it again later
    del optimizer

# -----------------------------------------------------------------

#print(adapter)

# Close the database
#adapter.close()

# -----------------------------------------------------------------
