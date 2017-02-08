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
from pts.evolve.engine import GeneticEngine
from pts.core.tools import filesystem as fs
from pts.core.tools import tables
from pts.core.tools import time
from pts.core.tools.random import load_state, save_state

# -----------------------------------------------------------------

last_generation = None

# Check the index of the last generation
for name in fs.directories_in_path():
    if "reference" in name or "original" in name: continue
    generation = int(name.split("Generation ")[1])
    if last_generation is None or generation > last_generation: last_generation = generation

# -----------------------------------------------------------------

if last_generation is None:
    generation_path = fs.cwd()
    print("Current generation: the initial population")
else:
    generation_path = fs.join(fs.cwd(), "Generation " + str(last_generation))
    print("Current generation: " + str(last_generation))

# -----------------------------------------------------------------

# Path to the current GA object
path = fs.join(generation_path, "ga.pickle")

# Path to the parameters table
parameters_path = fs.join(generation_path, "parameters.dat")

# Path to the chi squared table
chi_squared_path = fs.join(generation_path, "chi_squared.dat")

# Check whether the generation is scored
if not fs.is_file(chi_squared_path): raise RuntimeError("The last generation has not been scored yet!")

# Path to the random state
random_path = fs.join(generation_path, "rndstate.pickle")

# -----------------------------------------------------------------

# Load the random state
load_state(random_path)

# -----------------------------------------------------------------

# Load the parameters table
parameters = tables.from_file(parameters_path, format="ascii.ecsv")

# Load the chi squared table
chi_squared = tables.from_file(chi_squared_path, format="ascii.ecsv")

# -----------------------------------------------------------------

# Load the GA
ga = GeneticEngine.from_file(path)

# Check whether the chi-squared and parameter tables match
for i in range(len(parameters)):
    assert parameters["Unique name"][i] == chi_squared["Unique name"][i]

# Get the scores
scores = chi_squared["Chi-squared"]
check = parameters

# Set the scores
ga.set_scores(scores, check)

# -----------------------------------------------------------------

new_generation = last_generation + 1 if last_generation is not None else 0

# Path to the new generation
new_generation_path = fs.join(fs.cwd(), "Generation " + str(new_generation))
fs.create_directory(new_generation_path)

# path to the new GA instance
new_path = fs.join(new_generation_path, "ga.pickle")

# path to the new parameters table
new_parameters_path = fs.join(new_generation_path, "parameters.dat")

# -----------------------------------------------------------------

# Generate the new population
ga.generate_new_population()

# -----------------------------------------------------------------

name_column = []
par_a_column = []
par_b_column = []

for ind in ga.new_population:

    # Give the individual a unique name
    name = time.unique_name(precision="micro")
    name_column.append(name)
    par_a_column.append(ind.genomeList[0])
    par_b_column.append(ind.genomeList[1])

# Create the parameters table
data = [name_column, par_a_column, par_b_column]
names = ["Unique name", "Parameter a", "Parameter b"]
new_parameters_table = tables.new(data, names)

# -----------------------------------------------------------------

# Save the new parameters table
tables.write(new_parameters_table, new_parameters_path, format="ascii.ecsv")

# Dump the GA
ga.saveto(new_path)

# -----------------------------------------------------------------

# Save the state of the random generator
new_random_path = fs.join(new_generation_path, "rndstate.pickle")
save_state(new_random_path)

# -----------------------------------------------------------------

# Initialize evolution:
# self.initialize() # done in 'explore'
# self.internalPop.evaluate()
# self.internalPop.sort()


# Evaluate function of population is just loop over evaluate function of individuals
# Evaluate function of individual (genome) is just calculating the sum of scores of each target function


# newpop = ga.generate_new_population()
# for ind in newpop: print(ind.genomeList)

#pop = ga.get_population()
#print(pop)
#best = ga.bestIndividual()

#slope, intercept, r_value, p_value, std_err = stats.linregress(test_data_x, test_data_y)
#print("slope", best[0], "(real:", slope, ")")
#print("intercept", best[1], "(real:", intercept, ")")

# -----------------------------------------------------------------
