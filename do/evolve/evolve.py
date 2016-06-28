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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.evolve.simplega import SimpleGeneticAlgorithm
from pts.core.tools import filesystem as fs
from pts.core.tools import tables
from pts.core.tools import time
from pts.core.tools.random import load_state

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
    print("current generation: the initial population")
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
ga = SimpleGeneticAlgorithm.from_file(path)

print("Current generation according to GA: ", ga.currentGeneration)

pop_for_loop = ga.new_population if ga.new_population is not None else ga.internalPop

# Set scores of first population
index = 0
for ind in pop_for_loop:

    parameter_a = ind.genomeList[0]
    parameter_b = ind.genomeList[1]

    parameter_a_tab = parameters["Parameter a"][index]
    parameter_b_tab = parameters["Parameter b"][index]

    rel_diff_a = abs((parameter_a - parameter_a_tab) / parameter_a)
    rel_diff_b = abs((parameter_b - parameter_b_tab) / parameter_b)
    #print(rel_diff_a, rel_diff_b)

    #print(parameter_a, parameter_a_tab)
    #print(parameter_b, parameter_b_tab)

    assert np.isclose(parameter_a, parameter_a_tab, rtol=1e-11), rel_diff_a
    assert np.isclose(parameter_b, parameter_b_tab, rtol=1e-11), rel_diff_b

    assert parameters["Unique name"][index] == chi_squared["Unique name"][index]

    # Get the score
    score = chi_squared["Chi-squared"][index]

    # Set the score
    ind.score = score

    # Increment the index
    index += 1

# Replace
if ga.new_population is not None: ga.replace_internal_population()

# Sort the internal population
ga.internalPop.sort()

# Set new pop to None
ga.new_population = None

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

# Initialize evolution:
# self.initialize() # done in 'explore'
# self.internalPop.evaluate()
# self.internalPop.sort()


# Evaluate function of population is just loop over evaluate function of individuals
# Evaluate function of individual (genome) is just calculating the sum of scores of each target function


# newpop = ga.generate_new_population()
# for ind in newpop: print(ind.genomeList)

#pop = ga.getPopulation()
#print(pop)
#best = ga.bestIndividual()

#slope, intercept, r_value, p_value, std_err = stats.linregress(test_data_x, test_data_y)
#print("slope", best[0], "(real:", slope, ")")
#print("intercept", best[1], "(real:", intercept, ")")

# -----------------------------------------------------------------
