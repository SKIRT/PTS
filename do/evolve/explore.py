#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.explore Do first model exploration for GA test.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import tables
from pts.core.tools import filesystem as fs
from pts.evolve.engine import GAEngine, RawScoreCriteria
from pts.evolve.genomes.list1d import G1DList
from pts.evolve import mutators
from pts.evolve import initializators
from pts.evolve import constants
from pts.core.tools.logging import log
from pts.core.tools import time
from pts.core.tools.random import setup_prng, save_state
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Configuration
config = Configuration()
config.add_positional_optional("seed", int, "the random seed", 4357)
config.read()

# -----------------------------------------------------------------

seed = config.arguments.seed
prng = setup_prng(seed)

# -----------------------------------------------------------------

# path to the GA object
path = fs.join(fs.cwd(), "ga.pickle")

# Path to the parameter table
parameters_path = fs.join(fs.cwd(), "parameters.dat")

# -----------------------------------------------------------------

# Inform the user
log.info("Creating the GA engine ...")

# Genome instance
genome = G1DList(2)
genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
genome.initializator.set(initializators.G1DListInitializatorReal)
genome.mutator.set(mutators.G1DListMutatorRealGaussian)

#genome.evaluator.set(chi_squared_function)

# Inform the user
log.info("Creating the GA engine ...")

# Genetic algorithm instance
ga = GAEngine(genome)
ga.terminationCriteria.set(RawScoreCriteria)
ga.setMinimax(constants.minimaxType["minimize"])
ga.setGenerations(5)
ga.setCrossoverRate(0.5)
ga.setPopulationSize(100)
ga.setMutationRate(0.5)

# Initialize the genetic algorithm
ga.initialize()



name_column = []
par_a_column = []
par_b_column = []

pop = ga.internalPop
for ind in pop:

    # Give the individual a unique name
    name = time.unique_name(precision="micro")
    name_column.append(name)
    par_a_column.append(ind.genomeList[0])
    par_b_column.append(ind.genomeList[1])

# Create the parameters table
data = [name_column, par_a_column, par_b_column]
names = ["Unique name", "Parameter a", "Parameter b"]
parameters_table = tables.new(data, names)

# Save the genetic algorithm
ga.saveto(path)

#print("Current generation: ", ga.currentGeneration)

# Save the parameter table
tables.write(parameters_table, parameters_path, format="ascii.ecsv")

# -----------------------------------------------------------------

# Path to the random state
random_path = fs.join(fs.cwd(), "rndstate.pickle")

# Save the state of the random generator
save_state(random_path)

# -----------------------------------------------------------------
