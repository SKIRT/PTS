#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.example Example of using the 'evolve' subpackage for genetic algorithms.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import tables
from pts.core.tools import filesystem as fs
from pts.evolve.simplega import GSimpleGA, RawScoreCriteria
from pts.evolve.genomes.list1d import G1DList
from pts.evolve import mutators
from pts.evolve import initializators
from pts.evolve import constants

# -----------------------------------------------------------------

# path to the GA object
path = fs.join(fs.cwd(), "ga.pickle")

# Path to the parameter table
parameters_path = fs.join(fs.cwd(), "parameters.dat")

# -----------------------------------------------------------------

# Genome instance
genome = G1DList(2)
genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
genome.initializator.set(initializators.G1DListInitializatorReal)
genome.mutator.set(mutators.G1DListMutatorRealGaussian)

#genome.evaluator.set(chi_squared_function)

# Genetic Algorithm Instance
ga = GSimpleGA(genome)
ga.terminationCriteria.set(RawScoreCriteria)
ga.setMinimax(constants.minimaxType["minimize"])
ga.setGenerations(20)
ga.setCrossoverRate(0.5)
ga.setPopulationSize(100)
ga.setMutationRate(0.5)

# ga.initialize_evolution()

ga.initialize()

par_a_column = []
par_b_column = []
#par_c_column = []

pop = ga.internalPop
for ind in pop:

    #print(ind.genomeList)

    par_a_column.append(ind.genomeList[0])
    par_b_column.append(ind.genomeList[1])

#data = [par_a_column, par_b_column, par_c_column]
#names = ["Parameter a", "Parameter b", "Parameter c"]
data = [par_a_column, par_b_column]
names = ["Parameter a", "Parameter b"]

parameters_table = tables.new(data, names)

# Save the genetic algorithm
ga.saveto(path)

# Save the parameter table
tables.write(parameters_table, parameters_path, format="ascii.ecsv")

# -----------------------------------------------------------------
