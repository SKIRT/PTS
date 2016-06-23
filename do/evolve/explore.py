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

import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from pts.evolve.simplega import GSimpleGA, RawScoreCriteria
from pts.evolve.genomes.list1d import G1DList
from pts.evolve import mutators
from pts.evolve import initializators
from pts.evolve import selectors
from pts.evolve import constants
from scipy import stats
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

test_data_x = [20., 16., 19.79999924, 18.39999962, 17.10000038, 15.5, 14.69999981, 17.10000038, 15.39999962,
               16.20000076,
               15., 17.20000076, 16., 17., 14.39999962]
test_data_y = [88.59999847, 71.59999847, 93.30000305, 84.30000305, 80.59999847, 75.19999695, 69.69999695, 82.,
               69.40000153, 83.30000305, 79.59999847, 82.59999847, 80.59999847, 83.5, 76.30000305]

# -----------------------------------------------------------------

def chi_squared_function(chromosome):
    chi_squared = 0.0
    for i in range(len(test_data_x)):
        x = test_data_x[i]
        y = test_data_y[i]
        chromosome_y = chromosome[0] * x + chromosome[1]
        chi_squared += (y - chromosome_y) ** 2.
    chi_squared /= 2.0
    return chi_squared

# -----------------------------------------------------------------

# path to the GA object
path = fs.join(fs.home(), "ga.pickle")

# Genome instance
genome = G1DList(2)
genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
genome.initializator.set(initializators.G1DListInitializatorReal)
genome.mutator.set(mutators.G1DListMutatorRealGaussian)

genome.evaluator.set(chi_squared_function)

# Genetic Algorithm Instance
ga = GSimpleGA(genome)
ga.terminationCriteria.set(RawScoreCriteria)
ga.setMinimax(constants.minimaxType["minimize"])
ga.setGenerations(20)
ga.setCrossoverRate(0.5)
ga.setPopulationSize(100)
ga.setMutationRate(0.5)

#ga.evolve(freq_stats=50)
#exit()

# ga.initialize_evolution()

ga.initialize()

#pop = ga.internalPop
#for ind in pop: print(ind.genomeList)

# Set scores of first population
#for ind in self.internalPop:

# Dump
ga.saveto(path)

# -----------------------------------------------------------------
