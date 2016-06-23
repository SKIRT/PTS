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

from pts.evolve.GSimpleGA import GSimpleGA, RawScoreCriteria
from pts.evolve import G1DList
from pts.evolve import Mutators, Initializators
from pts.evolve import Selectors
from pts.evolve import Consts
from scipy import stats

import matplotlib.pyplot as plt

import math

from pts.core.tools import serialization
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

config = Configuration()
config.add_required("dump_or_load", str, "dump or load")
config.read()

# -----------------------------------------------------------------

test_data_x = [20., 16., 19.79999924, 18.39999962, 17.10000038, 15.5, 14.69999981, 17.10000038, 15.39999962,
               16.20000076,
               15., 17.20000076, 16., 17., 14.39999962]
test_data_y = [88.59999847, 71.59999847, 93.30000305, 84.30000305, 80.59999847, 75.19999695, 69.69999695, 82.,
               69.40000153, 83.30000305, 79.59999847, 82.59999847, 80.59999847, 83.5, 76.30000305]

# -----------------------------------------------------------------

# plt.figure
# plt.scatter(test_data_x, test_data_y)
# plt.show()

def chi_squared_function(chromosome):

    chi_squared = 0.0

    for i in range(len(test_data_x)):

        x = test_data_x[i]
        y = test_data_y[i]

        chromosome_y = chromosome[0] * x + chromosome[1]

        chi_squared += (y - chromosome_y) ** 2.

    chi_squared /= 2.0

    return chi_squared

# path to the GA object
path = fs.join(fs.home(), "ga.pickle")


if config.arguments.dump_or_load == "dump":

    # Genome instance
    genome = G1DList.G1DList(2)
    genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
    genome.initializator.set(Initializators.G1DListInitializatorReal)
    genome.mutator.set(Mutators.G1DListMutatorRealGaussian)

    genome.evaluator.set(chi_squared_function)

    # Genetic Algorithm Instance
    ga = GSimpleGA(genome)
    ga.terminationCriteria.set(RawScoreCriteria)
    ga.setMinimax(Consts.minimaxType["minimize"])
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



    # Dump
    ga.saveto(path)

elif config.arguments.dump_or_load == "load":

    ga = GSimpleGA.from_file(path)

    pop = ga.internalPop
    for ind in pop: print(ind.genomeList)

# newpop = ga.generate_new_population()
# for ind in newpop: print(ind.genomeList)

#pop = ga.getPopulation()
#print(pop)
#best = ga.bestIndividual()

#slope, intercept, r_value, p_value, std_err = stats.linregress(test_data_x, test_data_y)
#print("slope", best[0], "(real:", slope, ")")
#print("intercept", best[1], "(real:", intercept, ")")

# -----------------------------------------------------------------
