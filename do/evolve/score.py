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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.evolve.simplega import GSimpleGA, RawScoreCriteria
from pts.core.tools import serialization
from pts.core.tools import filesystem as fs
from pts.core.tools import tables

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
path = fs.join(fs.cwd(), "ga.pickle")

# Path to the parameters table
parameters_path = fs.join(fs.cwd(), "parameters.dat")

# -----------------------------------------------------------------

# Loda the GA
ga = GSimpleGA.from_file(path)

# Load the parameters table
table = tables.from_file(parameters_path,format="ascii.ecsv")

# Set scores of first population
index = 0
for ind in ga.internalPop:

    parameter_a = ind.genomeList[0]
    parameter_b = ind.genomeList[1]

    parameter_a_tab = table["Parameter a"][index]
    parameter_b_tab = table["Parameter b"][index]

    #rel_diff_a = abs((parameter_a - parameter_a_tab) / parameter_a)
    #rel_diff_b = abs((parameter_b - parameter_b_tab) / parameter_b)
    #print(rel_diff_a, rel_diff_b)

    assert np.isclose(parameter_a, parameter_a_tab, rtol=1e-11)
    assert np.isclose(parameter_b, parameter_b_tab, rtol=1e-11)

    #assert float(parameter_a) == float(parameter_a_tab), parameter_a - parameter_a_tab
    #assert float(parameter_b) == float(parameter_b_tab), parameter_b - parameter_b_tab

    # Set the score
    ind.score = chi_squared_function(ind)

    # Increment the index
    index += 1

# Dump
#ga.saveto(path)


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
