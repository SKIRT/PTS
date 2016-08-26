#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.original Reference, original implementation.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Import the relevant PTS classes and modules
from pts.evolve.engine import GAEngine, RawScoreCriteria
from pts.evolve.genomes.list1d import G1DList
from pts.evolve import mutators
from pts.evolve import initializators
from pts.evolve import constants
from pts.core.tools import filesystem as fs
from pts.core.tools.random import setup_prng
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("seed", int, "the random seed", 4357)

# Get configuration
setter = ArgumentConfigurationSetter("reference")
config = setter.run(definition)

# -----------------------------------------------------------------

x = np.linspace(12,25,100)

test_data_x = [20., 16., 19.79999924, 18.39999962, 17.10000038, 15.5, 14.69999981, 17.10000038, 15.39999962,
               16.20000076,
               15., 17.20000076, 16., 17., 14.39999962]
test_data_y = [88.59999847, 71.59999847, 93.30000305, 84.30000305, 80.59999847, 75.19999695, 69.69999695, 82.,
               69.40000153, 83.30000305, 79.59999847, 82.59999847, 80.59999847, 83.5, 76.30000305]

# -----------------------------------------------------------------

def fit_function(x, a, b):

    """
    This function ...
    :param a:
    :param b:
    :param x:
    :return:
    """

    return a * x + b

# -----------------------------------------------------------------

def chi_squared_function(chromosome):

   """
   This function calculates the chi-squared value for a certain set of parameters (chromosome)
   :param chromosome:
   :return:
   """

   chi_squared = 0.0
   for i in range(len(test_data_x)):
      x = test_data_x[i]
      y = test_data_y[i]
      chromosome_y = fit_function(x, chromosome[0], chromosome[1])
      chi_squared += (y - chromosome_y) ** 2.
   chi_squared /= 2.0
   return chi_squared

# -----------------------------------------------------------------

#seed = 4357
seed = config.seed
prng = setup_prng(seed)

# -----------------------------------------------------------------

# Genome instance
genome = G1DList(2)
genome.setParams(rangemin=0., rangemax=50., bestrawscore=0.00, rounddecimal=2)
genome.initializator.set(initializators.G1DListInitializatorReal)
genome.mutator.set(mutators.G1DListMutatorRealGaussian)

# Set the evaluator function
genome.evaluator.set(chi_squared_function)

# Genetic algorithm instance
ga = GAEngine(genome)
ga.terminationCriteria.set(RawScoreCriteria)
ga.setMinimax(constants.minimaxType["minimize"])
ga.setGenerations(5)
ga.setCrossoverRate(0.5)
ga.setPopulationSize(100)
ga.setMutationRate(0.5)

# Evolve
#ga.evolve(freq_stats=False)
ga.evolve()

print("Final generation:", ga.currentGeneration)

# -----------------------------------------------------------------

# Determine the path to the reference directory
ref_path = fs.join(fs.cwd(), "reference")
fs.create_directory(ref_path)

# -----------------------------------------------------------------

best = ga.bestIndividual()

best_parameter_a = best.genomeList[0]
best_parameter_b = best.genomeList[1]

best_path = fs.join(ref_path, "best.dat")

with open(best_path, 'w') as best_file:
    best_file.write("Parameter a: " + str(best_parameter_a) + "\n")
    best_file.write("Parameter b: " + str(best_parameter_b) + "\n")

popt, pcov = curve_fit(fit_function, test_data_x, test_data_y)
parameter_a_real = popt[0]
parameter_b_real = popt[1]

print("Best parameter a:", best_parameter_a, " REAL:", parameter_a_real)
print("Best parameter b:", best_parameter_b, " REAL:", parameter_b_real)

plt.figure()
plt.scatter(test_data_x, test_data_y)
plt.plot(x, [fit_function(x_i, best_parameter_a, best_parameter_b) for x_i in x])
plt.plot(x, [fit_function(x_i, parameter_a_real, parameter_b_real) for x_i in x])
plt.ylim(65, 95)
plt.xlim(12,22)

# Save the figure
plot_path = fs.join(ref_path, "best.pdf")
plt.savefig(plot_path)

# -----------------------------------------------------------------
