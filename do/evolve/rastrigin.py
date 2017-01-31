#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

from pts.evolve.engine import GAEngine, RawScoreCriteria
from pts.evolve.genomes.list1d import G1DList
from pts.evolve import mutators, initializators
from pts.evolve import selectors
from pts.evolve import constants
import math

# This is the Rastrigin Function, a deception function
def rastrigin(genome):
    
   n = len(genome)
   total = 0
   for i in xrange(n):
      total += genome[i]**2 - 10*math.cos(2*math.pi*genome[i])
      
   return (10*n) + total

def run_main():
    
   # Genome instance
   genome = G1DList(20)
   genome.setParams(rangemin=-5.2, rangemax=5.30, bestrawscore=0.00, rounddecimal=2)
   genome.initializator.set(initializators.G1DListInitializatorReal)
   genome.mutator.set(mutators.G1DListMutatorRealGaussian)

   genome.evaluator.set(rastrigin)

   # Genetic Algorithm Instance
   ga = GAEngine(genome)
   ga.terminationCriteria.set(RawScoreCriteria)
   ga.setMinimax(constants.minimaxType["minimize"])
   ga.setGenerations(3000)
   ga.setCrossoverRate(0.8)
   ga.setPopulationSize(100)
   ga.setMutationRate(0.06)

   ga.evolve(freq_stats=50)

   best = ga.bestIndividual()
   print best

# -----------------------------------------------------------------

if __name__ == "__main__":
   run_main()

# -----------------------------------------------------------------