#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization Contains the InputInitializer

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import bisect
import numpy as np
#import pyevolve

from pyevolve import G1DList
from pyevolve import GSimpleGA

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import serialization

# -----------------------------------------------------------------

# INSPYPRED:

# https://pypi.python.org/pypi/inspyred
# https://aarongarrett.github.io/inspyred/
#
# EXAMPLE:
#
# import random
# import time
# import inspyred
#
# def generate_binary(random, args):
#    bits = args.get('num_bits', 8)
#    return [random.choice([0, 1]) for i in range(bits)]
#
# @inspyred.ec.evaluators.evaluator
# def evaluate_binary(candidate, args):
#    return int("".join([str(c) for c in candidate]), 2)
#
#rand = random.Random()
#rand.seed(int(time.time()))
#ga = inspyred.ec.GA(rand)
#ga.observer = inspyred.ec.observers.stats_observer
#ga.terminator = inspyred.ec.terminators.evaluation_termination
#final_pop = ga.evolve(evaluator=evaluate_binary,
#                      generator=generate_binary,
#                      max_evaluations=1000,
#                      num_elites=1,
#                      pop_size=100,
#                      num_bits=10)
#final_pop.sort(reverse=True)
#for ind in final_pop:
#    print(str(ind))
#

#

# -----------------------------------------------------------------

# GPLEARN:

# http://gplearn.readthedocs.io/en/latest/index.html
# NO, THIS IS FOR MATHEMATICAL SYMBOLIC REGRESSION

# -----------------------------------------------------------------

# PYEVOLVE:

# http://pyevolve.sourceforge.net/graphs.html#graphs-screens

# https://github.com/perone/Pyevolve

# http://pyevolve.sourceforge.net/0_6rc1/  # LATEST RC (what is installed through pip)

## Create genome:

# Genome instance
#genome = G1DList.G1DList(20)

# The evaluator function (objective function)
#genome.evaluator.set(eval_func)

# This will create an instance of the G1DList.G1DList class (which resides in the G1DList module) with the
# list n-size of 20 and sets the evaluation function of the genome to the evaluation function “eval_func”
# that we have created before.

## Set ranges to the genome:

# genome.setParams(rangemin=0, rangemax=10)

## Create the GA Engine

# ga = GSimpleGA.GSimpleGA(genome)

# By default, the GA will evolve for 100 generations with a population size of 80 individuals,
# it will use the mutation rate of 2% and a crossover rate of 80%, the default selector is the Ranking Selection
# (Selectors.GRankSelector()) method. Those default parameters was not random picked, they are all based on the commom
# used properties.

## Evolve:

# Do the evolution, with stats dump
# frequency of 10 generations
#ga.evolve(freq_stats=10)

# Best individual
#print ga.bestIndividual()

###


# It is important to note that in Pyevolve, we have raw score and fitness score,
# the raw score is the return of the evaluation function and the fitness score is the scaled score or the raw
# score in absence of a scaling scheme.


# In our case genome is a list of 3 values (FUV young, FUV ionizing, dust mass)

genome = G1DList.G1DList(3)

ga = GSimpleGA.GSimpleGA(genome)


# Pyevolve have introduced the concept of the Interactive Mode in the course of evolution. When you are evolving,
# and the Interactive Mode is enabled, you can press the ESC Key anytime in the evolution process.
# By pressing that key, you will enter in the interactive mode, with a normal python prompt and the Interaction
# module exposed to you as the “it” module.


#ga_engine: The GSimpleGA.GSimpleGA instance, the GA Engine.
#it: The Interaction module, with the utilities and graph plotting functions.
#population: The current population.
#pyevolve: The main namespace, the pyevolve module.


# Plotting the current population raw scores histogram
# it.plotHistPopScore(population)

# Plotting the current population raw scores distribution
# it.plotPopScore(population)

# Get all the population raw scores
# popScores = it.getPopScores(population)
# popScores
# [17.0, 17.0, 16.0, 15.0, 13.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 9.0,
# 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
# 8.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.0, 5.0,
# 5.0, 5.0, 5.0, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0]

# This all comes from http://pyevolve.sourceforge.net/0_6rc1/getstarted.html#first-example


# EXAMPLE CHROMOSOME EVALUATION FUNCTION
# IN OUR CASE, THIS IS COMPLEX !! CREATING SKI FILE, SUBMITTING SKIRT SIMULATION TO CLUSTER, RETURNING SIMULATION OUTPUT,
# CALCULATING FLUXES FOR DIFFERENT BANDS, CHI SQUARED COMPARED TO OBSERVED FLUXES
#

def eval_func(chromosome):
   score = 0.0
   # iterate over the chromosome elements (items)
   for value in chromosome:
      if value==0:
         score += 1.0

   return score


# WHAT HAPPENS DURING C0NSTRUCTOR OF GSIMPLEGA:

# self.internalPop  = GPopulation(genome) # internal population is constructed from genome

    # THE CONSTRUCTOR OF GPOPULATION:

    # self.oneSelfGenome  = genome
    # self.internalPop    = []
    # self.internalPopRaw = []
    # self.popSize       = 0
    # ... among other things

# further things during constructor of GSIMPLEGA, setting properties:

#self.nGenerations = Consts.CDefGAGenerations
#self.pMutation = Consts.CDefGAMutationRate
#self.pCrossover = Consts.CDefGACrossoverRate
#self.nElitismReplacement = Consts.CDefGAElitismReplacement
#self.setPopulationSize(Consts.CDefGAPopulationSize)
#self.minimax = Consts.minimaxType["maximize"]
#self.elitism = True

# ==> a population consisting of one genome is constructed (=self.internalPop)

# WHAT HAPPENS WHEN EXECUTING .EVOLVE():

# while True:
#   self.step(): Just do one step in evolution, one generation

# WHAT HAPPENS WHEN EXECUTING .STEP():

# newPop = GPopulation(self.internalPop) # clone population
# for i in xrange(0, size_iterate, 2):

#     genomeMom = self.select(popID=self.currentGeneration)  # selects one individual from the population
#     genomeDad = self.select(popID=self.currentGeneration)

## Create sister and brother from mom and dad by applying cross-over

#     if not crossover_empty and self.pCrossover >= 1.0:
#     for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2):
#       (sister, brother) = it

## Mutate sister and brother genomes

#     sister.mutate(pmut=self.pMutation, ga_engine=self)
#     brother.mutate(pmut=self.pMutation, ga_engine=self)

#####




test_data_x = [20., 16., 19.79999924, 18.39999962, 17.10000038, 15.5, 14.69999981, 17.10000038, 15.39999962, 16.20000076,
               15., 17.20000076, 16., 17., 14.39999962]
test_data_y = [88.59999847, 71.59999847, 93.30000305, 84.30000305, 80.59999847, 75.19999695, 69.69999695, 82.,
               69.40000153, 83.30000305, 79.59999847, 82.59999847, 80.59999847, 83.5, 76.30000305]


def run_test():

    """
    This function ...
    :return:
    """

    genome = G1DList.G1DList(20)

    genome.evaluator.set(eval_func)

    ga = GSimpleGA.GSimpleGA(genome)

    ga.nGenerations = 10

    ga.evolve(freq_stats=1)

    #print(ga.bestIndividual())

# -----------------------------------------------------------------

class GeneticAlgorithm(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @classmethod
    def from_file(self, path):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_new_population(self):

        """
        This function ...
        :return:
        """

        # ga.step() : Just do one step in evolution, one generation

        # ga.evolve(self, freq_stats=0):
        # Do all the generations until the termination criteria, accepts
        # the freq_stats (default is 0) to dump statistics at n-generation


# -----------------------------------------------------------------
