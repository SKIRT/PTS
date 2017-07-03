#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.interaction In this module, you will find the funcionality for the :term:`Interactive mode`.
#  When you enter in the Interactive Mode, Pyevolve will automatic import this module
#  and exposes to you in the name space called "it".
#
#To use this mode, the parameter *interactiveMode* must be enabled in the
#:class:`GSimpleGA.GSimpleGA`.
#
#You can use the manual method to enter in the Interactive Mode at specific
#generation using the :meth:`GSimpleGA.GSimpleGA.setInteractiveGeneration` method.
#

# -----------------------------------------------------------------

# Import standard modules
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy

# -----------------------------------------------------------------

def get_population_scores(population, fitness=False):

   """
   Returns a list of population scores
   Example:
      >>> lst = get_population_scores(population)

   :param population: population object (:class:`GPopulation.GPopulation`)
   :param fitness: if True, the fitness score will be used, otherwise, the raw.
   :rtype: list of population scores
   """

   score_list = []
   for individual in population:
      score_list.append(individual.fitness if fitness else individual.score)
   return score_list

# -----------------------------------------------------------------

def plot_population_score(population, fitness=False):

   """
   Plot the population score distribution
   Example:
      >>> plot_population_score(population)
   :param population: population object (:class:`GPopulation.GPopulation`)
   :param fitness: if True, the fitness score will be used, otherwise, the raw.
   :rtype: None
   """

   score_list = get_population_scores(population, fitness)
   plt.plot(score_list, 'o')
   plt.title("Plot of population score distribution")
   plt.xlabel('Individual')
   plt.ylabel('Score')
   plt.grid(True)
   plt.show()

# -----------------------------------------------------------------

def plot_histogram_population_score(population, fitness=False):

   """
   Population score distribution histogram
   Example:
      >>> plot_histogram_population_score(population)
   :param population: population object (:class:`GPopulation.GPopulation`)
   :param fitness: if True, the fitness score will be used, otherwise, the raw.
   :rtype: None
   """

   score_list = get_population_scores(population, fitness)
   n, bins, patches = plt.hist(score_list, 50, facecolor='green', alpha=0.75, normed=1)
   plt.plot(bins, mlab.normpdf(bins, numpy.mean(score_list), numpy.std(score_list)), 'r--')
   plt.xlabel('Score')
   plt.ylabel('Frequency')
   plt.grid(True)
   plt.title("Plot of population score distribution")
   plt.show()

# -----------------------------------------------------------------
