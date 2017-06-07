#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.selectors This module have the *selection methods*, like roulette wheel, tournament, ranking,
#  etc.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from ...core.tools.random import prng

# -----------------------------------------------------------------

def GRankSelector(population, **args):

    """
    The Rank Selector - This selector will pick the best individual of
    the population every time.
    """

    from . import constants

    return_key = args.pop("return_key", False)

    count = 0

    if args["popID"] != GRankSelector.cachePopID:

      if population.sortType == constants.sortType["scaled"]:

         best_fitness = population.bestFitness().fitness
         for index in xrange(1, len(population.internalPop)):
            if population[index].fitness == best_fitness:
               count += 1
      else:

         best_raw = population.bestRaw().score
         for index in xrange(1, len(population.internalPop)):
            if population[index].score == best_raw:
               count += 1

      GRankSelector.cachePopID = args["popID"]
      GRankSelector.cacheCount = count

    else: count = GRankSelector.cacheCount

    if count == 0: index = 0 #individual = population[0]
    else: index = prng.randint(0, count + 1) #individual = population[prng.randint(0, count + 1)] # HERE IT SHOULD BE INCLUSIVE

    if return_key: return population.keys[index]
    else: return population[index]

    #return individual

# -----------------------------------------------------------------

GRankSelector.cachePopID = None
GRankSelector.cacheCount = None

# -----------------------------------------------------------------

def GUniformSelector(population, **args):

   """
   The Uniform Selector
   """

   return_key = args.pop("return_key", False)

   random_index = prng.randint(0, len(population))

   if return_key: return population.keys[random_index]
   else: return population[random_index]

# -----------------------------------------------------------------

def GTournamentSelector(population, **args):

   """
   The Tournament Selector
   It accepts the *tournamentPool* population parameter.
   .. note::
      the Tournament Selector uses the Roulette Wheel to
      pick individuals for the pool
   .. versionchanged:: 0.6
      Changed the parameter `poolSize` to the `tournamentPool`, now the selector
      gets the pool size from the population.
   """

   from . import constants

   return_key = args.pop("return_key", False)

   choosen = None
   should_minimize = population.minimax == "minimize" #constants.minimaxType["minimize"]
   minimax_operator = min if should_minimize else max

   poolSize = population.getParam("tournamentPool", constants.CDefTournamentPoolSize)

   # Pick individuals for the tournament pool
   args["return_key"] = return_key
   #tournament_pool = [GRouletteWheel(population, **args) for i in xrange(poolSize)]
   tournament_keys = [GRouletteWheel(population, **args) for i in xrange(poolSize)]

   if population.sortType == constants.sortType["scaled"]:
      #choosen = minimax_operator(tournament_pool, key=lambda ind: ind.fitness)
      choosen_key = minimax_operator(tournament_keys, key=lambda key: population[key].fitness)
   else:
      #choosen = minimax_operator(tournament_pool, key=lambda ind: ind.score)
      choosen_key = minimax_operator(tournament_keys, key=lambda key: population[key].score)

   # Return the choosen one
   #return choosen
   if return_key:  return choosen_key
   else: return population[choosen_key]

# -----------------------------------------------------------------

def GTournamentSelectorAlternative(population, **args):

   """
   The alternative Tournament Selector
   This Tournament Selector doesn't use the Roulette Wheel
   It accepts the *tournamentPool* population parameter.
   .. versionadded: 0.6
      Added the GTournamentAlternative function.
   """

   from . import constants

   return_key = args.pop("return_key", False)

   pool_size = population.getParam("tournamentPool", constants.CDefTournamentPoolSize)
   len_pop = len(population)

   should_minimize = population.minimax == "minimize" #constants.minimaxType["minimize"]

   minimax_operator = min if should_minimize else max

   # Pick random individuals to form the pool
   # Pick random indices
   random_indices = [prng.randint(0, len_pop) for i in xrange(pool_size)]
   tournament_pool = [population[ri] for ri in random_indices]

   if population.sortType == constants.sortType["scaled"]:
      #choosen = minimax_operator(tournament_pool, key=lambda ind: ind.fitness)
      choosen_index = minimax_operator(range(len(tournament_pool)), key=lambda index: tournament_pool[index].fitness)
   else:
      #choosen = minimax_operator(tournament_pool, key=lambda ind: ind.score)
      choosen_index = minimax_operator(range(len(tournament_pool)), key=lambda index: tournament_pool[index].score)

   choosen_individual_index = random_indices[choosen_index]

   if return_key: return population.keys[choosen_individual_index]
   else: return population[choosen_individual_index]

   #return choosen

# -----------------------------------------------------------------

def GRouletteWheel(population, **args):

   """
   The Roulette Wheel selector
   """

   return_key = args.pop("return_key", False)

   psum = None
   if args["popID"] != GRouletteWheel.cachePopID:
      GRouletteWheel.cachePopID = args["popID"]
      psum = GRouletteWheel_PrepareWheel(population)
      GRouletteWheel.cacheWheel = psum
   else:
      psum = GRouletteWheel.cacheWheel

   cutoff = prng.random_sample()
   lower = 0
   upper = len(population) - 1
   while(upper >= lower):
      i = lower + ((upper - lower) / 2)
      if psum[i] > cutoff:
         upper = i - 1
      else:
         lower = i + 1

   lower = min(len(population) - 1, lower)
   lower = max(0, lower)

   # Return key or individual itself
   if return_key: return population.best_fitness_key(lower)
   else: return population.bestFitness(lower)

# -----------------------------------------------------------------

GRouletteWheel.cachePopID = None
GRouletteWheel.cacheWheel = None

# -----------------------------------------------------------------

def GRouletteWheel_PrepareWheel(population):

   """
   A preparation for Roulette Wheel selection
   """

   from . import constants

   len_pop = len(population)

   psum = [i for i in xrange(len_pop)]

   population.statistics()

   if population.sortType == constants.sortType["scaled"]:

      pop_fitMax = population.stats["fitMax"]
      pop_fitMin = population.stats["fitMin"]

      if pop_fitMax == pop_fitMin:
         for index in xrange(len_pop):
            psum[index] = (index + 1) / float(len_pop)
      elif (pop_fitMax > 0 and pop_fitMin >= 0) or (pop_fitMax <= 0 and pop_fitMin < 0):

         population.sort()

         #if population.minimax == constants.minimaxType["maximize"]:
         if population.minimax == "maximize":

            psum[0] = population[0].fitness
            for i in xrange(1, len_pop):
               psum[i] = population[i].fitness + psum[i - 1]
            for i in xrange(len_pop):
               psum[i] /= float(psum[len_pop - 1])

         #else:
         elif population.minimax == "minimize":

            psum[0] = -population[0].fitness + pop_fitMax + pop_fitMin
            for i in xrange(1, len_pop):
               psum[i] = -population[i].fitness + pop_fitMax + pop_fitMin + psum[i - 1]
            for i in xrange(len_pop):
               psum[i] /= float(psum[len_pop - 1])

         # Invalid
         else: raise ValueError("Invalid minimax: " + str(population.minimax))

   else:

      pop_rawMax = population.stats["rawMax"]
      pop_rawMin = population.stats["rawMin"]

      if pop_rawMax == pop_rawMin:
         for index in xrange(len_pop):
            psum[index] = (index + 1) / float(len_pop)

      elif (pop_rawMax > 0 and pop_rawMin >= 0) or (pop_rawMax <= 0 and pop_rawMin < 0):

         population.sort()

         #if population.minimax == constants.minimaxType["maximize"]:
         if population.minimax == "maximize":

            psum[0] = population[0].score
            for i in xrange(1, len_pop):
               psum[i] = population[i].score + psum[i - 1]
            for i in xrange(len_pop):
               psum[i] /= float(psum[len_pop - 1])

         #else:
         elif population.minimax == "minimize":

            psum[0] = - population[0].score + pop_rawMax + pop_rawMin
            for i in xrange(1, len_pop):
               psum[i] = - population[i].score + pop_rawMax + pop_rawMin + psum[i - 1]
            for i in xrange(len_pop):
               psum[i] /= float(psum[len_pop - 1])

         else: raise ValueError("Invalid minimax: " + str(population.minimax))

   return psum

# -----------------------------------------------------------------
