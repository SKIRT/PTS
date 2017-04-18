#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.population This module contains the :class:`GPopulation.GPopulation` class, which is reponsible
#  to keep the population and the statistics.
#
# Default Parameters
# -------------------------------------------------------------
#  *Sort Type*
#    >>> Consts.sortType["scaled"]
#    The scaled sort type
#  *Minimax*
#    >>> Consts.minimaxType["maximize"]
#    Maximize the evaluation function
#  *Scale Method*
#    :func:`Scaling.LinearScaling`
#    The Linear Scaling scheme
#

# -----------------------------------------------------------------

# Import standard modules
from math import sqrt as math_sqrt
from functools import partial

# Import other evolve modules
import pts.evolve.core.constants as constants
import pts.evolve.core.utils as utils
from pts.evolve.core.functionslot import FunctionSlot
from pts.evolve.core.statistics import Statistics

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

try:
   from multiprocessing import cpu_count, Pool
   CPU_COUNT = cpu_count()
   MULTI_PROCESSING = True if CPU_COUNT > 1 else False
   log.debug("You have %d CPU cores, so the multiprocessing state is %s", CPU_COUNT, MULTI_PROCESSING)
except ImportError:
   MULTI_PROCESSING = False
   log.debug("You don't have multiprocessing support")

# -----------------------------------------------------------------

def key_raw_score(individual):

   """ A key function to return raw score
   :param individual: the individual instance
   :rtype: the individual raw score
   .. note:: this function is used by the max()/min() python functions
   """

   return individual.score

# -----------------------------------------------------------------

def key_fitness_score(individual):

   """ A key function to return fitness score, used by max()/min()
   :param individual: the individual instance
   :rtype: the individual fitness score
   .. note:: this function is used by the max()/min() python functions
   """

   return individual.fitness

# -----------------------------------------------------------------

def multiprocessing_eval(ind, **kwargs):

   """
   Internal used by the multiprocessing
   """

   ind.evaluate(**kwargs)
   return ind.score

# -----------------------------------------------------------------

def multiprocessing_eval_full(ind, **kwargs):

   """
   Internal used by the multiprocessing (full copy)
   """

   ind.evaluate(**kwargs)
   return ind

# -----------------------------------------------------------------

class Population(object):

    """ Population Class - The container for the population

    **Examples**
      Get the population from the :class:`GSimpleGA.GSimpleGA` (GA Engine) instance
         pop = ga_engine.get_population()

      Get the best fitness individual
         bestIndividual = pop.bestFitness()

      Get the best raw individual
         bestIndividual = pop.bestRaw()

      Get the statistics from the :class:`Statistics.Statistics` instance
         stats = pop.getStatistics()
         print stats["rawMax"]
         10.4

      Iterate, get/set individuals
         for ind in pop:
            print ind
         (...)

         for i in xrange(len(pop)):
            print pop[i]
         (...)

         pop[10] = newGenome
         pop[10].fitness
         12.5

    :param genome: the :term:`Sample genome`, or a GPopulation object, when cloning.

    """

    def __init__(self, genome):

        """
        The constructor ...
        :param genome:
        """

        if isinstance(genome, Population):

            self.oneSelfGenome = genome.oneSelfGenome
            self.internalPop = []
            self.internalPopRaw = []
            self.popSize = genome.popSize
            self.sortType = genome.sortType
            self.sorted = False
            self.minimax = genome.minimax
            self.scaleMethod = genome.scaleMethod
            self.allSlots = [self.scaleMethod]

            self.internalParams = genome.internalParams
            self.multiProcessing = genome.multiProcessing

            self.statted = False
            self.stats = Statistics()

        else:

            # Debugging
            log.debug("New population instance, %s class genomes", genome.__class__.__name__)

            self.oneSelfGenome = genome
            self.internalPop = []
            self.internalPopRaw = []
            self.popSize = 0
            self.sortType = constants.CDefPopSortType
            self.sorted = False
            self.minimax = constants.CDefPopMinimax
            self.scaleMethod = FunctionSlot("Scale Method")
            self.scaleMethod.set(constants.CDefPopScale)
            self.allSlots = [self.scaleMethod]

            self.internalParams = {}
            self.multiProcessing = (False, False, None)

            # Statistics
            self.statted = False
            self.stats = Statistics()

    # -----------------------------------------------------------------

    def setMultiProcessing(self, flag=True, full_copy=False, max_processes=None):

        """
        This function sets the flag to enable/disable the use of python multiprocessing module.
        Use this option when you have more than one core on your CPU and when your
        evaluation function is very slow.
        The parameter "full_copy" defines where the individual data should be copied back
        after the evaluation or not. This parameter is useful when you change the
        individual in the evaluation function.

        :param flag: True (default) or False
        :param full_copy: True or False (default)
        :param max_processes: None (default) or an integer value

        .. warning:: Use this option only when your evaluation function is slow, se you
                   will get a good tradeoff between the process communication speed and the
                   parallel evaluation.

        .. versionadded:: 0.6
         The `setMultiProcessing` method.

        """
        self.multiProcessing = (flag, full_copy, max_processes)

    # -----------------------------------------------------------------

    def setMinimax(self, minimax):

        """ Sets the population minimax
        Example:
         pop.setMinimax(Consts.minimaxType["maximize"])
        :param minimax: the minimax type
        """

        self.minimax = minimax

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        Returns the string representation of the population
        """

        ret = "- Population\n"
        ret += "\tPopulation Size:\t %d\n" % (self.popSize,)
        ret += "\tSort Type:\t\t %s\n" % (constants.sortType.keys()[constants.sortType.values().index(self.sortType)].capitalize(),)
        ret += "\tMinimax Type:\t\t %s\n" % (constants.minimaxType.keys()[constants.minimaxType.values().index(self.minimax)].capitalize(),)
        for slot in self.allSlots:
         ret += "\t" + slot.__repr__()
        ret += "\n"
        ret += self.stats.__repr__()
        return ret

    # -----------------------------------------------------------------

    def __len__(self):

        """
        Return the length of population
        """

        return len(self.internalPop)

    # -----------------------------------------------------------------

    def __getitem__(self, key):

        """
        Returns the specified individual from population
        """

        return self.internalPop[key]

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        Returns the iterator of the population
        """

        return iter(self.internalPop)

    # -----------------------------------------------------------------

    def __setitem__(self, key, value):

        """
        Set an individual of population
        """

        self.internalPop[key] = value
        self.clearFlags()

    # -----------------------------------------------------------------

    def clearFlags(self):

        """
        Clear the sorted and statted internal flags
        """

        self.sorted = False
        self.statted = False

    # -----------------------------------------------------------------

    def getStatistics(self):

        """
        Return a Statistics class for statistics
        :rtype: the :class:`Statistics.Statistics` instance
        """

        self.statistics()
        return self.stats

    # -----------------------------------------------------------------

    def statistics(self):

        """
        Do statistical analysis of population and set 'statted' to True
        """

        if self.statted:
         return
        log.debug("Running statistical calculations ...")
        raw_sum = 0

        len_pop = len(self)
        for ind in xrange(len_pop):
         raw_sum += self[ind].score

        self.stats["rawMax"] = max(self, key=key_raw_score).score
        self.stats["rawMin"] = min(self, key=key_raw_score).score
        self.stats["rawAve"] = raw_sum / float(len_pop)

        tmpvar = 0.0
        for ind in xrange(len_pop):
         s = self[ind].score - self.stats["rawAve"]
         s *= s
         tmpvar += s

        tmpvar /= float((len(self) - 1))

        try: self.stats["rawDev"] = math_sqrt(tmpvar)
        except ValueError: self.stats["rawDev"] = 0.0

        self.stats["rawVar"] = tmpvar

        self.statted = True

    # -----------------------------------------------------------------

    def bestFitness(self, index=0):

        """ Return the best scaled fitness individual of population
        :param index: the *index* best individual
        :rtype: the individual
        """

        self.sort()
        return self.internalPop[index]

    # -----------------------------------------------------------------

    def worstFitness(self):

        """
        Return the worst scaled fitness individual of the population
        :rtype: the individual
        """

        self.sort()
        return self.internalPop[-1]

    # -----------------------------------------------------------------

    def bestRaw(self, index=0):

        """
        Return the best raw score individual of population
        :param index: the *index* best raw individual
        :rtype: the individual
        .. versionadded:: 0.6
         The parameter `index`.
        """

        if self.sortType == constants.sortType["raw"]: return self.internalPop[index]
        else:

            self.sort()
            return self.internalPopRaw[index]

    # -----------------------------------------------------------------

    def worstRaw(self):

        """
        Return the worst raw score individual of population
        :rtype: the individual
        .. versionadded:: 0.6
         The parameter `index`.
        """

        if self.sortType == constants.sortType["raw"]: return self.internalPop[-1]
        else:
            self.sort()
            return self.internalPopRaw[-1]

    # -----------------------------------------------------------------

    def sort(self):

        """
        Sort the population
        """

        if self.sorted:
         return
        rev = (self.minimax == constants.minimaxType["maximize"])

        if self.sortType == constants.sortType["raw"]:
         self.internalPop.sort(cmp=utils.cmp_individual_raw, reverse=rev)
        else:
         self.scale()
         self.internalPop.sort(cmp=utils.cmp_individual_scaled, reverse=rev)
         self.internalPopRaw = self.internalPop[:]
         self.internalPopRaw.sort(cmp=utils.cmp_individual_raw, reverse=rev)

        self.sorted = True

    # -----------------------------------------------------------------

    def set_ranges(self, minima, maxima):

        """
        This function ...
        :param minima:
        :param maxima:
        :return:
        """

        # Create dictionary
        params = {"minima": minima, "maxima": maxima}

        # Update the parameters 'minima' and 'maxima' for each individual
        self.set_params_for_all_individuals(**params)

    # -----------------------------------------------------------------

    def setPopulationSize(self, size):

        """
        This function sets the population size
        :param size: the population size
        """

        self.popSize = size

    # -----------------------------------------------------------------

    def setSortType(self, sort_type):

        """ Sets the sort type
        Example:
         pop.setSortType(Consts.sortType["scaled"])
        :param sort_type: the Sort Type
        """

        self.sortType = sort_type

    # -----------------------------------------------------------------

    def create(self, **args):

        """
        Clone the example genome to fill the population
        """

        self.minimax = args["minimax"]
        self.internalPop = [self.oneSelfGenome.clone() for i in xrange(self.popSize)]
        self.clearFlags()

    # -----------------------------------------------------------------

    def __findIndividual(self, individual, end):

        """
        This function ...
        :param individual:
        :param end:
        :return:
        """

        for i in xrange(end):
         if individual.compare(self.internalPop[i]) == 0:
            return True

    # -----------------------------------------------------------------

    def set_params_for_all_individuals(self, **params):

        """
        This function ...
        :param params
        :return:
        """

        # Loop over all individuals in the internal population
        for i in xrange(len(self.internalPop)):

            # Set the parameters of this individual
            curr = self.internalPop[i]
            curr.setParams(**params)

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        Initialize all individuals of population,
        this calls the initialize() of individuals
        :param kwargs: arguemnts passed to the initialize function of the individuals
        """

        # Inform the user
        log.info("Initializing the population ...")

        if self.oneSelfGenome.getParam("full_diversity", True) and hasattr(self.oneSelfGenome, "compare"):

         for i in xrange(len(self.internalPop)):
            curr = self.internalPop[i]
            curr.initialize(**kwargs)
            while self.__findIndividual(curr, i):
               curr.initialize(**kwargs)

        else:

         for gen in self.internalPop:
            gen.initialize(**kwargs)

        # Clear
        self.clearFlags()

    # -----------------------------------------------------------------

    def evaluate(self, silent, **kwargs):

        """
        Evaluate all individuals in population, calls the evaluate() method of individuals
        :param silent:
        :param kwargs: this params are passed to the evaluation function
        """

        # Inform the user
        if not silent: log.info("Evaluating the new population ...")

        # We have multiprocessing
        if self.multiProcessing[0] and MULTI_PROCESSING:

            log.debug("Evaluating the population using the multiprocessing method")
            proc_pool = Pool(processes=self.multiProcessing[2])

            # Multiprocessing full_copy parameter
            if self.multiProcessing[1]:

                #results = proc_pool.map(multiprocessing_eval_full, self.internalPop)
                results = proc_pool.map(partial(multiprocessing_eval_full, **kwargs), self.internalPop)
                proc_pool.close()
                proc_pool.join()
                for i in xrange(len(self.internalPop)): self.internalPop[i] = results[i]

            else:

                #results = proc_pool.map(multiprocessing_eval, self.internalPop)
                results = proc_pool.map(partial(multiprocessing_eval, **kwargs), self.internalPop)
                proc_pool.close()
                proc_pool.join()
                for individual, score in zip(self.internalPop, results): individual.score = score

        else: # No multiprocessing: basically just a loop over evaluate() of the individuals

            # Evaluate each individual
            for ind in self.internalPop: ind.evaluate(**kwargs)

        # Clear flags
        self.clearFlags()

    # -----------------------------------------------------------------

    def scale(self, **args):

        """
        Scale the population using the scaling method
        :param args: this parameter is passed to the scale method
        """

        for it in self.scaleMethod.applyFunctions(self, **args):
         pass

        fit_sum = 0
        for ind in xrange(len(self)):
         fit_sum += self[ind].fitness

        self.stats["fitMax"] = max(self, key=key_fitness_score).fitness
        self.stats["fitMin"] = min(self, key=key_fitness_score).fitness
        self.stats["fitAve"] = fit_sum / float(len(self))

        self.sorted = False

    # -----------------------------------------------------------------

    def printStats(self):

        """
        Print statistics of the current population
        """

        message = ""
        if self.sortType == constants.sortType["scaled"]:
         message = "Max/Min/Avg Fitness(Raw) [%(fitMax).2f(%(rawMax).2f)/%(fitMin).2f(%(rawMin).2f)/%(fitAve).2f(%(rawAve).2f)]" % self.stats
        else:
         message = "Max/Min/Avg Raw [%(rawMax).2f/%(rawMin).2f/%(rawAve).2f]" % self.stats
        log.info(message)
        print message
        return message

    # -----------------------------------------------------------------

    def copy(self, pop):

        """ Copy current population to 'pop'
        :param pop: the destination population
        .. warning:: this method do not copy the individuals, only the population logic
        """

        pop.popSize = self.popSize
        pop.sortType = self.sortType
        pop.minimax = self.minimax
        pop.scaleMethod = self.scaleMethod
        pop.internalParams = self.internalParams
        pop.multiProcessing = self.multiProcessing

    # -----------------------------------------------------------------

    def getParam(self, key, nvl=None):

        """
        Gets an internal parameter
        Example:
         population.getParam("tournamentPool")
         5
        :param key: the key of param
        :param nvl: if the key doesn't exist, the nvl will be returned
        """

        return self.internalParams.get(key, nvl)

    # -----------------------------------------------------------------

    def get_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.internalParams.copy()

    # -----------------------------------------------------------------

    def setParams(self, **args):

        """
        Sets an internal parameter
        Example:
         population.setParams(tournamentPool=5)
        :param args: parameters to set
        .. versionadded:: 0.6
         The `setParams` method.
        """

        self.internalParams.update(args)

    # -----------------------------------------------------------------

    def clear(self):

        """
        Remove all individuals from population
        """

        del self.internalPop[:]
        del self.internalPopRaw[:]
        self.clearFlags()

    # -----------------------------------------------------------------

    def clone(self):

        """
        Return a brand-new cloned population
        """

        newpop = Population(self.oneSelfGenome)
        self.copy(newpop)
        return newpop

# -----------------------------------------------------------------
