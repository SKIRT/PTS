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
from abc import ABCMeta, abstractproperty, abstractmethod
from math import sqrt as math_sqrt
from functools import partial

# Import other evolve modules
import pts.evolve.core.constants as constants
import pts.evolve.core.utils as utils
from pts.evolve.core.functionslot import FunctionSlot
from pts.evolve.core.statistics import Statistics
from pts.core.basics.containers import NamedList

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools.parallelization import MULTI_PROCESSING, Pool
from ...core.tools import strings
from ...core.tools import sequences

# -----------------------------------------------------------------

def key_raw_score(individual):

   """
   A key function to return raw score
   :param individual: the individual instance
   :rtype: the individual raw score
   .. note:: this function is used by the max()/min() python functions
   """

   return individual.score

# -----------------------------------------------------------------

def key_fitness_score(individual):

   """
   A key function to return fitness score, used by max()/min()
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

class PopulationBase(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, genome=None, **kwargs):

        """
        The constructor
        :param genome:
        :param kwargs:
        """

        # Set the genome
        self.oneSelfGenome = genome

        # The internal population representation
        self.internalPop = None
        self.internalPopRaw = None

        self.internalParams = {}
        self.multiProcessing = (False, False, None)

        # Statistics
        self.statted = False
        self.stats = Statistics()

        # Properties
        self.sorted = False
        self.minimax = kwargs.pop("minimax", constants.CDefPopMinimax)
        self.popSize = kwargs.pop("size", None)
        self.sortType = kwargs.pop("sort_type", constants.CDefPopSortType)
        self.scaleMethod = FunctionSlot("Scale Method")
        self.scaleMethod.set(kwargs.pop("scale_method", constants.CDefPopScale))
        self.allSlots = [self.scaleMethod]

    # -----------------------------------------------------------------

    @property
    def nindividuals(self):

        """
        This function ...
        :return:
        """

        return len(self.internalPop)

    # -----------------------------------------------------------------

    @abstractproperty
    def keys(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    def setMinimax(self, minimax):

        """
        Sets the population minimax
        Example:
         pop.setMinimax(Consts.minimaxType["maximize"])
        :param minimax: the minimax type
        """

        self.minimax = minimax

    # -----------------------------------------------------------------

    def __len__(self):

        """
        Return the length of population
        """

        return len(self.internalPop)

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

    @abstractproperty
    def individuals(self):

        """
        THis function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        Returns the iterator of the population
        """

        return iter(self.internalPop)

    # -----------------------------------------------------------------

    def __getitem__(self, key):

        """
        Returns the specified individual from population
        """

        return self.internalPop[key]

    # -----------------------------------------------------------------

    def __setitem__(self, key, value):

        """
        Set an individual of population
        """

        self.internalPop[key] = value
        self.clearFlags()

    # -----------------------------------------------------------------

    def statistics(self):

        """
        Do statistical analysis of population and set 'statted' to True
        """

        if self.statted: return
        log.debug("Running statistical calculations ...")
        raw_sum = 0

        len_pop = len(self)
        for ind in xrange(len_pop): raw_sum += self[ind].score

        # Set maximum, minimum and average
        self.stats["rawMax"] = max(self.individuals, key=key_raw_score).score
        self.stats["rawMin"] = min(self.individuals, key=key_raw_score).score
        self.stats["rawAve"] = raw_sum / float(len_pop)

        # Calculate the variance
        tmpvar = 0.0
        #for ind in xrange(len_pop):
        for ind in self:
            #s = self[ind].score - self.stats["rawAve"]
            s = ind.score - self.stats["rawAve"]
            s *= s
            tmpvar += s
        tmpvar /= float((len(self) - 1))

        # Set the standard deviation
        try: self.stats["rawDev"] = math_sqrt(tmpvar)
        except ValueError: self.stats["rawDev"] = 0.0

        # Set the variance
        self.stats["rawVar"] = tmpvar

        # Set statted flag
        self.statted = True

    # -----------------------------------------------------------------

    def scale(self, **args):

        """
        Scale the population using the scaling method
        :param args: this parameter is passed to the scale method
        """

        for it in self.scaleMethod.applyFunctions(self, **args): pass

        fit_sum = 0
        #for ind in xrange(len(self)): fit_sum += self[ind].fitness
        for ind in self: fit_sum += ind.fitness

        # Calculate max, min and average fitness
        self.stats["fitMax"] = max(self.individuals, key=key_fitness_score).fitness
        self.stats["fitMin"] = min(self.individuals, key=key_fitness_score).fitness
        self.stats["fitAve"] = fit_sum / float(len(self))

        # Set sorted flag to False
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

        """
        Copy current population to 'pop'
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

    @abstractmethod
    def clear(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def clone(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def clone_population(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @property
    def raw_sorting(self):

        """
        This function ...
        :return:
        """

        return self.sortType == constants.sortType["raw"]

    # -----------------------------------------------------------------

    @property
    def scaled_sorting(self):

        """
        This function ...
        :return:
        """

        return self.sortType == constants.sortType["scaled"]

    # -----------------------------------------------------------------

    def sort(self):

        """
        Sort the population
        """

        # Already sorted?
        if self.sorted: return

        #rev = (self.minimax == constants.minimaxType["maximize"])

        # Reverse or not
        if self.minimax == "minimize": rev = False
        elif self.minimax == "maximize": rev = True
        else: raise ValueError("Wrong minimax type: must be 'maximize' or 'minimize'")

        if self.raw_sorting: self.internalPop.sort(cmp=utils.cmp_individual_raw, reverse=rev)
        elif self.scaled_sorting:

            self.scale()
            self.internalPop.sort(cmp=utils.cmp_individual_scaled, reverse=rev)
            self.set_raw_from_internal()
            self.internalPopRaw.sort(cmp=utils.cmp_individual_raw, reverse=rev)

        else: raise ValueError("Invalid state of the sort type")

        # Set sorted flag
        self.sorted = True

    # -----------------------------------------------------------------

    def check_sorted(self):

        """
        This function ...
        :return:
        """

        # Reverse or not
        if self.minimax == "minimize": rev = False
        elif self.minimax == "maximize": rev = True
        else: raise ValueError("Wrong minimax type: must be 'maximize' or 'minimize'")

        if self.sortType == constants.sortType["raw"]:

            if not sequences.is_sorted([ind.score for ind in self.internalPop], invert=rev): raise RuntimeError("Not sorted")

        else:

            if not sequences.is_sorted([ind.fitness for ind in self.internalPop], invert=rev): raise RuntimeError("Not sorted")

            if not sequences.is_sorted([ind.score for ind in self.internalPopRaw], invert=rev): raise RuntimeError("Not sorted")

    # -----------------------------------------------------------------

    @abstractmethod
    def set_raw_from_internal(self):

        """
        This function ...
        :return: 
        """

        pass

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

    def check_population_size(self, size):

        """
        This function ...
        :param size: 
        :return: 
        """

        if self.popSize is None: raise ValueError("Population size hasn't been set yet")
        elif self.popSize != size: raise ValueError("The population size (" + str(self.popSize) + ") doesn't match the passed population size (" + str(size) + ")")

    # -----------------------------------------------------------------

    def set_or_check_population_size(self, size):

        """
        This function ...
        :param size: 
        :return: 
        """

        if self.popSize is None: self.popSize = size
        elif self.popSize != size: raise ValueError("The population has an incorrect number of individuals: is " + str(self.popSize) + " but should be " + str(size) + " according to genetic engine")

    # -----------------------------------------------------------------

    def setSortType(self, sort_type):

        """ Sets the sort type
        Example:
         pop.setSortType(Consts.sortType["scaled"])
        :param sort_type: the Sort Type
        """

        self.sortType = sort_type

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

    @abstractmethod
    def create(self):

        """
        This function ...
        :return: 
        """

        pass

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
            for gen in self.internalPop: gen.initialize(**kwargs)

        # Clear
        self.clearFlags()

    # -----------------------------------------------------------------

    def bestFitness(self, index=0):

        """
        Return the best scaled fitness individual of population
        :param index: the *index* best individual
        :rtype: the individual
        """

        self.sort()
        return self.internalPop[index]

    # -----------------------------------------------------------------

    @abstractmethod
    def best_fitness_key(self, index=0):

        """
        This function ...
        :param index:
        :return:
        """

        pass

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
                results = proc_pool.map(partial(multiprocessing_eval_full, **kwargs), self.individuals)
                proc_pool.close()
                proc_pool.join()
                for i in xrange(len(self)): self.internalPop[i] = results[i]

            else:

                #results = proc_pool.map(multiprocessing_eval, self.internalPop)
                results = proc_pool.map(partial(multiprocessing_eval, **kwargs), self.individuals)
                proc_pool.close()
                proc_pool.join()
                for individual, score in zip(self.individuals, results): individual.score = score

        else: # No multiprocessing: basically just a loop over evaluate() of the individuals

            # Evaluate each individual
            for ind in self.individuals: ind.evaluate(**kwargs)

        # Clear flags
        self.clearFlags()

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

    @abstractmethod
    def append(self, genome):

        """
        This function ...
        :param genome: 
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def __repr__(self):

        """
        THis function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def replace(self, key, new_key, value):

        """
        This function ...
        :param key: index or name
        :param new_key: index or name
        :param value:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def scores(self):

        """
        This function ...
        :return:
        """

        return [ind.score for ind in self]

    # -----------------------------------------------------------------

    @property
    def fitnesses(self):

        """
        This function ...
        :return:
        """

        return [ind.fitness for ind in self]

# -----------------------------------------------------------------

class NamedPopulation(PopulationBase):

    """
    The NamedPopulation class: represents a population where each individual has a name
    """

    def __init__(self, genome=None, **kwargs):

        """
        This function ...
        :param genome: 
        :param kwargs:
        """

        # Call the constructor of the base class
        super(NamedPopulation, self).__init__(genome, **kwargs)

        # The containers of individuals
        self.internalPop = NamedList()
        self.internalPopRaw = NamedList()

        # Create name iterator
        self.name_iterator = strings.alphabet_strings_iterator()

    # -----------------------------------------------------------------

    @classmethod
    def from_population(cls, population):

        """
        This function ...
        :param population: 
        :return: 
        """

        # Create new population
        pop = cls()

        # Set attributes
        pop.oneSelfGenome = population.oneSelfGenome
        pop.internalPop = NamedList()
        pop.internalPopRaw = NamedList()
        pop.popSize = population.popSize
        pop.sortType = population.sortType
        pop.sorted = False
        pop.minimax = population.minimax
        pop.scaleMethod = population.scaleMethod
        pop.allSlots = [pop.scaleMethod]

        pop.internalParams = population.internalParams
        pop.multiProcessing = population.multiProcessing

        pop.statted = False
        pop.stats = Statistics()

        # Set the state of the name iterator
        pop.name_iterator = population.name_iterator

        # Return the population
        return pop

    # -----------------------------------------------------------------

    @property
    def individuals(self):

        """
        This function ...
        :return: 
        """

        return self.internalPop.values

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return: 
        """

        return self.internalPop.names

    # -----------------------------------------------------------------

    @property
    def items(self):

        """
        This function ...
        :return: 
        """

        return self.internalPop.items

    # -----------------------------------------------------------------

    def clone(self):

        """
        Return a brand-new cloned population
        """

        newpop = NamedPopulation(self.oneSelfGenome)
        self.copy(newpop)
        return newpop

    # -----------------------------------------------------------------

    def clone_population(self):

        """
        This function ...
        :return: 
        """

        return NamedPopulation.from_population(self)

    # -----------------------------------------------------------------

    def set_raw_from_internal(self):

        """
        THis function ...
        :return: 
        """

        self.internalPopRaw = self.internalPop.copy()

    # -----------------------------------------------------------------

    def create(self, **args):

        """
        Clone the example genome to fill the population
        """

        # Set minimax attribute
        self.minimax = args["minimax"]

        # Generate individuals
        for _ in range(self.popSize):

            # Generate name
            name = self.name_iterator.next()

            # Generate genome
            genome = self.oneSelfGenome.clone()

            # Add to the population
            self.append(genome, name=name)

        # Clear all flags
        self.clearFlags()

    # -----------------------------------------------------------------

    def append(self, genome, name=None):

        """
        This function ...
        :param genome: 
        :param name:
        :return: 
        """

        # Get the unique name for this individual
        if name is None: name = self.name_iterator.next()
        elif name in self.names: raise ValueError("Already an individual with this name")

        # Add to the internal population
        self.internalPop.append(name, genome)

        # Return the name
        return name

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return: 
        """

        self.internalPop = NamedList()
        self.internalPopRaw = NamedList()
        self.clearFlags()

    # -----------------------------------------------------------------

    @property
    def keys(self):

        """
        This function ...
        :return: 
        """

        return self.names

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        Returns the string representation of the population
        """

        ret = "- NamedPopulation\n"
        ret += "\tPopulation Size:\t %d\n" % (self.popSize,)
        ret += "\tSort Type:\t\t %s\n" % (constants.sortType.keys()[constants.sortType.values().index(self.sortType)].capitalize(),)
        #ret += "\tMinimax Type:\t\t %s\n" % (constants.minimaxType.keys()[constants.minimaxType.values().index(self.minimax)].capitalize(),)
        ret += "\tMinimax Type:\t\t %s\n" % (self.minimax.capitalize(),)
        for slot in self.allSlots:
         ret += "\t" + slot.__repr__()
        ret += "\n"
        ret += self.stats.__repr__()
        return ret

    # -----------------------------------------------------------------

    def best_fitness_key(self, index=0):

        """
        This function ...
        :param index:
        :return:
        """

        self.sort()
        return self.names[index]

    # -----------------------------------------------------------------

    def replace(self, key, new_key, value):

        """
        This function ...
        :param key:
        :param new_key:
        :param value:
        :return:
        """

        # Replace
        return self.internalPop.replace(key, value, new_name=new_key)

# -----------------------------------------------------------------

class Population(PopulationBase):

    """
    Population Class - The container for the population

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

    def __init__(self, genome=None, **kwargs):

        """
        The constructor ...
        :param genome:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Population, self).__init__(genome, **kwargs)

        # The containers of individuals
        self.internalPop = []
        self.internalPopRaw = []

    # -----------------------------------------------------------------

    @classmethod
    def from_population(cls, population):

        """
        This function ...
        :param population: 
        :return: 
        """

        # Create new population
        pop = cls()

        # Set attributes
        pop.oneSelfGenome = population.oneSelfGenome
        pop.internalPop = []
        pop.internalPopRaw = []
        pop.popSize = population.popSize
        pop.sortType = population.sortType
        pop.sorted = False
        pop.minimax = population.minimax
        pop.scaleMethod = population.scaleMethod
        pop.allSlots = [pop.scaleMethod]

        pop.internalParams = population.internalParams
        pop.multiProcessing = population.multiProcessing

        pop.statted = False
        pop.stats = Statistics()

        # Return the population
        return pop

    # -----------------------------------------------------------------

    @property
    def keys(self):

        """
        This function ...
        :return: 
        """

        return range(len(self))

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        Returns the string representation of the population
        """

        ret = "- Population\n"
        ret += "\tPopulation Size:\t %d\n" % (self.popSize,)
        ret += "\tSort Type:\t\t %s\n" % (constants.sortType.keys()[constants.sortType.values().index(self.sortType)].capitalize(),)
        #ret += "\tMinimax Type:\t\t %s\n" % (constants.minimaxType.keys()[constants.minimaxType.values().index(self.minimax)].capitalize(),)
        ret += "\tMinimax Type:\t\t %s\n" % (self.minimax.capitalize(),)
        for slot in self.allSlots:
         ret += "\t" + slot.__repr__()
        ret += "\n"
        ret += self.stats.__repr__()
        return ret

    # -----------------------------------------------------------------

    @property
    def individuals(self):

        """
        This function ...
        :return: 
        """

        return self.internalPop

    # -----------------------------------------------------------------

    def set_raw_from_internal(self):

        """
        This function ...
        :return: 
        """

        self.internalPopRaw = self.internalPop[:]

    # -----------------------------------------------------------------

    def create(self, **args):

        """
        Clone the example genome to fill the population
        """

        self.minimax = args["minimax"]
        self.internalPop = [self.oneSelfGenome.clone() for i in xrange(self.popSize)]
        self.clearFlags()

    # -----------------------------------------------------------------

    def append(self, genome):

        """
        This function ...
        :param genome: 
        :return: 
        """

        self.internalPop.append(genome)
        return self.nindividuals - 1

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

    def clone_population(self):

        """
        This function ...
        :return: 
        """

        return Population.from_population(self)

    # -----------------------------------------------------------------

    def best_fitness_key(self, index=0):

        """
        This function ...
        :param index:
        :return:
        """

        self.sort()
        return index

    # -----------------------------------------------------------------

    def replace(self, key, new_key, value):

        """
        This function ...
        :param key:
        :param value:
        :return
        """

        # Remove the key'th element
        old = self.internalPop.pop(key)

        # Inser in the population list
        self.internalPop.insert(new_key, value)

        # Return the old individual
        return old

# -----------------------------------------------------------------
