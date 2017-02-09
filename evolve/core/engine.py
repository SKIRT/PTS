#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.engine This module contains the GA Engine, the GA Engine class is responsible
#  for all the evolutionary process. It contains the GA Engine related
#  functions, like the Termination Criteria functions for convergence analysis, etc.
#
#  Default Parameters:
#*Number of Generations*
#   Default is 100 generations
#*Mutation Rate*
#   Default is 0.02, which represents 2%
#*Crossover Rate*
#   Default is 0.9, which represents 90%
#*Elitism Replacement*
#   Default is 1 individual
#*Population Size*
#   Default is 80 individuals
#*Minimax*
#   >>> constants.minimaxType["maximize"]
#   Maximize the evaluation function
#*DB Adapter*
#   Default is **None**
#*Migration Adapter*
#   Default is **None**
#*Interactive Mode*
#   Default is **True**
#*Selector (Selection Method)*
#   :func:`Selectors.GRankSelector`
#   The Rank Selection method
#

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import division, print_function

# Import standard modules
from copy import deepcopy
import numpy as np
from time import time
from types import BooleanType
from sys import stdout as sys_stdout

# Import other evolve modules
from pts.evolve.core.population import Population
from pts.evolve.core.functionslot import FunctionSlot
from pts.evolve.core.genome import GenomeBase
from pts.evolve.core.dbadapters import DBBaseAdapter
import pts.evolve.core.constants as constants
import pts.evolve.core.utils as utils

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import serialization
from ...core.tools.random import prng

# -----------------------------------------------------------------

def RawScoreCriteria(ga_engine):

    """
    Terminate the evolution using the **bestrawscore** and **rounddecimal**
    parameter obtained from the individual
    Example:
       genome.setParams(bestrawscore=0.00, rounddecimal=2)
       (...)
       >>> ga_engine.terminationCriteria.set(RawScoreCriteria)
    """

    ind = ga_engine.bestIndividual()
    bestRawScore = ind.getParam("bestrawscore")
    roundDecimal = ind.getParam("rounddecimal")

    if bestRawScore is None:
        utils.raiseException("you must specify the bestrawscore parameter", ValueError)

    if ga_engine.getMinimax() == constants.minimaxType["maximize"]:
        if roundDecimal is not None:
            return round(bestRawScore, roundDecimal) <= round(ind.score, roundDecimal)
        else:
            return bestRawScore <= ind.score
    else:
        if roundDecimal is not None:
            return round(bestRawScore, roundDecimal) >= round(ind.score, roundDecimal)
        else:
            return bestRawScore >= ind.score

# -----------------------------------------------------------------

def ConvergenceCriteria(ga_engine):

    """
    Terminate the evolution when the population have converged
    Example: ga_engine.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
    """

    pop = ga_engine.get_population()
    return pop[0] == pop[len(pop) - 1]

# -----------------------------------------------------------------

def RawStatsCriteria(ga_engine):

    """
    Terminate the evolution based on the raw stats
    Example: ga_engine.terminationCriteria.set(GSimpleGA.RawStatsCriteria)
    """

    stats = ga_engine.getStatistics()
    if stats["rawMax"] == stats["rawMin"]:
        if stats["rawAve"] == stats["rawMax"]:
            return True
    return False

# -----------------------------------------------------------------

def FitnessStatsCriteria(ga_engine):

    """
    Terminate the evoltion based on the fitness stats
    Example: ga_engine.terminationCriteria.set(GSimpleGA.FitnessStatsCriteria)
    """

    stats = ga_engine.getStatistics()
    if stats["fitMax"] == stats["fitMin"]:
        if stats["fitAve"] == stats["fitMax"]:
            return True
    return False

# -----------------------------------------------------------------

class GeneticEngine(object):

    """
    This class represents the Genetic Engine

    Example:
    # ga = GeneticEngine(genome)
    # ga.selector.set(Selectors.GRouletteWheel)
    # ga.setGenerations(120)
    # ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

    :param genome: the :term:`Sample Genome`
    :param interactiveMode: this flag enables the Interactive Mode, the default is True
    .. note:: if you use the same random seed, all the runs of algorithm will be the same
    """

    selector = None
    """ This is the function slot for the selection method
    if you want to change the default selector, you must do this: ::

       ga_engine.selector.set(Selectors.GRouletteWheel) """

    stepCallback = None
    """ This is the :term:`step callback function` slot,
    if you want to set the function, you must do this: ::

       def your_func(ga_engine):
          # Here you have access to the GA Engine
          return False

       ga_engine.stepCallback.set(your_func)

    now *"your_func"* will be called every generation.
    When this function returns True, the GA Engine will stop the evolution and show
    a warning, if False, the evolution continues.
    """

    terminationCriteria = None
    """ This is the termination criteria slot, if you want to set one
    termination criteria, you must do this: ::

       ga_engine.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

    Now, when you run your GA, it will stop when the population converges.

    There are those termination criteria functions: :func:`GSimpleGA.RawScoreCriteria`,
    :func:`GSimpleGA.ConvergenceCriteria`, :func:`GSimpleGA.RawStatsCriteria`, :func:`GSimpleGA.FitnessStatsCriteria`

    But you can create your own termination function, this function receives
    one parameter which is the GA Engine, follows an example: ::

       def ConvergenceCriteria(ga_engine):
          pop = ga_engine.get_population()
          return pop[0] == pop[len(pop)-1]

    When this function returns True, the GA Engine will stop the evolution and show
    a warning, if False, the evolution continues, this function is called every
    generation.
    """

    def __init__(self, genome, interactive=True):

        """
        Initializator of GSimpleGA
        :param genome:
        :param interactive:
        """

        if type(interactive) != BooleanType:
            utils.raiseException("Interactive Mode option must be True or False", TypeError)

        if not isinstance(genome, GenomeBase):
            utils.raiseException("The genome must be a GenomeBase subclass", TypeError)

        # Initialize things
        self.internalPop = Population(genome)
        self.nGenerations = constants.CDefGAGenerations
        self.pMutation = constants.CDefGAMutationRate
        self.pCrossover = constants.CDefGACrossoverRate
        self.nElitismReplacement = constants.CDefGAElitismReplacement
        self.setPopulationSize(constants.CDefGAPopulationSize)
        self.minimax = constants.minimaxType["maximize"]
        self.elitism = True

        # The new population
        self.new_population = None

        # Adapters
        self.dbAdapter = None
        self.migrationAdapter = None

        self.time_init = None
        self.max_time = None
        self.interactiveMode = interactive
        self.interactiveGen = -1
        self.GPMode = False

        self.selector = FunctionSlot("Selector")
        self.stepCallback = FunctionSlot("Generation Step Callback")
        self.terminationCriteria = FunctionSlot("Termination Criteria")
        self.selector.set(constants.CDefGASelector)
        self.allSlots = (self.selector, self.stepCallback, self.terminationCriteria)

        self.internalParams = {}

        self.currentGeneration = 0

        # Kwargs for the different functions
        self.kwargs = dict()

        # GP Testing
        for classes in constants.CDefGPGenomes:
            if isinstance(self.internalPop.oneSelfGenome, classes):
                self.setGPMode(True)
                break

        # Debugging
        log.debug("A genetic algorithm engine was created with " + str(self.nGenerations) + " generations")

        # Path of the engine object file
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Loading the genetic algorithm engine from '" + path + "' ...")

        # Load the GA object from file
        ga = serialization.load(path)

        # Set the path of the GA file
        ga.path = path

        # Return the GA
        return ga

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Save to the current path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the genetic algorithm engine to '" + path + "' ...")

        # Set the new path as the current path and save
        self.path = path
        serialization.dump(self, path, protocol=2)

    # -----------------------------------------------------------------

    @property
    def evaluator_kwargs(self):

        """
        This function ...
        :return:
        """

        if "evaluator" in self.kwargs: return self.kwargs["evaluator"]
        else: return dict() # return empty dict

    # -----------------------------------------------------------------

    @property
    def initializator_kwargs(self):

        """
        This function ...
        :return:
        """

        # Get the pre-defined mutator kwargs or create a new dictionary
        if "initializator" in self.kwargs: kwargs = deepcopy(self.kwargs["initializator"])
        else: kwargs = dict()

        kwargs["ga_engine"] = self

        # Return the kwargs
        return kwargs

    # -----------------------------------------------------------------

    @property
    def mutator_kwargs(self):

        """
        This function ...
        :return:
        """

        # Get the pre-defined mutator kwargs or create a new dictionary
        if "mutator" in self.kwargs: kwargs = deepcopy(self.kwargs["mutator"])
        else: kwargs = dict()

        # pmut=self.pMutation, ga_engine=self
        kwargs["pmut"] = self.pMutation
        kwargs["ga_engine"] = self

        # Return the kwargs
        return kwargs

    # -----------------------------------------------------------------

    @property
    def crossover_kwargs(self):

        """
        This function ...
        :return:
        """

        if "crossover" in self.kwargs: return self.kwargs["crossover"]
        else: return dict()

    # -----------------------------------------------------------------

    @property
    def callback_kwargs(self):

        """
        This function ...
        :return:
        """

        if "callback" in self.kwargs: return self.kwargs["callback"]
        else: return dict() # return empty dict

    # -----------------------------------------------------------------

    def setGPMode(self, bool_value):

        """
        Sets the Genetic Programming mode of the GA Engine
        :param bool_value: True or False
        """

        self.GPMode = bool_value

    # -----------------------------------------------------------------

    def getGPMode(self):

        """ Get the Genetic Programming mode of the GA Engine
        :rtype: True or False
        """

        return self.GPMode

    # -----------------------------------------------------------------

    def __call__(self, *args, **kwargs):

        """ A method to implement a callable object
        Example: ga_engine(freq_stats=10)
        .. versionadded:: 0.6
           The callable method.
        """

        if kwargs.get("freq_stats", None):
            return self.evolve(kwargs.get("freq_stats"))
        else: return self.evolve()

    # -----------------------------------------------------------------

    def setParams(self, **args):

        """ Set the internal params
        Example: ga.setParams(gp_terminals=['x', 'y'])
        :param args: params to save
        ..versionaddd:: 0.6
           Added the *setParams* method.
        """

        self.internalParams.update(args)

    # -----------------------------------------------------------------

    def getParam(self, key, nvl=None):

        """ Gets an internal parameter
        Example: ga.getParam("gp_terminals")
           ['x', 'y']
        :param key: the key of param
        :param nvl: if the key doesn't exist, the nvl will be returned
        ..versionaddd:: 0.6
           Added the *getParam* method.
        """

        return self.internalParams.get(key, nvl)

    # -----------------------------------------------------------------

    def setInteractiveGeneration(self, generation):

        """ Sets the generation in which the GA must enter in the
        Interactive Mode
        :param generation: the generation number, use "-1" to disable
        .. versionadded::0.6
           The *setInteractiveGeneration* method.
        """

        if generation < -1:
            utils.raiseException("Generation must be >= -1", ValueError)
        self.interactiveGen = generation

    # -----------------------------------------------------------------

    def getInteractiveGeneration(self):

        """ returns the generation in which the GA must enter in the
        Interactive Mode
        :rtype: the generation number or -1 if not set
        .. versionadded::0.6
           The *getInteractiveGeneration* method.
        """

        return self.interactiveGen

    # -----------------------------------------------------------------

    def setElitismReplacement(self, numreplace):

        """ Set the number of best individuals to copy to the next generation on the elitism
        :param numreplace: the number of individuals
        .. versionadded:: 0.6
           The *setElitismReplacement* method.
        """

        if numreplace < 1:
            utils.raiseException("Replacement number must be >= 1", ValueError)
        self.nElitismReplacement = numreplace

    # -----------------------------------------------------------------

    def setInteractiveMode(self, flag=True):

        """ Enable/disable the interactive mode
        :param flag: True or False
        .. versionadded: 0.6
           The *setInteractiveMode* method.
        """

        if type(flag) != BooleanType:
            utils.raiseException("Interactive Mode option must be True or False", TypeError)
        self.interactiveMode = flag

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        The string representation of the GA Engine
        """

        minimax_type = constants.minimaxType.keys()[constants.minimaxType.values().index(self.minimax)]
        ret = "- GSimpleGA\n"
        ret += "\tGP Mode:\t\t %s\n" % self.getGPMode()
        ret += "\tPopulation Size:\t %d\n" % self.internalPop.popSize
        ret += "\tGenerations:\t\t %d\n" % self.nGenerations
        ret += "\tCurrent Generation:\t %d\n" % self.currentGeneration
        ret += "\tMutation Rate:\t\t %.2f\n" % self.pMutation
        ret += "\tCrossover Rate:\t\t %.2f\n" % self.pCrossover
        ret += "\tMinimax Type:\t\t %s\n" % minimax_type.capitalize()
        ret += "\tElitism:\t\t %s\n" % self.elitism
        ret += "\tElitism Replacement:\t %d\n" % self.nElitismReplacement
        ret += "\tDB Adapter:\t\t %s\n" % self.dbAdapter
        for slot in self.allSlots:
            ret += "\t" + slot.__repr__()
        ret += "\n"

        # Return the string
        return ret

    # -----------------------------------------------------------------

    def setMultiProcessing(self, flag=True, full_copy=False, max_processes=None):

        """
        Sets the flag to enable/disable the use of python multiprocessing module.
        Use this option when you have more than one core on your CPU and when your
        evaluation function is very slow.

        Pyevolve will automaticly check if your Python version has **multiprocessing**
        support and if you have more than one single CPU core. If you don't have support
        or have just only one core, Pyevolve will not use the **multiprocessing**
        feature.

        Pyevolve uses the **multiprocessing** to execute the evaluation function over
        the individuals, so the use of this feature will make sense if you have a
        truly slow evaluation function (which is commom in GAs).

        The parameter "full_copy" defines where the individual data should be copied back
        after the evaluation or not. This parameter is useful when you change the
        individual in the evaluation function.

        :param flag: True (default) or False
        :param full_copy: True or False (default)
        :param max_processes: None (default) or an integer value

        .. warning:: Use this option only when your evaluation function is slow, so you'll
                     get a good tradeoff between the process communication speed and the
                     parallel evaluation. The use of the **multiprocessing** doesn't means
                     always a better performance.

        .. note:: To enable the multiprocessing option, you **MUST** add the *__main__* check
                  on your application, otherwise, it will result in errors. See more on the
                  `Python Docs <http://docs.python.org/library/multiprocessing.html#multiprocessing-programming>`__
                  site.

        .. versionadded:: 0.6
           The `setMultiProcessing` method.
        """

        if type(flag) != BooleanType:
            utils.raiseException("Multiprocessing option must be True or False", TypeError)

        if type(full_copy) != BooleanType:
            utils.raiseException("Multiprocessing 'full_copy' option must be True or False", TypeError)

        self.internalPop.setMultiProcessing(flag, full_copy, max_processes)

    # -----------------------------------------------------------------

    def setMigrationAdapter(self, migration_adapter=None):

        """
        Sets the Migration Adapter
        .. versionadded:: 0.6
           The `setMigrationAdapter` method.
        """

        self.migrationAdapter = migration_adapter
        if self.migrationAdapter is not None:
            self.migrationAdapter.set_engine(self)

    # -----------------------------------------------------------------

    def setDBAdapter(self, dbadapter=None):

        """
        Sets the DB Adapter of the GA Engine
        :param dbadapter: one of the :mod:`DBAdapters` classes instance
        .. warning:: the use the of a DB Adapter can reduce the speed performance of the
                     Genetic Algorithm.
        """

        if (dbadapter is not None) and (not isinstance(dbadapter, DBBaseAdapter)):
            utils.raiseException("The DB Adapter must be a DBBaseAdapter subclass", TypeError)
        self.dbAdapter = dbadapter

    # -----------------------------------------------------------------

    def set_kwargs(self, target, kwargs):

        """
        This function ...
        :param target:
        :param kwargs:
        :return:
        """

        self.kwargs[target] = kwargs

    # -----------------------------------------------------------------

    def set_ranges(self, minima, maxima):

        """
        This function ...
        :param minima:
        :param maxima:
        :return:
        """

        # Set the ranges of the internal population (is this sufficient??)
        self.internalPop.set_ranges(minima, maxima)

    # -----------------------------------------------------------------

    def setPopulationSize(self, size):

        """
        Sets the population size, calls setPopulationSize() of GPopulation
        :param size: the population size
        .. note:: the population size must be >= 2
        """

        if size < 2:
            utils.raiseException("population size must be >= 2", ValueError)
        self.internalPop.setPopulationSize(size)

    # -----------------------------------------------------------------

    def setSortType(self, sort_type):

        """
        Sets the sort type, constants.sortType["raw"]/constants.sortType["scaled"]
        Example: ga_engine.setSortType(constants.sortType["scaled"])
        :param sort_type: the Sort Type
        """

        if sort_type not in constants.sortType.values():
            utils.raiseException("sort type must be a constants.sortType type", TypeError)
        self.internalPop.sortType = sort_type

    # -----------------------------------------------------------------

    def setMutationRate(self, rate):

        """
        Sets the mutation rate, between 0.0 and 1.0
        :param rate: the rate, between 0.0 and 1.0
        """

        if (rate > 1.0) or (rate < 0.0):
            utils.raiseException("Mutation rate must be >= 0.0 and <= 1.0", ValueError)
        self.pMutation = rate

    # -----------------------------------------------------------------

    def setCrossoverRate(self, rate):

        """
        Sets the crossover rate, between 0.0 and 1.0
        :param rate: the rate, between 0.0 and 1.0
        """

        if (rate > 1.0) or (rate < 0.0):
            utils.raiseException("Crossover rate must be >= 0.0 and <= 1.0", ValueError)
        self.pCrossover = rate

    # -----------------------------------------------------------------

    def setGenerations(self, num_gens):

        """
        Sets the number of generations to evolve
        :param num_gens: the number of generations
        """

        if num_gens < 1:
            utils.raiseException("Number of generations must be >= 1", ValueError)
        self.nGenerations = num_gens

    # -----------------------------------------------------------------

    def getGenerations(self):

        """
        Return the number of generations to evolve
        :rtype: the number of generations
        .. versionadded:: 0.6
           Added the *getGenerations* method
        """

        return self.nGenerations

    # -----------------------------------------------------------------

    def getMinimax(self):

        """ Gets the minimize/maximize mode
        :rtype: the constants.minimaxType type
        """

        return self.minimax

    # -----------------------------------------------------------------

    def setMinimax(self, mtype):

        """
        Sets the minimize/maximize mode, use constants.minimaxType
        :param mtype: the minimax mode, from constants.minimaxType
        """

        if mtype not in constants.minimaxType.values():
            utils.raiseException("Minimax must be maximize or minimize", TypeError)
        self.minimax = mtype

    # -----------------------------------------------------------------

    def getCurrentGeneration(self):

        """ Gets the current generation
        :rtype: the current generation
        """

        return self.currentGeneration

    # -----------------------------------------------------------------

    def setElitism(self, flag):

        """
        Sets the elitism option, True or False
        :param flag: True or False
        """

        if type(flag) != BooleanType:
            utils.raiseException("Elitism option must be True or False", TypeError)
        self.elitism = flag

    # -----------------------------------------------------------------

    def getDBAdapter(self):

        """ Gets the DB Adapter of the GA Engine
        :rtype: a instance from one of the :mod:`DBAdapters` classes
        """

        return self.dbAdapter

    # -----------------------------------------------------------------

    def setMaxTime(self, seconds):

        """ Sets the maximun evolve time of the GA Engine
        :param seconds: maximum time in seconds
        """

        self.max_time = seconds

    # -----------------------------------------------------------------

    def getMaxTime(self):

        """ Get the maximun evolve time of the GA Engine
        :rtype: True or False
        """

        return self.max_time

    # -----------------------------------------------------------------

    def bestIndividual(self):

        """ Returns the population best individual
        :rtype: the best individual
        """

        return self.internalPop.bestRaw()

    # -----------------------------------------------------------------

    def worstIndividual(self):

        """
        Returns the population worst individual
        :rtype: the best individual
        """

        return self.internalPop.worstRaw()

    # -----------------------------------------------------------------

    def __gp_catch_functions(self, prefix):

        """
        Internally used to catch functions with some specific prefix
        as non-terminals of the GP core
        """

        import __main__ as mod_main

        function_set = {}

        main_dict = mod_main.__dict__
        for obj, addr in main_dict.items():
            if obj[0:len(prefix)] == prefix:
                try:
                    op_len = addr.func_code.co_argcount
                except:
                    continue
                function_set[obj] = op_len

        if len(function_set) <= 0:
            utils.raiseException("No function set found using function prefix '%s' !" % prefix, ValueError)

        self.setParams(gp_function_set=function_set)

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function initializes the GA Engine. Create and initialize the first population
        """

        # Inform the user
        log.info("Initializing the GA engine ...")

        # Keep track of the time passed
        self.time_init = time()

        # Create the first population
        self.internalPop.create(minimax=self.minimax)

        # Initialize the population (initializes all individuals of the population)
        #self.internalPop.initialize(ga_engine=self)
        self.internalPop.initialize(**self.initializator_kwargs)

    # -----------------------------------------------------------------

    def set_scores(self, scores, check=None):

        """
        This function ...
        :param scores:
        :param check:
        :return:
        """

        # Set the scores for the initial population
        if self.is_initial_generation: self.set_scores_for_population(self.internalPop, scores, check)

        # Set the scores for the new population
        else:

            # Set scores
            self.set_scores_for_population(self.new_population, scores, check)

            # Replace
            if self.new_population is not None: self.replace_internal_population()

            # Increment the current generation number
            self.currentGeneration += 1

        # Sort the internal population
        self.internalPop.sort()

        # Set new pop to None
        self.new_population = None

    # -----------------------------------------------------------------

    def set_scores_for_population(self, population, scores, check=None):

        """
        This function ...
        :param population:
        :param scores:
        :param check:
        :return:
        """

        index = 0
        for individual in population:

            if check is not None:

                # Get the parameter values for this individual
                parameter_a = individual.genomeList[0]
                parameter_b = individual.genomeList[1]
                parameter_c = individual.genomeList[2]

                parameter_a_check = check[0][index]
                parameter_b_check = check[1][index]
                parameter_c_check = check[2][index]
                rel_diff_a = abs((parameter_a - parameter_a_check) / parameter_a)
                rel_diff_b = abs((parameter_b - parameter_b_check) / parameter_b)
                rel_diff_c = abs((parameter_c - parameter_c_check) / parameter_c)
                assert np.isclose(parameter_a, parameter_a_check, rtol=1e-11), rel_diff_a
                assert np.isclose(parameter_b, parameter_b_check, rtol=1e-11), rel_diff_b
                assert np.isclose(parameter_c, parameter_c_check, rtol=1e-11), rel_diff_c

            # Set the score
            individual.score = scores[index]

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def get_population(self):

        """ Return the internal population of GA Engine
        :rtype: the population (:class:`GPopulation.GPopulation`)
        """

        return self.internalPop

    # -----------------------------------------------------------------

    def getStatistics(self):

        """
        Gets the Statistics class instance of current generation
        :rtype: the statistics instance (:class:`Statistics.Statistics`)
        """

        return self.internalPop.getStatistics()

    # -----------------------------------------------------------------

    def generate_new_population(self, silent=False):

        """
        This function ...
        :return:
        """

        # Inform the user
        if not silent: log.info("Creating generation " + str(self.currentGeneration) + " ...")

        # Clone the current internal population
        new_population = Population(self.internalPop)
        log.debug("Population was cloned.")

        #size_iterate = len(self.internalPop)
        size_iterate = self.internalPop.popSize

        # Odd population size
        if size_iterate % 2 != 0: raise ValueError("The population size cannot be odd")
            #size_iterate -= 1

        crossover_empty = self.select(popID=self.currentGeneration).crossover.isEmpty()

        mutator_kwargs = self.mutator_kwargs

        # Loop over the population
        for i in xrange(0, size_iterate, 2):

            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)

            if not crossover_empty and self.pCrossover >= 1.0:
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2, **self.crossover_kwargs):
                    (sister, brother) = it
            else:
                if not crossover_empty and utils.randomFlipCoin(self.pCrossover):
                    for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2, **self.crossover_kwargs):
                        (sister, brother) = it
                else:
                    sister = genomeMom.clone()
                    brother = genomeDad.clone()

            #sister.mutate(pmut=self.pMutation, ga_engine=self)
            #brother.mutate(pmut=self.pMutation, ga_engine=self)
            sister.mutate(**mutator_kwargs)
            brother.mutate(**mutator_kwargs)

            new_population.internalPop.append(sister)
            new_population.internalPop.append(brother)

        if len(self.internalPop) % 2 != 0:

            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)

            if utils.randomFlipCoin(self.pCrossover):
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=1, **self.crossover_kwargs):
                    (sister, brother) = it
            else:
                sister = prng.choice([genomeMom, genomeDad])
                sister = sister.clone()
                #sister.mutate(pmut=self.pMutation, ga_engine=self)
                sister.mutate(**mutator_kwargs)

            new_population.internalPop.append(sister)

        # Set the new population
        self.new_population = new_population

    # -----------------------------------------------------------------

    def step(self, silent=False):

        """
        This function performs one step in the evolution, i.e. one generation
        :param silent:
        """

        # Inform the user
        if not silent: log.info("Performing step in the evolutionary process ...")

        # Generate the new population
        self.generate_new_population(silent=silent)

        # Evaluate
        #print(self.generation_description)
        self.new_population.evaluate(silent, **self.evaluator_kwargs)

        # Replace population
        self.replace_internal_population()

        # Sort the population
        self.internalPop.sort()

        # Set new population to None
        self.new_population = None

        # Inform the user
        #log.success("The generation %d was finished.", self.currentGeneration)

        # Increment the current generation number
        self.currentGeneration += 1

        if self.max_time:
           total_time = time() - self.time_init
           if total_time > self.max_time:
              return True

        # Return whether the desired number of generations has been reached
        return self.currentGeneration == self.nGenerations

    # -----------------------------------------------------------------

    def replace_internal_population(self):

        """
        This function ...
        :return:
        """

        # Elitism:
        if self.elitism: self.do_elitism(self.new_population)

        # Set the new population as the internal population and sort it
        self.internalPop = self.new_population

    # -----------------------------------------------------------------

    def do_elitism(self, new_population):

        """
        This function ...
        :param new_population:
        :return:
        """

        # Elitism concept has different meanings for metaheuristic and particular for GA. In general, elitism
        # is related with memory: "remember the best solution found" (kind of greedying). In the most traditional way,
        # for evolutionary algorithms (GA, EA, DE...), elitism implies the best solution found is used for to build the
        # next generation. Particularly, for GA it means keeping the best individual found intact for next generation.
        # In GA, there is something like a  proof for multiobjective optimization (villalobos, coello and hernandez, 2005)
        # that linked convergence and elitism which is a kind of cool stuff.
        #
        # Elitism involves copying a small proportion of the fittest candidates, unchanged, into the next generation.
        # This can sometimes have a dramatic impact on performance by ensuring that the EA does not waste time
        # re-discovering previously discarded partial solutions. Candidate solutions that are preserved unchanged
        # through elitism remain eligible for selection as parents when breeding the remainder of the next generation.

        log.debug("Doing elitism ...")

        if self.getMinimax() == constants.minimaxType["maximize"]:

            for i in xrange(self.nElitismReplacement):

                ##re-evaluate before being sure this is the best
                # self.internalPop.bestRaw(i).evaluate(**self.evaluator_kwargs) # IS THIS REALLY NECESSARY ?

                if self.internalPop.bestRaw(i).score > new_population.bestRaw(i).score:
                    new_population[len(new_population) - 1 - i] = self.internalPop.bestRaw(i)

        elif self.getMinimax() == constants.minimaxType["minimize"]:

            for i in xrange(self.nElitismReplacement):

                ##re-evaluate before being sure this is the best
                # self.internalPop.bestRaw(i).evaluate(**self.evaluator_kwargs) # IS THIS REALLY NECESSARY ?

                if self.internalPop.bestRaw(i).score < new_population.bestRaw(i).score:
                    new_population[len(new_population) - 1 - i] = self.internalPop.bestRaw(i)

    # -----------------------------------------------------------------

    def printStats(self):

        """
        Print generation statistics
        :rtype: the printed statistics as string
        .. versionchanged:: 0.6
           The return of *printStats* method.
        """

        percent = self.currentGeneration * 100 / float(self.nGenerations)
        message = "Gen. %d (%.2f%%):" % (self.currentGeneration, percent)
        log.info(message)
        print(message,)
        sys_stdout.flush()
        self.internalPop.statistics()
        stat_ret = self.internalPop.printStats()
        return message + stat_ret

    # -----------------------------------------------------------------

    def printTimeElapsed(self):

        """
        Shows the time elapsed since the begin of evolution
        """

        total_time = time() - self.time_init
        print("Total time elapsed: %.3f seconds." % total_time)
        return total_time

    # -----------------------------------------------------------------

    def dumpStatsDB(self):

        """
        Dumps the current statistics to database adapter
        """

        log.debug("Dumping stats to the DB Adapter")
        self.internalPop.statistics()
        self.dbAdapter.insert(self)

    # -----------------------------------------------------------------

    def evolve(self, freq_stats=0, progress_bar=False):

        """ Do all the generations until the termination criteria, accepts
        the freq_stats (default is 0) to dump statistics at n-generation

        Example: ga_engine.evolve(freq_stats=10)

        :param freq_stats: if greater than 0, the statistics will be
                           printed every freq_stats generation.
        :param progress_bar: show progress bar
        :rtype: returns the best individual of the evolution

        .. versionadded:: 0.6
           the return of the best individual

        """

        # 1. Initialize
        self.initialize_evolution()

        # 2. Do the evolution loop
        self.evolve_loop(freq_stats, progress_bar=progress_bar)

        # 3. Finish evolution, return best individual
        return self.finish_evolution()

    # -----------------------------------------------------------------

    @property
    def is_initial_generation(self):

        """
        This function ...
        :return:
        """

        if self.new_population is None:
            if self.currentGeneration > 0: raise RuntimeError("Inconsistent state: 'new_population' does exist but 'currentGeneration' >0")
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def generation_description(self):

        """
        This function ...
        :return:
        """

        if self.is_initial_generation: return "Initial generation"
        else: return "Generation " + str(self.currentGeneration)

    # -----------------------------------------------------------------

    def initialize_evolution(self):

        """
        This function ...
        :return:
        """

        self.initialize()

        # Inform the user ...
        log.info("Evaluating and sorting the initial population ...")

        # Evaluate and sort the internal population
        self.internalPop.evaluate(True, **self.evaluator_kwargs)
        self.sort_internal_population()

    # -----------------------------------------------------------------

    def sort_internal_population(self):

        """
        This function ...
        :return:
        """

        self.internalPop.sort()

    # -----------------------------------------------------------------

    def evolve_loop(self, freq_stats, progress_bar=False):

        """
        This function ...
        :param freq_stats:
        :param progress_bar:
        :return:
        """

        # Inform the user ...
        log.info("Starting the evolutionary algorithm ...")

        if progress_bar:

            from ...core.tools.progress import Bar, BAR_FILLED_CHAR, BAR_EMPTY_CHAR

            freq_stats = 10000
            silent = True

            counter = 1

            # Progress bar
            with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                     filled_char=BAR_FILLED_CHAR, expected_size=self.nGenerations, every=1, add_datetime=True) as bar:

                # Loop
                while True:

                    bar.show(counter)
                    if not self.evolve_generation(freq_stats, silent=silent): break
                    counter += 1
        else:

            while True:
                if not self.evolve_generation(freq_stats): break

    # -----------------------------------------------------------------

    def finish_evolution(self, silent=True):

        """
        This function ...
        :return:
        """

        if not silent:

            self.printStats()
            self.printTimeElapsed()

        if self.dbAdapter:
            log.debug("Closing the DB Adapter")
            if not (self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0):
                self.dumpStatsDB()
            self.dbAdapter.commitAndClose()

        if self.migrationAdapter:
            log.debug("Closing the Migration Adapter")
            self.migrationAdapter.stop()

        return self.bestIndividual()

    # -----------------------------------------------------------------

    def evolve_generation(self, freq_stats, silent=False):

        """
        This function ...
        :param freq_stats:
        :param silent:
        :return:
        """

        stopFlagCallback = False
        stopFlagTerminationCriteria = False

        if self.migrationAdapter:
            log.debug("Migration adapter: exchange")
            self.migrationAdapter.exchange()
            self.internalPop.clearFlags()
            self.internalPop.sort()

        if not self.stepCallback.isEmpty():
            if "callback" in self.kwargs:
                for it in self.stepCallback.applyFunctions(self, **self.kwargs["callback"]):
                    stopFlagCallback = it
            else:
                for it in self.stepCallback.applyFunctions(self):
                    stopFlagCallback = it

        if not self.terminationCriteria.isEmpty():
            for it in self.terminationCriteria.applyFunctions(self):
                stopFlagTerminationCriteria = it

        if freq_stats and not silent:
            if (self.currentGeneration % freq_stats == 0) or (self.getCurrentGeneration() == 0):
                self.printStats()

        if self.dbAdapter:
            if self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0:
                self.dumpStatsDB()

        if stopFlagTerminationCriteria:
            log.debug("Evolution stopped by the Termination Criteria !")
            if freq_stats and not silent:
                print("\n\tEvolution stopped by Termination Criteria function !\n")
            return False

        if stopFlagCallback:
            log.debug("Evolution stopped by Step Callback function !")
            if freq_stats and not silent:
                print("\n\tEvolution stopped by Step Callback function !\n")
            return False

        # Perform step (generation)
        if self.step(silent=silent): return False

        # Stopping criterion
        return True

    # -----------------------------------------------------------------

    def select(self, **args):

        """
        Select one individual from population
        :param args: this parameters will be sent to the selector
        """

        for it in self.selector.applyFunctions(self.internalPop, **args):
            return it

# -----------------------------------------------------------------
