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
from time import time
from types import BooleanType
from sys import stdout as sys_stdout

# Import other evolve modules
from pts.evolve.core.population import Population, NamedPopulation, PopulationBase
from pts.evolve.core.functionslot import FunctionSlot
from pts.evolve.core.genome import GenomeBase
from pts.evolve.core.adapters import DataBaseAdapter
import pts.evolve.core.constants as constants
import pts.evolve.core.utils as utils

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import serialization
from ...core.tools.random import prng
from ...core.tools.stringify import tostr
from ...core.basics.containers import DefaultOrderedDict

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

    #if ga_engine.getMinimax() == constants.minimaxType["maximize"]:
    if ga_engine.getMinimax() == "maximize":

        if roundDecimal is not None:
            return round(bestRawScore, roundDecimal) <= round(ind.score, roundDecimal)
        else:
            return bestRawScore <= ind.score

    #else:
    elif ga_engine.getMinimax() == "minimize":

        if roundDecimal is not None:
            return round(bestRawScore, roundDecimal) >= round(ind.score, roundDecimal)
        else:
            return bestRawScore >= ind.score

    else: raise RuntimeError("Invalid minimax type: " + str(ga_engine.getMinimax()))

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

    def __init__(self, genome_or_pop, interactive=True, named_individuals=False):

        """
        Initializator of GSimpleGA
        :param genome_or_pop:
        :param interactive:
        """

        if type(interactive) != BooleanType:
            utils.raiseException("Interactive Mode option must be True or False", TypeError)

        # Genome is passed
        if isinstance(genome_or_pop, GenomeBase):

            # Set named flag
            self.named_individuals = named_individuals

            # Create the internal population
            if self.named_individuals: self.internalPop = NamedPopulation(genome_or_pop)
            else: self.internalPop = Population(genome_or_pop)

            # Set the flag that we still have to do the population initialization
            self.initialized_population = False

            # Initialize
            self.setPopulationSize(constants.CDefGAPopulationSize)

        # Initial population is passed
        elif isinstance(genome_or_pop, PopulationBase):

            # Check whether the population is valid
            if genome_or_pop.popSize < 2: raise ValueError("The population size cannot be smaller than 2")
            if genome_or_pop.popSize % 2 != 0: raise ValueError("The population size cannot be odd")

            # Checks passed: set the internal (initial) population
            self.internalPop = genome_or_pop

            # Check whether the population uses named or unnamed individuals
            if isinstance(genome_or_pop, NamedPopulation): self.named_individuals = True
            elif isinstance(genome_or_pop, Population): self.named_individuals = False
            else: raise ValueError("Invalid population object")

            # Set the flag that we don't have to do the population initialization anymore
            self.initialized_population = True

        # Unknown input
        else: raise ValueError("Genome or population should be passed")

        # Initialize things
        self.nGenerations = constants.CDefGAGenerations
        self.pMutation = constants.CDefGAMutationRate
        self.pCrossover = constants.CDefGACrossoverRate
        self.nElitismReplacement = constants.CDefGAElitismReplacement
        self.minimax = constants.CDefPopMinimax
        self.elitism = True

        # The new population
        self.new_population = None

        # Adapters
        self.database_adapters = []
        self.migrationAdapter = None

        # Properties
        self.time_init = None
        self.max_time = None
        self.interactiveMode = interactive
        self.interactiveGen = -1
        self.GPMode = False

        # Function slots
        self.selector = FunctionSlot("Selector")
        self.stepCallback = FunctionSlot("Generation Step Callback")
        self.terminationCriteria = FunctionSlot("Termination Criteria")
        self.selector.set(constants.CDefGASelector)
        self.allSlots = (self.selector, self.stepCallback, self.terminationCriteria)

        # The internal parameters
        self.internalParams = {}

        # A counter that keeps track of which generation we are running
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

        # Remove db adapter
        #self.dbAdapter = None

        # Remove the adapters
        self.database_adapters = []

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

    def get_parameters(self):

        """
        This function ...
        :return: 
        """

        return self.internalParams.copy()

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

        #minimax_type = constants.minimaxType.keys()[constants.minimaxType.values().index(self.minimax)]
        ret = "- GSimpleGA\n"
        ret += "\tGP Mode:\t\t %s\n" % self.getGPMode()
        ret += "\tPopulation Size:\t %d\n" % self.internalPop.popSize
        ret += "\tGenerations:\t\t %d\n" % self.nGenerations
        ret += "\tCurrent Generation:\t %d\n" % self.currentGeneration
        ret += "\tMutation Rate:\t\t %.2f\n" % self.pMutation
        ret += "\tCrossover Rate:\t\t %.2f\n" % self.pCrossover
        ret += "\tMinimax Type:\t\t %s\n" % self.minimax.capitalize()
        ret += "\tElitism:\t\t %s\n" % self.elitism
        ret += "\tElitism Replacement:\t %d\n" % self.nElitismReplacement
        #ret += "\tDB Adapter:\t\t %s\n" % self.dbAdapter
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

    def add_database_adapter(self, adapter):

        """
        This function ...
        :param adapter:
        :return:
        """

        # Check
        if not isinstance(adapter, DataBaseAdapter): raise ValueError("The adapter must be a DataBaseAdapter subclass")

        # Add an adapter
        self.database_adapters.append(adapter)

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

        if self.initialized_population: self.check_population_size(size)
        else:

            # Check the value
            if size < 2: raise ValueError("Population size must be >= 2") #utils.raiseException("population size must be >= 2", ValueError)

            # Other check
            if size % 2 != 0: raise ValueError("The population size cannot be odd")

            # Set the value
            self.internalPop.setPopulationSize(size)

    # -----------------------------------------------------------------

    def check_population_size(self, size):

        """
        This function ...
        :param size: 
        :return: 
        """

        # Check the value
        if size < 2: raise ValueError("Population size must be >= 2")

        # Other check
        if size % 2 != 0: raise ValueError("The population size cannot be odd")

        # Check
        self.internalPop.check_population_size(size)

    # -----------------------------------------------------------------

    def set_or_check_population_size(self, size):

        """
        This function ...
        :param size: 
        :return: 
        """

        # Check the value
        if size < 2: raise ValueError("Population size must be >= 2")

        # Other check
        if size % 2 != 0: raise ValueError("The population size cannot be odd")

        # Set or check
        self.internalPop.set_or_check_population_size(size)

    # -----------------------------------------------------------------

    def setSortType(self, sort_type):

        """
        Sets the sort type, constants.sortType["raw"]/constants.sortType["scaled"]
        Example: ga_engine.setSortType(constants.sortType["scaled"])
        :param sort_type: the Sort Type
        """

        # Check the value
        if sort_type not in constants.sortType.values(): raise ValueError("Sort type must be a constants.sortType type") #utils.raiseException("sort type must be a constants.sortType type", TypeError)

        # Set the value
        self.internalPop.sortType = sort_type

    # -----------------------------------------------------------------

    def setMutationRate(self, rate):

        """
        Sets the mutation rate, between 0.0 and 1.0
        :param rate: the rate, between 0.0 and 1.0
        """

        if (rate > 1.0) or (rate < 0.0): utils.raiseException("Mutation rate must be >= 0.0 and <= 1.0", ValueError)
        self.pMutation = rate

    # -----------------------------------------------------------------

    def setCrossoverRate(self, rate):

        """
        Sets the crossover rate, between 0.0 and 1.0
        :param rate: the rate, between 0.0 and 1.0
        """

        if (rate > 1.0) or (rate < 0.0): utils.raiseException("Crossover rate must be >= 0.0 and <= 1.0", ValueError)
        self.pCrossover = rate

    # -----------------------------------------------------------------

    def setGenerations(self, num_gens):

        """
        Sets the number of generations to evolve
        :param num_gens: the number of generations
        """

        if num_gens < 1: utils.raiseException("Number of generations must be >= 1", ValueError)
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

    @property
    def minimize(self):

        """
        This function ...
        :return:
        """

        return self.minimax == "minimize"

    # -----------------------------------------------------------------

    @property
    def maximize(self):

        """
        This function ...
        :return:
        """

        return self.minimax == "maximize"

    # -----------------------------------------------------------------

    def setMinimax(self, mtype):

        """
        Sets the minimize/maximize mode, use constants.minimaxType
        :param mtype: the minimax mode, from constants.minimaxType
        """

        #if mtype not in constants.minimaxType.values():
        #    utils.raiseException("Minimax must be maximize or minimize", TypeError)

        if mtype not in ["minimize", "maximize"]: raise ValueError("Minimax must be 'maximize' or 'minimize'")

        # Set
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

    def setMaxTime(self, seconds):

        """
        Sets the maximun evolve time of the GA Engine
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

    #def __gp_catch_functions(self, prefix):

        #"""
        #Internally used to catch functions with some specific prefix
        #as non-terminals of the GP core
        #"""

        #import __main__ as mod_main

        #function_set = {}

        #main_dict = mod_main.__dict__
        #for obj, addr in main_dict.items():
        #    if obj[0:len(prefix)] == prefix:
        #        try:
        #            op_len = addr.func_code.co_argcount
        #        except:
        #            continue
        #        function_set[obj] = op_len
        #
        #if len(function_set) <= 0:
        #    utils.raiseException("No function set found using function prefix '%s' !" % prefix, ValueError)
        #
        #self.setParams(gp_function_set=function_set)

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function initializes the GA Engine. Create and initialize the first population
        """

        # Inform the user
        log.info("Initializing the GA engine ...")

        # Keep track of the time passed
        self.time_init = time()

        # Initialize population
        if not self.initialized_population: self.initialize_population()

    # -----------------------------------------------------------------

    def initialize_population(self):

        """
        This function initilizes the population
        :return: 
        """

        # Create the first population
        self.internalPop.create(minimax=self.minimax)

        # Initialize the population (initializes all individuals of the population)
        #self.internalPop.initialize(ga_engine=self)
        self.internalPop.initialize(**self.initializator_kwargs)

    # -----------------------------------------------------------------

    def set_scores(self, scores, names=None, check=None, rtol=1e-5, atol=1e-8, binary_parameters=None):

        """
        This function ...
        :param scores:
        :param names:
        :param check:
        :param rtol:
        :param atol:
        :param binary_parameters:
        :return:
        """

        # Elitism data
        elitism_data = None

        # Set the scores for the initial population
        if self.is_initial_generation: self.set_scores_for_population(self.internalPop, scores, names, check, rtol=rtol, atol=atol, binary_parameters=binary_parameters)

        # Set the scores for the new population
        else:

            # Set scores
            self.set_scores_for_population(self.new_population, scores, names, check, rtol=rtol, atol=atol, binary_parameters=binary_parameters)

            # Replace
            if self.new_population is not None: elitism_data = self.replace_internal_population()

            # Increment the current generation number
            self.currentGeneration += 1

        # Sort the internal population
        self.internalPop.sort()

        # Set new pop to None
        self.new_population = None

        # Return the information about the elitism replacements as a result of the scoring
        return elitism_data

    # -----------------------------------------------------------------

    @property
    def scores(self):

        """
        This fucntion ...
        :return:
        """

        the_scores = dict()

        # Loop over the individual keys
        for key in self.population.keys:

            # Get the score
            score = self.population[key].score

            # Add to the dictionary
            the_scores[key] = score

        # Return the dictinoary
        return the_scores

    # -----------------------------------------------------------------

    def set_scores_for_population(self, population, scores, names=None, check=None, rtol=1e-5, atol=1e-8, binary_parameters=None):

        """
        This function ...
        :param population:
        :param scores:
        :param names:
        :param check:
        :param rtol:
        :param atol:
        :param binary_parameters:
        :return:
        """

        index = 0
        for individual_key in population.keys:

            # Get individual
            individual = population[individual_key]

            # Check based on name
            if names is not None:

                # Check whether the individual key corresponds to the name
                if names[index] != individual_key: raise ValueError("Check failed: individual has name '" + individual_key + "' but name check is '" + names[index] + "'")

            # Check based on genome, only if no name checks were available
            elif check is not None:

                # Import functions
                from pts.evolve.optimize.parameters import equal_genomes, show_genome_differences

                genome_a = individual.genomeList
                genome_b = check[index]

                # Check
                if not equal_genomes(genome_a, genome_b, rtol=rtol, atol=atol, binary_parameters=binary_parameters):

                    # Show error message
                    log.error("Failed check for individual " + str(index) + ":")

                    # Show the differences between the genomes
                    show_genome_differences(genome_a, genome_b, rtol=rtol, atol=atol, binary_parameters=binary_parameters)

                    # Raise error since check failed
                    raise ValueError("Check failed: individual = " + str(genome_a) + " and check = " + str(genome_b))

            # No checks
            else: pass # no checks have to be done

            # Set the score
            individual.score = scores[index]

            # Increment the index of the individuals
            index += 1

    # -----------------------------------------------------------------

    @property
    def population(self):

        """
        This function ...
        :return:
        """

        return self.internalPop

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

        # NEW
        self.dump_statistics_adapters()

        # Inform the user
        if not silent: log.info("Creating generation " + str(self.currentGeneration) + " ...")

        # Clone the current internal population
        new_population = self.internalPop.clone_population()
        log.debug("Population was cloned")

        size_iterate = self.internalPop.popSize

        # Odd population size
        if size_iterate % 2 != 0: raise ValueError("The population size cannot be odd")

        # Check whether a crossover function is set
        crossover_empty = self.select(popID=self.currentGeneration).crossover.isEmpty()

        # Set the mutator arguments
        mutator_kwargs = self.mutator_kwargs

        # Make a list to contain the crossover data
        crossover_data = []

        # Determine the generation index
        current_generation_index = self.currentGeneration
        new_generation_index = self.currentGeneration + 1

        # Loop over the population
        for i in xrange(0, size_iterate, 2):

            # Perform crossover
            mother_key, father_key, sister, brother, applied, details = self.perform_crossover(crossover_empty=crossover_empty)

            # Mutate
            sister.mutate(**mutator_kwargs)
            brother.mutate(**mutator_kwargs)

            # Add the sister and brother to the population
            sister_key = new_population.append(sister)
            brother_key = new_population.append(brother)

            # Add entry to the crossover data
            entry = [current_generation_index, mother_key, father_key, sister_key, brother_key, applied, details]
            crossover_data.append(entry)

        # Generate one more individual if the population size is not even
        if len(self.internalPop) % 2 != 0:

            # Perform crossover
            mother_key, father_key, sister, brother, applied, details = self.perform_crossover(crossover_empty=crossover_empty)

            # CHoose one child
            child = prng.choice([sister, brother])

            # NEW: mutate here
            child.mutate(**mutator_kwargs)

            # Add the sister to the new population
            sister_key = new_population.append(child)
            brother_key = None

            # Add entry to the crossover data
            entry = [current_generation_index, mother_key, father_key, sister_key, brother_key, applied, details]
            crossover_data.append(entry)

        # Set the new population
        self.new_population = new_population

        # Return the crossover data
        return crossover_data

    # -----------------------------------------------------------------

    def step(self, silent=False):

        """
        This function performs one step in the evolution, i.e. one generation
        :param silent:
        """

        # Inform the user
        if not silent: log.info("Performing step in the evolutionary process ...")

        # Debugging
        log.debug("Generating new population ...")

        # Generate the new population
        self.generate_new_population(silent=silent)

        # Debugging
        log.debug("Evaluating generation '" + self.generation_description + "' ...")

        # Evaluate
        self.new_population.evaluate(silent, **self.evaluator_kwargs)

        # Replace population
        elitism_data = self.replace_internal_population()

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

        # Elitism
        if self.elitism: elitism_data = self.do_elitism(self.new_population)
        else: elitism_data = None

        # Set the new population as the internal population and sort it
        self.internalPop = self.new_population

        # Return the elitism data
        return elitism_data

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

        # Inform the user
        log.info("Performing elitism ...")

        # Elitism data
        data = DefaultOrderedDict(list)

        # Determine the generation index
        new_generation_index = self.currentGeneration + 1

        # Debugging
        log.debug("Elitism is performed on the new population that will afterwards become generation " + str(new_generation_index-1) + " ...")

        # Loop over the number of elitism replacements
        for i in xrange(self.nElitismReplacement):

            ##re-evaluate before being sure this is the best
            # self.internalPop.bestRaw(i).evaluate(**self.evaluator_kwargs) # IS THIS REALLY NECESSARY ?

            # Get best individual (raw score) of old population
            old_best = self.internalPop.bestRaw(i)
            old_key = self.internalPop.keys[i]

            # Check if sorted
            self.internalPop.check_sorted()

            # Get best individual (raw score) of new population
            new_best = new_population.bestRaw(i)
            #new_key = new_population.keys[i]

            # Check if sorted
            new_population.check_sorted()

            # Determine ID of the old individual
            if isinstance(self.internalPop, NamedPopulation): old_id = old_key
            else: old_id = str(old_key) # index

            # Get scores
            old_best_raw = old_best.score
            new_best_raw = new_best.score

            # Get fitness
            old_best_fitness = old_best.fitness
            new_best_fitness = new_best.fitness

            # NEW: VERIFY THAT HE BEST RAW SCORES ARE INDEED THE MINIMUM OR MAXIMUM OF ALL THE SCORES IN THEIR GENERATION (DEBUGGING)
            if self.nElitismReplacement == 1:

                if self.minimize:
                    if old_best_raw != min(self.internalPop.scores):
                        log.error("Something went wrong picking the best individual of the old (internal) population: score of " + str(old_best_raw) + " but lowest score is " + str(min(self.internalPop.scores)))
                        log.error("All scores: " + tostr(self.internalPop.scores))
                        log.error("Fitness: " + str(old_best_fitness) + " and lowest fitness is " + str(min(self.internalPop.fitnesses)))
                        log.error("All fitnesses: " + tostr(self.internalPop.fitnesses))
                        exit()
                    if new_best_raw != min(new_population.scores):
                        log.error("Something went wrong picking the best individual of the new population: score of " + str(new_best_raw) + " but lowest score is " + str(min(new_population.scores)))
                        log.error("All scores: " + tostr(new_population.scores))
                        log.error("Fitness: " + str(new_best_fitness) + " and lowest fitness is " + str(min(new_population.fitnesses)))
                        log.error("All fitnesses: " + tostr(new_population.fitnesses))
                        exit()
                elif self.maximize:
                    if old_best_raw != max(self.internalPop.scores):
                        log.error("Something went wrong picking the best individual of the old (internal) population: score of " + str(old_best_raw) + " but highest score is " + str(max(self.internalPop.scores)))
                        log.error("All scores: " + tostr(self.internalPop.scores))
                        log.error("Fitness: " + str(old_best_fitness) + " and highest fitness is " + str(max(self.internalPop.fitnesses)))
                        log.error("All fitnesses: " + tostr(self.internalPop.fitnesses))
                        exit()
                    if new_best_raw != max(new_population.scores):
                        log.error("Something went wrong picking the best individual of the new population: score of " + str(new_best_raw) + " but highest score is " + str(max(new_population.score)))
                        log.error("All scores: " + tostr(new_population.scores))
                        log.error("Fitness: " + str(new_best_fitness) + " and highest fitness is " + str(max(new_population.fitnesses)))
                        log.error("All fitnesses: " + tostr(new_population.fitnesses))
                        exit()
                else: raise ValueError("Invalid state of 'minimax': must be 'maximize' or 'minimize")

            # Check condition, depending on the min max type
            if self.minimize: condition = old_best_raw < new_best_raw # old_best_raw > new_best_raw WAS WRONG?!
            elif self.maximize: condition = old_best_raw > new_best_raw #old_best_raw < new_best_raw WAS WRONG?!
            else: raise ValueError("Invalid state of 'minimax': must be 'maximize' or 'minimize'")

            # Determine the index of the individual to be replaced
            replacement_index = len(new_population) - 1 - i

            # NEW: VERIFY THAT THE RAW SCORE OF THE REPLACEMENT INDIVIDUAL IS INDEED THE WORST
            if self.nElitismReplacement == 1:

                replaced = new_population[replacement_index]
                if self.minimize:
                    if replaced.score != max(new_population.scores):
                        log.error("Something went wrong picking the worst individual of the new population: score of " + str(replaced.score) + " but highest score is " + str(max(new_population.scores)))
                        log.error("Fitness: " + str(replaced.fitness) + " and highest fitness is " + str(max(new_population.fitnesses)))
                        exit()
                elif self.maximize:
                    if replaced.score != min(new_population.scores):
                        log.error("Something went wrong picking the worst individual of the new population: score of " + str(replaced.score) + " but lowest score is " + str(min(new_population.scores)))
                        log.error("Fitness: " + str(replaced.fitness) + " and lowest fitness is " + str(min(new_population.fitnesses)))
                        exit()

            # Determine the individual ID
            if isinstance(new_population, NamedPopulation): individual_id = new_population.names[replacement_index]
            else: individual_id = str(replacement_index)

            # Replace the individual, if the condition is met
            if condition:

                # Determine the score of the individual to be replaced
                replaced = new_population[replacement_index]
                replaced_raw = replaced.score
                replaced_fitness = replaced.fitness

                # Debugging
                log.debug("Replacing the " + str(replacement_index) + "th individual [" + individual_id + "] (raw score=" + str(replaced_raw) + ", fitness=" + str(replaced_fitness) + ") with the " + str(i) + "th individual from the parent population [" + old_id + "] (raw score=" + str(old_best_raw) + ", fitness=" + str(old_best_fitness) + ")...")

                # Replace
                new_population.replace(replacement_index, old_key, old_best) # Replace so the individual ID is changed, but the order is preserved

            # No replacement
            else: replaced_raw = replaced_fitness = None

            # Set the data
            data["Generation"].append(new_generation_index)
            data["Elitism replacement"].append(i)
            data["Min_or_max"].append(self.getMinimax())
            data["Old best raw score"].append(old_best_raw)
            data["Old best fitness"].append(old_best_fitness)
            data["Old individual ID"].append(old_id)
            data["New best raw score"].append(new_best_raw)
            data["New best fitness"].append(new_best_fitness)
            data["Individual ID"].append(individual_id)
            data["Elitism performed"].append(condition)
            data["Replaced raw score"].append(replaced_raw)
            data["Replaced fitness"].append(replaced_fitness)

        # RECALCULATE THE FITNESSES AND/OR RESORT?????!!
        new_population.scale()
        new_population.sort()

        # Return the elitism data
        return data

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

    def dump_statistics_adapters(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Dumping data to the database adapters ...")

        # Calculate the statistics of the internal populations
        self.internalPop.statistics()

        # Loop over the adapters, add the state of the engine
        for adapter in self.database_adapters:

            # Check the desired frequency
            if self.currentGeneration % adapter.getStatsGenFreq() == 0: adapter.insert(self)

    # -----------------------------------------------------------------

    def evolve(self, freq_stats=0, progress_bar=False):

        """
        Do all the generations until the termination criteria, accepts
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
            if self.currentGeneration > 0: raise RuntimeError("Inconsistent state: 'new_population' does exist but 'currentGeneration' > 0")
            return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def is_first_generation(self):

        """
        This function ...
        :return: 
        """

        if self.currentGeneration == 1:
            if self.new_population is None: raise RuntimeError("Inconsistent state: current generation == 1 but 'new_population' does not exist")
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

        # Inform the user
        log.info("Finishing the evolution ...")

        if not silent:
            self.printStats()
            self.printTimeElapsed()

        # Close the DB adapter ==> NO, NOW THIS IS THE RESPONSIBILITY OF THE MODULE THAT CREATES THE DBADAPTER
        #if self.dbAdapter:
        #    if not (self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0): self.dumpStatsDB()
        #    #self.dbAdapter.commit_and_close()

        # Finish the database adapters
        #for adapter in self.database_adapters:
        #    if not (self.currentGeneration % adapter.getStatsGenFreq() == 0): adapter.insert(self)

        self.dump_statistics_adapters()

        # Stop the migration adapter ==> NO, NOW THIS IS THE RESPONSIBILITY OF THE MODULE THAT CREATES IT
        #if self.migrationAdapter:
        #    log.debug("Closing the Migration Adapter")
        #    self.migrationAdapter.stop()

        # Return the best individual
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

    def perform_crossover(self, crossover_empty=False):

        """
        This function ...
        :return:
        """

        applied = False
        details = None

        # Select mother and father
        #genomeMom = self.select(popID=self.currentGeneration)
        #genomeDad = self.select(popID=self.currentGeneration)

        mother_key = self.select(popID=self.currentGeneration, return_key=True)
        father_key = self.select(popID=self.currentGeneration, return_key=True)

        # Get mother and father genome
        genomeMom = self.internalPop[mother_key]
        genomeDad = self.internalPop[father_key]

        # Always crossover
        if not crossover_empty and self.pCrossover >= 1.0:

            # Apply
            # DIDN'T MAKE SENSE: SISTER AND BROTHER WERE JUST OVERWRITTEN FOR EVERY NEW CROSSOVER FUNCTION THAT WAS CALLED!
            #for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2, **self.crossover_kwargs):
            #    (sister, brother) = it

            # TODO: what to do with multiple crossover functions now?

            # Apply the crossover
            sister, brother, details = genomeMom.crossover.apply(0, mom=genomeMom, dad=genomeDad, count=2, return_details=True, **self.crossover_kwargs)

            # Set flag
            applied = True

        # Crossover with some probability
        else:

            # Apply crossover function(s)
            if not crossover_empty and utils.randomFlipCoin(self.pCrossover):

                # Apply
                # DIDN'T MAKE SENSE: SISTER AND BROTHER WERE JUST OVERWRITTEN FOR EVERY NEW CROSSOVER FUNCTION THAT WAS CALLED!
                #for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2, **self.crossover_kwargs):
                #    (sister, brother) = it

                # TODO: what to do with multiple crossover functions now?

                # Apply the crossover
                sister, brother, details = genomeMom.crossover.apply(0, mom=genomeMom, dad=genomeDad, count=2, return_details=True, **self.crossover_kwargs)

                # Set flag
                applied = True

            # Make clones of the mother and father
            else:

                sister = genomeMom.clone()
                brother = genomeDad.clone()

        # Return
        return mother_key, father_key, sister, brother, applied, details

# -----------------------------------------------------------------
