#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.reference This module contains a function that
#  calls the original pyevolve implementation, for benchmarking.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from ..core.tools.logging import log
from ..core.basics.range import RealRange, IntegerRange

# -----------------------------------------------------------------

def call(settings, input_dict):

    """
    This function ...
    :param settings:
    :param input_dict:
    :return:
    """

    # Inform the user
    log.info("Calling Pyevolve ...")

    # Import pyevolve modules
    #from pyevolve import logEnable
    from pyevolve.G1DList import G1DList
    from pyevolve.GSimpleGA import GSimpleGA, RawScoreCriteria
    from pyevolve import Initializators
    from pyevolve import Mutators
    from pyevolve import Consts

    #logEnable(None, "DEBUG")

    # Get the settings
    output = settings["output"]
    nparameters = settings["nparameters"] if "nparameters" in settings else None
    nindividuals = settings["nindividuals"] if "nindividuals" in settings else None
    parameter_range = settings["parameter_range"] if "parameter_range" in settings else None
    best_raw_score = settings["best_raw_score"] if "best_raw_score" in settings else None
    round_decimal = settings["round_decimal"] if "round_decimal" in settings else None
    ngenerations = settings["ngenerations"] if "ngenerations" in settings else None
    mutation_rate = settings["mutation_rate"] if "mutation_rate" in settings else None
    crossover_rate = settings["crossover_rate"] if "crossover_rate" in settings else None
    stats_freq = settings["stats_freq"] if "stats_freq" in settings else None
    mutation_method = settings["mutation_method"] if "mutation_method" in settings else None
    min_or_max = settings["min_or_max"] if "min_or_max" in settings else None

    # No messages
    stats_freq = 0

    # Get the input
    evaluator = input_dict["evaluator"] # this is required!
    genome = input_dict["genome"].clone() if "genome" in input_dict else None
    initializator = input_dict["initializator"] if "initializator" in input_dict else None
    mutator = input_dict["mutator"] if "mutator" in input_dict else None
    crossover = input_dict["crossover"] if "crossover" in input_dict else None
    callback = input_dict["callback"] if "callback" in input_dict else None

    # Create genome instance
    if genome is not None: genome = G1DList(nparameters)
    if parameter_range is not None: genome.setParams(rangemin=parameter_range.min, rangemax=parameter_range.max)
    if best_raw_score is not None: genome.setParams(bestrawscore=best_raw_score)
    if round_decimal is not None: genome.setParams(rounddecimal=round_decimal)

    # Set initializator
    if initializator is not None: genome.initializator.set(initializator)
    else:
        if isinstance(parameter_range, IntegerRange): genome.initializator.set(Initializators.G1DListInitializatorInteger)
        elif isinstance(parameter_range, RealRange): genome.initializator.set(Initializators.G1DListInitializatorReal)
        else: raise ValueError("Invalid parameter range")

    # Set mutation method
    if mutator is not None: genome.mutator.set(mutator)
    else:
        if isinstance(parameter_range, IntegerRange):
            if mutation_method == "range": genome.mutator.set(Mutators.G1DListMutatorIntegerRange)
            elif mutation_method == "gaussian": genome.mutator.set(Mutators.G1DListMutatorIntegerGaussian)
            elif mutation_method == "binary": genome.mutator.set(Mutators.G1DListMutatorIntegerBinary)
            else: raise ValueError("Mutation method '" + mutation_method + "' not recognized")
        elif isinstance(parameter_range, RealRange):
            if mutation_method == "range": genome.mutator.set(Mutators.G1DListMutatorRealRange)
            elif mutation_method == "gaussian": genome.mutator.set(Mutators.G1DListMutatorRealGaussian)
            else: raise ValueError("Mutation method '" + mutation_method + "' not valid for genome of real values")
        else: raise ValueError("Invalid parameter range")

    # Set crossover
    if crossover is not None: genome.crossover.set(crossover)

    # The evaluator function (objective function)
    genome.evaluator.set(evaluator)

    # Genetic Algorithm Instance
    ga = GSimpleGA(genome)
    ga.terminationCriteria.set(RawScoreCriteria)
    ga.setMinimax(Consts.minimaxType[min_or_max])
    ga.setGenerations(ngenerations)
    if crossover_rate is not None: ga.setCrossoverRate(crossover_rate)
    ga.setPopulationSize(nindividuals)
    ga.setMutationRate(mutation_rate)

    # Set the callback function
    if callback is not None: ga.stepCallback.set(callback)

    # -----------------------------------------------------------------

    # Do the evolution
    ga.evolve(freq_stats=stats_freq)

    # -----------------------------------------------------------------

    # Get the best individual
    best = ga.bestIndividual()

    # Return the best individual
    return best

# -----------------------------------------------------------------
