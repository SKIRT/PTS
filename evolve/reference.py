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

def call(settings, evaluator):

    """
    This function ...
    :param settings:
    :param evaluator:
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
    nparameters = settings["nparameters"]
    nindividuals = settings["nindividuals"]
    parameter_range = settings["parameter_range"]
    best_raw_score = settings["best_raw_score"]
    round_decimal = settings["round_decimal"]
    ngenerations = settings["ngenerations"]
    mutation_rate = settings["mutation_rate"]
    crossover_rate = settings["crossover_rate"]
    stats_freq = settings["stats_freq"]
    mutation_method = settings["mutation_method"]
    min_or_max = settings["min_or_max"]

    # No messages
    stats_freq = 0

    # Create genome instance
    genome = G1DList(nparameters)
    genome.setParams(rangemin=parameter_range.min, rangemax=parameter_range.max)
    if best_raw_score is not None: genome.setParams(bestrawscore=best_raw_score)
    if round_decimal is not None: genome.setParams(rounddecimal=round_decimal)

    # Set initializator
    if isinstance(parameter_range, IntegerRange): genome.initializator.set(Initializators.G1DListInitializatorInteger)
    elif isinstance(parameter_range, RealRange): genome.initializator.set(Initializators.G1DListInitializatorReal)
    else: raise ValueError("Invalid parameter range")

    # Set mutation method
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

    # Do the evolution
    ga.evolve(freq_stats=stats_freq)

    # -----------------------------------------------------------------

    # Get the best individual and its score
    best = ga.bestIndividual()
    score = best.getRawScore()
    #print("best", best)
    #print("score", score)

    # Return the best individual
    return best

# -----------------------------------------------------------------
