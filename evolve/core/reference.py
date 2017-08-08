#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.reference This module contains a function that
#  calls the original pyevolve implementation, for benchmarking.

# -----------------------------------------------------------------

# Import standard modules
from functools import partial

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.range import RealRange, IntegerRange
from ...core.tools.progress import Bar, BAR_FILLED_CHAR, BAR_EMPTY_CHAR

# -----------------------------------------------------------------

def print_generation(ga_engine):

    """
    This function ...
    :param ga_engine:
    :return:
    """

    #percent = self.currentGeneration * 100 / float(self.nGenerations)
    #message = "Gen. %d (%.2f%%):" % (self.currentGeneration, percent)

    print("Generation " + str(ga_engine.currentGeneration) + "/" + str(ga_engine.nGenerations))

# -----------------------------------------------------------------

def refresh_progress_bar(ga_engine, **kwargs):

    """
    This function ...
    :param ga_engine:
    :return:
    """

    bar = kwargs.pop("bar")

    counter = ga_engine.currentGeneration
    bar.show(counter)

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
    from pyevolve.G1DList import G1DList
    from pyevolve.G2DList import G2DList
    from pyevolve.G1DBinaryString import G1DBinaryString
    from pyevolve.G2DBinaryString import G2DBinaryString
    from pyevolve.GSimpleGA import GSimpleGA, RawScoreCriteria
    from pyevolve.Initializators import G1DListInitializatorInteger, G1DListInitializatorReal
    from pyevolve.Mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary, G1DListMutatorRealRange, G1DListMutatorRealGaussian
    from pyevolve.Consts import minimaxType

    # Get the settings
    output = settings["output"]
    genome_type = settings["genome_type"] if "genome_type" in settings else "list"
    genome_dimension = settings["genome_dimension"] if "genome_dimension" in settings else 1
    nparameters = settings["nparameters"] if "nparameters" in settings else None
    nparameters2 = settings["nparameters2"] if "nparameters2" in settings else None
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

    # Get kwargs for the functions
    evaluator_kwargs = input_dict["evaluator_kwargs"] if "evaluator_kwargs" in input_dict else None
    initializator_kwargs = input_dict["initializator_kwargs"] if "initializator_kwargs" in input_dict else None
    mutator_kwargs = input_dict["mutator_kwargs"] if "mutator_kwargs" in input_dict else None
    crossover_kwargs = input_dict["crossover_kwargs"] if "crossover_kwargs" in input_dict else None
    callback_kwargs = input_dict["callback_kwargs"] if "callback_kwargs" in input_dict else None

    # Create partial functions with kwargs set
    if evaluator_kwargs is not None: evaluator = partial(evaluator, **evaluator_kwargs)
    if initializator is not None and initializator_kwargs is not None: initializator = partial(initializator, **initializator_kwargs)
    if mutator is not None and mutator_kwargs is not None: mutator = partial(mutator, **mutator_kwargs)
    if crossover is not None and crossover_kwargs is not None: crossover = partial(crossover, **crossover_kwargs)
    if callback is not None and callback_kwargs is not None: callback = partial(callback, **callback_kwargs)

    # Create genome instance
    if genome is None:

        # List genome
        if genome_type == "list":

            # Create 1D genome
            if genome_dimension == 1: genome = G1DList(nparameters)
            elif genome_dimension == 2: genome = G2DList(nparameters, nparameters2)
            else: raise ValueError("Dimensions > 2 are not supported")

        # Binary string genome
        elif genome_type == "binary_string":

            # 1D or 2D
            if genome_dimension == 1: genome = G1DBinaryString(nparameters)
            elif genome_dimension == 2: genome = G2DBinaryString(nparameters, nparameters2)
            else: raise ValueError("Dimensions > 2 are not supported")

        # Invalid option
        else: raise ValueError("Genome type must be 'list' or 'binary_string")

    # Set genome properties
    if parameter_range is not None: genome.setParams(rangemin=parameter_range.min, rangemax=parameter_range.max)
    if best_raw_score is not None: genome.setParams(bestrawscore=best_raw_score)
    if round_decimal is not None: genome.setParams(rounddecimal=round_decimal)

    # Set initializator
    if initializator is not None: genome.initializator.set(initializator)
    else:
        if isinstance(parameter_range, IntegerRange): genome.initializator.set(G1DListInitializatorInteger)
        elif isinstance(parameter_range, RealRange): genome.initializator.set(G1DListInitializatorReal)
        else: raise ValueError("Invalid parameter range")

    # Set mutation method
    if mutator is not None: genome.mutator.set(mutator)
    else:
        if isinstance(parameter_range, IntegerRange):
            if mutation_method == "range": genome.mutator.set(G1DListMutatorIntegerRange)
            elif mutation_method == "gaussian": genome.mutator.set(G1DListMutatorIntegerGaussian)
            elif mutation_method == "binary": genome.mutator.set(G1DListMutatorIntegerBinary)
            else: raise ValueError("Mutation method '" + mutation_method + "' not recognized")
        elif isinstance(parameter_range, RealRange):
            if mutation_method == "range": genome.mutator.set(G1DListMutatorRealRange)
            elif mutation_method == "gaussian": genome.mutator.set(G1DListMutatorRealGaussian)
            else: raise ValueError("Mutation method '" + mutation_method + "' not valid for genome of real values")
        else: raise ValueError("Invalid parameter range")

    # Set crossover
    if crossover is not None: genome.crossover.set(crossover)

    # The evaluator function (objective function)
    genome.evaluator.set(evaluator)

    # Genetic Algorithm Instance
    ga = GSimpleGA(genome)
    ga.terminationCriteria.set(RawScoreCriteria)
    ga.setMinimax(minimaxType[min_or_max])
    ga.setGenerations(ngenerations)
    if crossover_rate is not None: ga.setCrossoverRate(crossover_rate)
    ga.setPopulationSize(nindividuals)
    ga.setMutationRate(mutation_rate)

    # Set the callback function
    if callback is not None: ga.stepCallback.set(callback)

    # -----------------------------------------------------------------

    # Do the evolution
    with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR, filled_char=BAR_FILLED_CHAR, expected_size=ga.nGenerations, every=1, add_datetime=True) as bar:

        # Create partial function for the callback
        refresh_progress_bar_partial = partial(refresh_progress_bar, **{"bar": bar})

        # Set callback function to print generation
        ga.stepCallback.set(refresh_progress_bar_partial)

        # Evolve
        ga.evolve(freq_stats=stats_freq)

    # -----------------------------------------------------------------

    # Get the best individual
    best = ga.bestIndividual()

    # Return the best individual
    return best

# -----------------------------------------------------------------
