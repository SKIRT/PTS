#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ..genomes.list1d import G1DList
from ..genomes.list2d import G2DList
from ..genomes.binarystring1d import G1DBinaryString
from ..genomes.binarystring2d import G2DBinaryString
from ..core.crossovers import G1DListCrossoverSinglePoint, G1DListCrossoverTwoPoint, G1DListCrossoverUniform, G1DListCrossoverOX
from ..core.crossovers import G1DListCrossoverEdge, G1DListCrossoverCutCrossfill, G1DListCrossoverRealSBX, G1DListCrossoverMix
from ..core.crossovers import G2DListCrossoverUniform, G2DListCrossoverSingleVPoint, G2DListCrossoverSingleHPoint
from ..core.crossovers import G1DBinaryStringXSinglePoint, G1DBinaryStringXTwoPoint, G1DBinaryStringXUniform
from ..core.crossovers import G2DBinaryStringXSingleHPoint, G2DBinaryStringXSingleVPoint, G2DBinaryStringXUniform
from ..core.mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary
from ..core.mutators import G1DListMutatorRealGaussian, G1DListMutatorRealRange
from ..core.mutators import HeterogeneousListMutatorRealRange, HeterogeneousListMutatorRealGaussian
from ..core.mutators import HeterogeneousListMutatorIntegerRange, HeterogeneousListMutatorIntegerGaussian
from ..core.mutators import G1DBinaryStringMutatorSwap, G1DBinaryStringMutatorFlip
from ..core.mutators import G2DBinaryStringMutatorSwap, G2DBinaryStringMutatorFlip
from ..core.crossovers import G1DListCrossoverSinglePoint_origins, G1DListCrossoverTwoPoint_origins
from ..core.crossovers import G1DListCrossoverUniform_origins, G1DListCrossoverOX_origins
from ..core.crossovers import G1DListCrossoverEdge_origins, G1DListCrossoverCutCrossfill_origins
from ..core.crossovers import G1DListCrossoverRealSBX_origins, G1DListCrossoverMix_origins
from ..core.crossovers import G2DListCrossoverUniform_origins, G2DListCrossoverSingleHPoint_origins, G2DListCrossoverSingleVPoint_origins
from ..core.crossovers import G1DBinaryStringXSinglePoint_origins, G1DBinaryStringXTwoPoint_origins, G1DBinaryStringXUniform_origins
from ..core.crossovers import G2DBinaryStringXUniform_origins, G2DBinaryStringXSingleVPoint_origins, G2DBinaryStringXSingleHPoint_origins
from ..core.selectors import GRankSelector, GUniformSelector, GTournamentSelector, GTournamentSelectorAlternative, GRouletteWheel
from ..core.scaling import LinearScaling, SigmaTruncScaling, PowerLawScaling, BoltzmannScaling, ExponentialScaling, SaturatedScaling
from ..core.initializators import G1DListInitializatorReal, G1DListInitializatorInteger, HeterogeneousListInitializerReal
from ..core.initializators import G1DBinaryStringInitializator, HeterogeneousListInitializerInteger

# -----------------------------------------------------------------

# CROSSOVER

# -----------------------------------------------------------------

list_crossovers_1d = dict()
list_crossovers_1d["single_point"] = G1DListCrossoverSinglePoint
list_crossovers_1d["two_point"] =  G1DListCrossoverTwoPoint
list_crossovers_1d["uniform"] = G1DListCrossoverUniform
list_crossovers_1d["OX"] = G1DListCrossoverOX
list_crossovers_1d["edge"] = G1DListCrossoverEdge
list_crossovers_1d["cut_crossfill"] = G1DListCrossoverCutCrossfill
list_crossovers_1d["real_SBX"] = G1DListCrossoverRealSBX
list_crossovers_1d["mix"] = G1DListCrossoverMix

# -----------------------------------------------------------------

list_crossovers_2d = dict()
list_crossovers_2d["uniform"] = G2DListCrossoverUniform
list_crossovers_2d["single_vertical_point"] = G2DListCrossoverSingleVPoint
list_crossovers_2d["single_horizontal_point"] = G2DListCrossoverSingleHPoint

# -----------------------------------------------------------------

list_crossovers = dict()
list_crossovers[1] = list_crossovers_1d
list_crossovers[2] = list_crossovers_2d

# -----------------------------------------------------------------

binary_string_crossovers_1d = dict()
binary_string_crossovers_1d["single_point"] = G1DBinaryStringXSinglePoint
binary_string_crossovers_1d["two_point"] = G1DBinaryStringXTwoPoint
binary_string_crossovers_1d["uniform"] = G1DBinaryStringXUniform

# -----------------------------------------------------------------

binary_string_crossovers_2d = dict()
binary_string_crossovers_2d["uniform"] = G2DBinaryStringXUniform
binary_string_crossovers_2d["single_vertical_point"] = G2DBinaryStringXSingleVPoint
binary_string_crossovers_2d["single_horizontal_point"] = G2DBinaryStringXSingleHPoint

# -----------------------------------------------------------------

binary_string_crossovers = dict()
binary_string_crossovers[1] = binary_string_crossovers_1d
binary_string_crossovers[2] = binary_string_crossovers_2d

# -----------------------------------------------------------------

crossovers = dict()
crossovers["list"] = list_crossovers
crossovers["binary_string"] = binary_string_crossovers

# -----------------------------------------------------------------

def get_crossover(genome_type, dimension, method, genome_size=None):

    """
    This function ...
    :param genome_type:
    :param dimension:
    :param method:
    :param genome_size:
    :return:
    """

    # Check
    if genome_size is not None and genome_type == "list" and dimension == 1 and genome_size == 1:

        # Give a warning that we are not using the user specificied option
        log.warning("Mix crossover will be used because the genome size is only one")
        return G1DListCrossoverMix

    # Check whether the genome size is not zero
    if genome_size is not None and genome_size == 0: raise ValueError("The genome size cannot be zero")

    # Return the crossover
    return crossovers[genome_type][dimension][method]

# -----------------------------------------------------------------

def get_crossover_method(crossover_function):

    """
    This function ...
    :param crossover_function:
    :return:
    """

    for genome_type in crossovers:
        for dimension in crossovers[genome_type]:
            for method in crossovers[genome_type][dimension]:
                if get_crossover(genome_type, dimension, method) == crossover_function: return method
    #else: return None

    raise ValueError("Could not determine the crossover method for " + str(crossover_function.__name__))

# -----------------------------------------------------------------

# CROSSOVER ORIGINS

# -----------------------------------------------------------------

list_crossover_origins_1d = dict()
list_crossover_origins_1d["single_point"] = G1DListCrossoverSinglePoint_origins
list_crossover_origins_1d["two_point"] = G1DListCrossoverTwoPoint_origins
list_crossover_origins_1d["uniform"] = G1DListCrossoverUniform_origins
list_crossover_origins_1d["OX"] = G1DListCrossoverOX_origins
list_crossover_origins_1d["edge"] = G1DListCrossoverEdge_origins
list_crossover_origins_1d["cut_crossfill"] = G1DListCrossoverCutCrossfill_origins
list_crossover_origins_1d["real_SBX"] = G1DListCrossoverRealSBX_origins
list_crossover_origins_1d["mix"] = G1DListCrossoverMix_origins

# -----------------------------------------------------------------

list_crossover_origins_2d = dict()
list_crossover_origins_2d["uniform"] = G2DListCrossoverUniform_origins
list_crossover_origins_2d["single_vertical_point"] = G2DListCrossoverSingleVPoint_origins
list_crossover_origins_2d["single_horizontal_point"] = G2DListCrossoverSingleHPoint_origins

# -----------------------------------------------------------------

list_crossover_origins = dict()
list_crossover_origins[1] = list_crossover_origins_1d
list_crossover_origins[2] = list_crossover_origins_2d

# -----------------------------------------------------------------

binary_string_crossover_origins_1d = dict()
binary_string_crossover_origins_1d["single_point"] = G1DBinaryStringXSinglePoint_origins
binary_string_crossover_origins_1d["two_point"] = G1DBinaryStringXTwoPoint_origins
binary_string_crossover_origins_1d["uniform"] = G1DBinaryStringXUniform_origins

# -----------------------------------------------------------------

binary_string_crossover_origins_2d = dict()
binary_string_crossover_origins_2d["uniform"] = G2DBinaryStringXUniform_origins
binary_string_crossover_origins_2d["single_vertical_point"] = G2DBinaryStringXSingleVPoint_origins
binary_string_crossover_origins_2d["single_horizontal_point"] = G2DBinaryStringXSingleHPoint_origins

# -----------------------------------------------------------------

binary_string_crossover_origins = dict()
binary_string_crossover_origins[1] = binary_string_crossover_origins_1d
binary_string_crossover_origins[2] = binary_string_crossover_origins_2d

# -----------------------------------------------------------------

crossover_origins = dict()
crossover_origins["list"] = list_crossover_origins
crossover_origins["binary_string"] = binary_string_crossover_origins

# -----------------------------------------------------------------

def get_crossover_origins(genome_type, dimension, method):

    """
    This function ...
    :return:
    """

    return crossover_origins[genome_type][dimension][method]

# -----------------------------------------------------------------

# MUTATORS

# -----------------------------------------------------------------

list_mutators_1d_integer = dict()
list_mutators_1d_integer["range"] = G1DListMutatorIntegerRange
list_mutators_1d_integer["gaussian"] = G1DListMutatorIntegerGaussian
list_mutators_1d_integer["binary"] = G1DListMutatorIntegerBinary

# -----------------------------------------------------------------

list_mutators_1d_integer_hetero = dict()
list_mutators_1d_integer_hetero["range"] = HeterogeneousListMutatorIntegerRange
list_mutators_1d_integer_hetero["gaussian"] = HeterogeneousListMutatorIntegerGaussian

# -----------------------------------------------------------------

list_mutators_1d_real = dict()
list_mutators_1d_real["range"] = G1DListMutatorRealRange
list_mutators_1d_real["gaussian"] = G1DListMutatorRealGaussian

# -----------------------------------------------------------------

list_mutators_1d_real_hetero = dict()
list_mutators_1d_real_hetero["range"] = HeterogeneousListMutatorRealRange
list_mutators_1d_real_hetero["gaussian"] = HeterogeneousListMutatorRealGaussian

# -----------------------------------------------------------------

list_mutators_1d = dict()
list_mutators_1d["integer"] = list_mutators_1d_integer
list_mutators_1d["real"] = list_mutators_1d_real

# -----------------------------------------------------------------

list_mutators_1d_hetero = dict()
list_mutators_1d_hetero["integer"] = list_mutators_1d_integer_hetero
list_mutators_1d_hetero["real"] = list_mutators_1d_real_hetero

# -----------------------------------------------------------------

list_mutators_2d_integer = dict()

list_mutators_2d_integer_hetero = dict()

list_mutators_2d_real = dict()

list_mutators_2d_real_hetero = dict()

# -----------------------------------------------------------------

list_mutators_2d = dict()
list_mutators_2d["integer"] = list_mutators_2d_integer
list_mutators_2d["real"] = list_mutators_2d_real

list_mutators_2d_hetero = dict()
list_mutators_2d_hetero["integer"] = list_mutators_2d_integer
list_mutators_2d_hetero["real"] = list_mutators_2d_real

# -----------------------------------------------------------------

list_mutators = dict()
list_mutators[1] = dict()
list_mutators[1][False] = list_mutators_1d
list_mutators[1][True] = list_mutators_1d_hetero
list_mutators[2] = dict()
list_mutators[2][False] = list_mutators_2d
list_mutators[2][True] = list_mutators_2d_hetero

# -----------------------------------------------------------------

binary_string_mutators_1d = dict()
binary_string_mutators_1d["swap"] = G1DBinaryStringMutatorSwap
binary_string_mutators_1d["flip"] = G1DBinaryStringMutatorFlip

# -----------------------------------------------------------------

binary_string_mutators_2d = dict()
binary_string_mutators_2d["swap"] = G2DBinaryStringMutatorSwap
binary_string_mutators_2d["flip"] = G2DBinaryStringMutatorFlip

# -----------------------------------------------------------------

binary_string_mutators = dict()
binary_string_mutators[1] = binary_string_mutators_1d
binary_string_mutators[2] = binary_string_mutators_2d

# -----------------------------------------------------------------

mutators = dict()
mutators["list"] = list_mutators
mutators["binary_string"] = binary_string_mutators

# -----------------------------------------------------------------

def get_mutator(genome_type, dimension, method, parameter_type=None, hetero=False):

    """
    This function ...
    :param genome_type:
    :param dimension:
    :param method:
    :param parameter_type:
    :param hetero:
    :return:
    """

    # List genome (real or integer)
    if genome_type == "list":

        # Check whether parameter type is specified
        if parameter_type is None: raise ValueError("Parameter type has to be specified")
        return list_mutators[dimension][hetero][parameter_type][method]

    # Binary string genome
    elif genome_type == "binary_string": return binary_string_mutators[dimension][method]

    # Invalid
    else: raise ValueError("Invalid genome type: '" + genome_type + "'")

# -----------------------------------------------------------------

# GENOMES

# -----------------------------------------------------------------

genomes_1d = dict()
genomes_1d["list"] = G1DList
genomes_1d["binary_string"] = G1DBinaryString

# -----------------------------------------------------------------

genomes_2d = dict()
genomes_2d["list"] = G2DList
genomes_2d["binary_string"] = G2DBinaryString

# -----------------------------------------------------------------

genomes = dict()
genomes[1] = genomes_1d
genomes[2] = genomes_2d

# -----------------------------------------------------------------

def get_genome_class(dimension, genome_type):

    """
    This function ...
    :return:
    """

    return genomes[dimension][genome_type]

# -----------------------------------------------------------------

def create_genome(dimension, genome_type, nparameters, nparameters2=None, nbits=None):

    """
    This function ...
    :param dimension:
    :param genome_type:
    :param nparameters:
    :param nparameters2:
    :param nbits:
    :return:
    """

    # List genome
    if genome_type == "list":

        # Create 1D genome
        if dimension == 1: genome = G1DList(nparameters)
        elif dimension == 2:
            if nparameters2 is None: raise ValueError("Nparameters2 should be specified")
            genome = G2DList(nparameters, nparameters2)
        else: raise ValueError("Dimensions > 2 are not supported")

    # Binary genome
    elif genome_type == "binary_string":

        # Number of bits has to be specified
        if nbits is None: raise ValueError("Total number of bits has to be specified")

        # 1D
        if dimension == 1: genome = G1DBinaryString(nbits)

        # 2D
        elif dimension == 2: raise NotImplementedError("Not implemented")

        # Invalid
        else: raise ValueError("Dimensions > 2 are not supported")

    # Invalid
    else: raise ValueError("Genome type must be 'list' or 'binary_string")

    # Return the genome
    return genome

# -----------------------------------------------------------------

def get_genome_type(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    if isinstance(genome, G1DList) or isinstance(genome, G2DList): return "list"
    elif isinstance(genome, G1DBinaryString) or isinstance(genome, G2DBinaryString): return "binary_string"
    else: raise ValueError("Genome type not recognized: " + str(type(genome)))

# -----------------------------------------------------------------

def is_list_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return get_genome_type(genome) == "list"

# -----------------------------------------------------------------

def is_binary_string_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return get_genome_type(genome) == "binary_string"

# -----------------------------------------------------------------

def get_genome_dimension(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return genome.dimension

# -----------------------------------------------------------------

def is_1d_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return get_genome_dimension(genome) == 1

# -----------------------------------------------------------------

def is_2d_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return get_genome_dimension(genome) == 2

# -----------------------------------------------------------------

# SELECTION

# -----------------------------------------------------------------

def get_selector(method):

    """
    This function ...
    :param method:
    :return:
    """

    # Rank selector
    if method == "rank": return GRankSelector

    # Uniform selector
    elif method == "uniform": return GUniformSelector

    # Tournament selector
    elif method == "tournament": return GTournamentSelector

    # Alternative tournament selector
    elif method == "tournament_alternative": return GTournamentSelectorAlternative

    # Roulette wheel selector
    elif method == "roulette_wheel": return GRouletteWheel

    # Invalid
    else: raise ValueError("Invalid selector method")

# -----------------------------------------------------------------

def get_selection_method(selector):

    """"
    This function ...
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

# SCALING

# -----------------------------------------------------------------

def get_scaling(method):

    """
    This function ...
    :param method:
    :return:
    """

    # Linear scaling
    if method == "linear": return LinearScaling

    # Sigma truncation scaling
    elif method == "sigma_truncation": return SigmaTruncScaling

    # Power law scaling
    elif method == "power_law": return PowerLawScaling

    # Boltzmann scaling
    elif method == "boltzmann": return BoltzmannScaling

    # Exponential scaling
    elif method == "exponential": return ExponentialScaling

    # Saturated scaling
    elif method == "saturated": return SaturatedScaling

    # Invalid
    else: raise ValueError("Invalid scaling method '" + method + "'")

# -----------------------------------------------------------------

# INITIALIZATION

# -----------------------------------------------------------------

def get_initializator(genome_type, parameter_type=None, hetero=False):

    """
    This function ...
    :param genome_type:
    :param parameter_type:
    :param hetero:
    :return:
    """

    # Inform the user
    log.info("Getting the initializator type ...")

    if genome_type == "list": return get_list_initializator(parameter_type, hetero=hetero)
    elif genome_type == "binary_string": return get_binary_string_initializator()
    else: raise ValueError("Genome type not recognized")

# -----------------------------------------------------------------

def get_list_initializator(parameter_type, hetero=False):

    """
    This function ...
    :param parameter_type:
    :param hetero:
    :return:
    """

    # Integer type
    if parameter_type == "integer":

        if hetero: return HeterogeneousListInitializerInteger
        else: return G1DListInitializatorInteger

    # Real type
    elif parameter_type == "real":

        if hetero: return HeterogeneousListInitializerReal
        else: return G1DListInitializatorReal

    # Invalid
    else: raise ValueError("Invalid parameter type")

# -----------------------------------------------------------------

def get_binary_string_initializator():

    """
    This function ...
    :return:
    """

    return G1DBinaryStringInitializator

# -----------------------------------------------------------------
