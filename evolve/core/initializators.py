#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.Initializers In this module we have the genetic operators of initialization for each
#  chromosome representation, the most part of initialization is done by
#  choosing random data.
#
#.. note:: In Pyevolve, the Initializator defines the data type that will
#          be used on the chromosome, for example, the :func:`G1DListInitializatorInteger`
#          will initialize the G1DList with Integers.
#

# -----------------------------------------------------------------

# Import standard modules
import math

# Import other evolve modules
from . import tree
from . import utils

# Import the relevant PTS classes and modules
from ...core.tools.random import prng

# -----------------------------------------------------------------

def G1DBinaryStringInitializator(genome, **args):

    """
    1D Binary String initializator
    """

    genome.genomeList = [prng.choice((0, 1)) for _ in xrange(genome.getListSize())]

# -----------------------------------------------------------------

def G2DBinaryStringInitializator(genome, **args):

    """ Integer initialization function of 2D Binary String
    .. versionadded:: 0.6
       The *G2DBinaryStringInitializator* function
    """

    genome.clearString()

    for i in xrange(genome.getHeight()):
        for j in xrange(genome.getWidth()):
            random_gene = prng.choice((0, 1))
            genome.setItem(i, j, random_gene)

# -----------------------------------------------------------------

def G1DListInitializatorAllele(genome, **args):

    """ Allele initialization function of G1DList
    To use this initializator, you must specify the *allele* genome parameter with the
    :class:`GAllele.GAlleles` instance.
    """

    allele = genome.getParam("allele", None)
    if allele is None:
        utils.raiseException("to use the G1DListInitializatorAllele, you must specify the 'allele' parameter")

    genome.genomeList = [allele[i].getRandomAllele() for i in xrange(genome.getListSize())]

# -----------------------------------------------------------------

def G1DListInitializatorInteger(genome, **args):

    """
    Integer initialization function of G1DList
    This initializator accepts the *rangemin* and *rangemax* genome parameters.
    """

    range_min = genome.getParam("rangemin", 0)
    range_max = genome.getParam("rangemax", 100)

    genome.genomeList = [prng.randint(range_min, range_max + 1) for i in xrange(genome.getListSize())] # HERE IT SHOULD BE INCLUSIVE

# -----------------------------------------------------------------

def G1DListInitializatorReal(genome, **args):

    """
    Real initialization function of G1DList
    This initializator accepts the *rangemin* and *rangemax* genome parameters.
    """

    range_min = genome.getParam("rangemin", 0)
    range_max = genome.getParam("rangemax", 100)

    genome.genomeList = [prng.uniform(range_min, range_max) for i in xrange(genome.getListSize())]

# -----------------------------------------------------------------

def HeterogeneousListInitializerReal(genome, **args):

    """
    :param genome:
    :param args:
    :return:
    """

    genome.genomeList = []

    for index in xrange(genome.getListSize()):

        range_min = genome.getParam("minima")[index]
        range_max = genome.getParam("maxima")[index]

        genome.genomeList.append(prng.uniform(range_min, range_max))

# -----------------------------------------------------------------

def HeterogeneousListInitializerInteger(genome, **args):

    """
    This function ...
    :param genome:
    :param args:
    :return:
    """

    genome.genomeList = []

    for index in xrange(genome.getListSize()):

        range_min = int(math.ceil(genome.getParam("minima")[index]))
        range_max = int(math.floor(genome.getParam("maxima")[index]))

        genome.genomeList.append(prng.randint(range_min, range_max + 1)) # HERE IT SHOULD BE INCLUSIVE

# -----------------------------------------------------------------

def G2DListInitializatorInteger(genome, **args):

    """ Integer initialization function of G2DList
    This initializator accepts the *rangemin* and *rangemax* genome parameters.
    """

    genome.clearList()

    for i in xrange(genome.getHeight()):
        for j in xrange(genome.getWidth()):
            randomInteger = prng.randint(genome.getParam("rangemin", 0), # HERE IT SHOULD BE INCLUSIVE
                                         genome.getParam("rangemax", 100) + 1) # HERE IT SHOULD BE INCLUSIVE
            genome.setItem(i, j, randomInteger)

# -----------------------------------------------------------------

def G2DListInitializatorReal(genome, **args):

    """ Integer initialization function of G2DList
    This initializator accepts the *rangemin* and *rangemax* genome parameters.
    """

    genome.clearList()

    for i in xrange(genome.getHeight()):
        for j in xrange(genome.getWidth()):
            randomReal = prng.uniform(genome.getParam("rangemin", 0), genome.getParam("rangemax", 100))
            genome.setItem(i, j, randomReal)

# -----------------------------------------------------------------

def G2DListInitializatorAllele(genome, **args):

    """ Allele initialization function of G2DList
    To use this initializator, you must specify the *allele* genome parameter with the
    :class:`GAllele.GAlleles` instance.
    .. warning:: the :class:`GAllele.GAlleles` instance must have the homogeneous flag enabled
    """

    allele = genome.getParam("allele", None)
    if allele is None:
        utils.raiseException("to use the G2DListInitializatorAllele, you must specify the 'allele' parameter")

    if not allele.homogeneous:
        utils.raiseException("to use the G2DListInitializatorAllele, the 'allele' must be homogeneous")

    genome.clearList()

    for i in xrange(genome.getHeight()):
        for j in xrange(genome.getWidth()):
            random_allele = allele[0].getRandomAllele()
            genome.setItem(i, j, random_allele)

# -----------------------------------------------------------------

def GTreeInitializatorInteger(genome, **args):

    """ Integer initialization function of GTree
    This initializator accepts the *rangemin* and *rangemax* genome parameters.
    It accepts the following parameters too:

    *max_depth*
       The max depth of the tree

    *max_siblings*
       The number of maximum siblings of an node

    *method*
       The method, accepts "grow", "full" or "ramped".

    .. versionadded:: 0.6
       The *GTreeInitializatorInteger* function.
    """

    max_depth = genome.getParam("max_depth", 5)
    max_siblings = genome.getParam("max_siblings", 2)

    range_min = genome.getParam("rangemin", 0)
    range_max = genome.getParam("rangemax", 100)

    lambda_generator = lambda: prng.randint(range_min, range_max + 1) # HERE IT SHOULD BE INCLUSIVE

    method = genome.getParam("method", "grow")

    if method == "grow":
        root = tree.buildGTreeGrow(0, lambda_generator, max_siblings, max_depth)
    elif method == "full":
        root = tree.buildGTreeFull(0, lambda_generator, max_siblings, max_depth)
    elif method == "ramped":
        if utils.randomFlipCoin(0.5):
            root = tree.buildGTreeGrow(0, lambda_generator, max_siblings, max_depth)
        else:
            root = tree.buildGTreeFull(0, lambda_generator, max_siblings, max_depth)
    else:
        utils.raiseException("Unknown tree initialization method [%s] !" % method)

    genome.setRoot(root)
    genome.processNodes()
    assert genome.getHeight() <= max_depth

# -----------------------------------------------------------------

def GTreeInitializatorAllele(genome, **args):

    """ Allele initialization function of GTree
    To use this initializator, you must specify the *allele* genome parameter with the
    :class:`GAllele.GAlleles` instance.
    .. warning:: the :class:`GAllele.GAlleles` instance **must** have the homogeneous flag enabled
    .. versionadded:: 0.6
       The *GTreeInitializatorAllele* function.
    """

    max_depth = genome.getParam("max_depth", 5)
    max_siblings = genome.getParam("max_siblings", 2)
    method = genome.getParam("method", "grow")

    allele = genome.getParam("allele", None)
    if allele is None:
        utils.raiseException("to use the GTreeInitializatorAllele, you must specify the 'allele' parameter")

    if not allele.homogeneous:
        utils.raiseException("to use the GTreeInitializatorAllele, the 'allele' must be homogeneous")

    if method == "grow":
        root = tree.buildGTreeGrow(0, allele[0].getRandomAllele, max_siblings, max_depth)
    elif method == "full":
        root = tree.buildGTreeFull(0, allele[0].getRandomAllele, max_siblings, max_depth)
    elif method == "ramped":
        if utils.randomFlipCoin(0.5):
            root = tree.buildGTreeGrow(0, allele[0].getRandomAllele, max_siblings, max_depth)
        else:
            root = tree.buildGTreeFull(0, allele[0].getRandomAllele, max_siblings, max_depth)
    else:
        utils.raiseException("Unknown tree initialization method [%s] !" % method)

    genome.setRoot(root)
    genome.processNodes()
    assert genome.getHeight() <= max_depth

# -----------------------------------------------------------------

def GTreeGPInitializator(genome, **args):

    """This initializator accepts the follow parameters:
    *max_depth*
       The max depth of the tree
    *method*
       The method, accepts "grow", "full" or "ramped"
    .. versionadded:: 0.6
       The *GTreeGPInitializator* function.
    """

    max_depth = genome.getParam("max_depth", 5)
    method = genome.getParam("method", "grow")
    ga_engine = args["ga_engine"]

    if method == "grow":
        root = tree.buildGTreeGPGrow(ga_engine, 0, max_depth)
    elif method == "full":
        root = tree.buildGTreeGPFull(ga_engine, 0, max_depth)
    elif method == "ramped":
        if utils.randomFlipCoin(0.5):
            root = tree.buildGTreeGPFull(ga_engine, 0, max_depth)
        else:
            root = tree.buildGTreeGPGrow(ga_engine, 0, max_depth)
    else:
        utils.raiseException("Unknown tree initialization method [%s] !" % method)

    genome.setRoot(root)
    genome.processNodes()
    assert genome.getHeight() <= max_depth

# -----------------------------------------------------------------
