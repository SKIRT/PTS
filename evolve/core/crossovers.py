#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.crossovers In this module we have the genetic operators of crossover (or recombination)
#  for each chromosome representation.

# -----------------------------------------------------------------

# Import standard modules
import math
import numpy as np

# Import other evolve modules
from . import utils

# Import the relevant PTS classes and modules
from ...core.tools.random import prng

# -----------------------------------------------------------------

def G1DBinaryStringXSinglePoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    cut = details

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size) # true means as father

    origins_sister[cut:] = False
    origins_brother[cut:] = False

    # Return
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DBinaryStringXSinglePoint(genome, **args):

    """
    The crossover of 1D Binary String, Single Point
    .. warning:: You can't use this crossover method for binary strings with length of 1.
    """

    return_details = args.pop("return_details", False)
    details = args.pop("details", None)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    if len(gMom) == 1: utils.raiseException("The Binary String have one element, can't use the Single Point Crossover method !", TypeError)

    if details is not None: cut = details
    else: cut = prng.randint(1, len(gMom))

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      sister[cut:] = gDad[cut:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      brother[cut:] = gMom[cut:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G1DBinaryStringXTwoPoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    cuts = details

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size)  # true means as father

    origins_sister[cuts[0]:cuts[1]] = False
    origins_brother[cuts[0]:cuts[1]] = False

    # Return the origins
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DBinaryStringXTwoPoint(genome, **args):

    """
    The 1D Binary String crossover, Two Point
    .. warning:: You can't use this crossover method for binary strings with length of 1.
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    if len(gMom) == 1:
      utils.raiseException("The Binary String have one element, can't use the Two Point Crossover method !", TypeError)

    cuts = [prng.randint(1, len(gMom)), prng.randint(1, len(gMom))]

    if cuts[0] > cuts[1]:
      utils.listSwapElement(cuts, 0, 1)

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      sister[cuts[0]:cuts[1]] = gDad[cuts[0]:cuts[1]]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      brother[cuts[0]:cuts[1]] = gMom[cuts[0]:cuts[1]]

    # Return the sister and brother
    if return_details: return sister, brother, cuts
    else: return sister, brother

# -----------------------------------------------------------------

def G1DBinaryStringXUniform_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    positions = details

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size)  # true means as father

    for position in positions:
        origins_sister[position] = False
        origins_brother[position] = False

    # Return the origins
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DBinaryStringXUniform(genome, **args):

    """
    The G1DList Uniform Crossover
    """

    return_details = args.pop("return_details", False)

    from . import constants

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    # Positions which were swapped between sister and brother
    positions = []

    for i in xrange(len(gMom)):

      if utils.randomFlipCoin(constants.CDefG1DBinaryStringUniformProb):
         positions.append(i)
         temp = sister[i]
         sister[i] = brother[i]
         brother[i] = temp

    # Return the sister and brother
    if return_details: return sister, brother, positions
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverSinglePoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    cut = details

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size)  # true means as father

    origins_sister[cut:] = False
    origins_brother[cut:] = False

    # Return
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DListCrossoverSinglePoint(genome, **args):

    """
    The crossover of G1DList, Single Point
    .. warning:: You can't use this crossover method for lists with just one element.
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    if len(gMom) == 1: raise RuntimeError("The 1D list has only one element, can't use the Single Point Crossover method")

    cut = prng.randint(1, len(gMom))

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      sister[cut:] = gDad[cut:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      brother[cut:] = gMom[cut:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverTwoPoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size)  # true means as father

    cuts = details

    origins_sister[cuts[0]:cuts[1]] = False
    origins_brother[cuts[0]:cuts[1]] = False

    # Return
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DListCrossoverTwoPoint(genome, **args):

    """
    The G1DList crossover, Two Point
    .. warning:: You can't use this crossover method for lists with just one element.
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    if len(gMom) == 1: raise RuntimeError("The 1D list has only one element, can't use the Two Point Crossover method")

    cuts = [prng.randint(1, len(gMom)), prng.randint(1, len(gMom))]

    if cuts[0] > cuts[1]:
      utils.listSwapElement(cuts, 0, 1)

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      sister[cuts[0]:cuts[1]] = gDad[cuts[0]:cuts[1]]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      brother[cuts[0]:cuts[1]] = gMom[cuts[0]:cuts[1]]

    # Return the sister and brother
    if return_details: return sister, brother, cuts
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverUniform_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    positions = details

    origins_sister = np.array([True] * size)  # true means as mother
    origins_brother = np.array([True] * size)  # true means as father

    for position in positions:
        origins_sister[position] = False
        origins_brother[position] = False

    # Return the origins
    return origins_sister, origins_brother

# -----------------------------------------------------------------

def G1DListCrossoverUniform(genome, **args):

    """
    The G1DList Uniform Crossover
    Each gene has a 50% chance of being swapped between mom and dad
    """

    return_details = args.pop("return_details", False)

    from . import constants

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    # Positions at which swapping occured
    positions = []

    for i in xrange(len(gMom)):
      if utils.randomFlipCoin(constants.CDefG1DListCrossUniformProb):
         positions.append(i)
         temp = sister[i]
         sister[i] = brother[i]
         brother[i] = temp

    # Return the sister and brother
    if return_details: return sister, brother, positions
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverMix_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    randoms = details

    #origins_sister = np.array([True] * size)  # true means as mother
    #origins_brother = np.array([True] * size)  # true means as father

    origins_sister = []
    origins_brother = []

    #for random in randoms:

    # Return the fractions of 'from mother' for sister and 'from father' for brother as arrays
    return np.array(randoms), np.array(randoms)

# -----------------------------------------------------------------

def G1DListCrossoverMix(genome, **args):

    """
    Each gene is a random linear combination between mom and dad
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    # The random numbers to make the mixes
    randoms = []

    for i in xrange(len(gMom)):

        # Generate random uniform number between 0 and 1
        ra = prng.uniform()

        # Add the random number
        randoms.append(ra)

        temp_s = sister[i]
        temp_b = brother[i]

        sister[i] = ra * temp_s + (1. - ra) * temp_b
        brother[i] = (1. - ra) * temp_s + ra * temp_b

    # Return sister and brother
    if return_details: return sister, brother, randoms
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverOX_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G1DListCrossoverOX(genome, **args):

    """
    The OX Crossover for G1DList  (order crossover)
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]
    listSize = len(gMom)

    c1, c2 = prng.randint(1, len(gMom)), prng.randint(1, len(gMom))

    while c1 == c2: c2 = prng.randint(1, len(gMom))

    if c1 > c2:
      h = c1
      c1 = c2
      c2 = h

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      P1 = [c for c in gMom[c2:] + gMom[:c2] if c not in gDad[c1:c2]]
      sister.genomeList = P1[listSize - c2:] + gDad[c1:c2] + P1[:listSize - c2]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      P2 = [c for c in gDad[c2:] + gDad[:c2] if c not in gMom[c1:c2]]
      brother.genomeList = P2[listSize - c2:] + gMom[c1:c2] + P2[:listSize - c2]

    assert listSize == len(sister)
    assert listSize == len(brother)

    # Return
    if return_details: return sister, brother, (c1,c2)
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverEdge_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G1DListCrossoverEdge(genome, **args):

    """ THe Edge Recombination crossover for G1DList (widely used for TSP problem)
    See more information in the `Edge Recombination Operator <http://en.wikipedia.org/wiki/Edge_recombination_operator>`_
    Wikipedia entry.
    """

    return_details = args.pop("return_details", False)

    gMom, sisterl = args["mom"], []
    gDad, brotherl = args["dad"], []

    mom_edges, dad_edges, merge_edges = utils.G1DListGetEdgesComposite(gMom, gDad)

    for c, u in (sisterl, set(gMom)), (brotherl, set(gDad)):

      curr = None

      for i in xrange(len(gMom)):
         curr = prng.choice(tuple(u)) if not curr else curr
         c.append(curr)
         u.remove(curr)
         d = [v for v in merge_edges.get(curr, []) if v in u]
         if d:
            curr = prng.choice(d)
         else:
            s = [v for v in mom_edges.get(curr, []) if v in u]
            s += [v for v in dad_edges.get(curr, []) if v in u]
            curr = prng.choice(s) if s else None

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    sister.genomeList = sisterl
    brother.genomeList = brotherl

    # Return the sister and brother
    if return_details: return sister, brother, (mom_edges, dad_edges, merge_edges)
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverCutCrossfill_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G1DListCrossoverCutCrossfill(genome, **args):

    """
    The crossover of G1DList, Cut and crossfill, for permutations
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    if len(gMom) == 1: utils.raiseException("The 1D List have one element, can't use the Single Point Crossover method !", TypeError)

    cut = prng.randint(1, len(gMom))

    if args["count"] >= 1:
      sister = gMom.clone()
      mother_part = gMom[0:cut]
      sister.resetStats()
      i = (len(sister) - cut)
      x = 0
      for v in gDad:
         if v in mother_part:
            continue
         if x >= i:
            break
         sister[cut + x] = v
         x += 1

    if args["count"] == 2:
      brother = gDad.clone()
      father_part = gDad[0:cut]
      brother.resetStats()
      i = (len(brother) - cut)
      x = 0
      for v in gMom:
         if v in father_part:
            continue
         if x >= i:
            break
         brother[cut + x] = v
         x += 1

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G1DListCrossoverRealSBX_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G1DListCrossoverRealSBX(genome, **args):

    """
    Experimental SBX Implementation - Follows the implementation in NSGA-II (Deb, et.al)
    Some implementation `reference <http://vision.ucsd.edu/~sagarwal/icannga.pdf>`_.
    And another reference to the `Simulated Binary Crossover <http://www.mitpressjournals.org/doi/abs/10.1162/106365601750190406>`_.
    .. warning:: This crossover method is Data Type Dependent, which means that
                must be used for 1D genome of real values.
    """

    return_details = args.pop("return_details", False)

    from . import constants

    EPS = constants.CDefG1DListSBXEPS
    # Crossover distribution index
    eta_c = constants.CDefG1DListSBXEtac

    gMom = args["mom"]
    gDad = args["dad"]

    # Get the variable bounds ('gDad' could have been used; but I love Mom:-))
    lb = gMom.getParam("rangemin", constants.CDefRangeMin)
    ub = gMom.getParam("rangemax", constants.CDefRangeMax)

    sister = gMom.clone()
    brother = gDad.clone()

    sister.resetStats()
    brother.resetStats()

    randoms = []
    swapped = []

    for i in range(0, len(gMom)):

      if math.fabs(gMom[i] - gDad[i]) > EPS:

         if gMom[i] > gDad[i]:
            #swap
            temp = gMom[i]
            gMom[i] = gDad[i]
            gDad[i] = temp

         #random number betwn. 0 & 1
         u = prng.random_sample()

         randoms.append(u)

         beta = 1.0 + 2 * (gMom[i] - lb) / (1.0 * (gDad[i] - gMom[i]))
         alpha = 2.0 - beta ** (-(eta_c + 1.0))

         if u <= (1.0 / alpha):
            beta_q = (u * alpha) ** (1.0 / ((eta_c + 1.0) * 1.0))
         else:
            beta_q = (1.0 / (2.0 - u * alpha)) ** (1.0 / (1.0 * (eta_c + 1.0)))

         brother[i] = 0.5 * ((gMom[i] + gDad[i]) - beta_q * (gDad[i] - gMom[i]))

         beta = 1.0 + 2.0 * (ub - gDad[i]) / (1.0 * (gDad[i] - gMom[i]))
         alpha = 2.0 - beta ** (-(eta_c + 1.0))

         if u <= (1.0 / alpha):
            beta_q = (u * alpha) ** (1.0 / ((eta_c + 1) * 1.0))
         else:
            beta_q = (1.0 / (2.0 - u * alpha)) ** (1.0 / (1.0 * (eta_c + 1.0)))

         sister[i] = 0.5 * ((gMom[i] + gDad[i]) + beta_q * (gDad[i] - gMom[i]))

         if brother[i] > ub:
            brother[i] = ub
         if brother[i] < lb:
            brother[i] = lb

         if sister[i] > ub:
            sister[i] = ub
         if sister[i] < lb:
            sister[i] = lb

         if prng.random_sample() > 0.5:

            swapped.append(i)

            # Swap
            temp = sister[i]
            sister[i] = brother[i]
            brother[i] = temp

      else:

         randoms.append(None)

         sister[i] = gMom[i]
         brother[i] = gDad[i]

    # Return the sister and brother
    if return_details: return sister, brother, (randoms, swapped)
    else: return sister, brother

# -----------------------------------------------------------------

def G2DListCrossoverUniform_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DListCrossoverUniform(genome, **args):

    """
    The G2DList Uniform Crossover
    """

    return_details = args.pop("return_details", False)

    from . import constants

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    h, w = gMom.getSize()

    swapped = []

    for i in xrange(h):
      for j in xrange(w):

         if utils.randomFlipCoin(constants.CDefG2DListCrossUniformProb):

            swapped.append((i, j))

            temp = sister.getItem(i, j)
            sister.setItem(i, j, brother.getItem(i, j))
            brother.setItem(i, j, temp)

    # Return the sister and brother
    if return_details: return sister, brother, swapped
    else: return sister, brother

# -----------------------------------------------------------------

def G2DListCrossoverSingleVPoint_origins(size, details):

    """
    THis function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DListCrossoverSingleVPoint(genome, **args):

    """
    The crossover of G2DList, Single Vertical Point
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    cut = prng.randint(1, gMom.getWidth())

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      for i in xrange(sister.getHeight()):
         sister[i][cut:] = gDad[i][cut:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      for i in xrange(brother.getHeight()):
         brother[i][cut:] = gMom[i][cut:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G2DListCrossoverSingleHPoint_origins(size, details):

    """
    THis function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DListCrossoverSingleHPoint(genome, **args):

    """
    The crossover of G2DList, Single Horizontal Point
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    cut = prng.randint(1, gMom.getHeight())

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      for i in xrange(cut, sister.getHeight()):
         sister[i][:] = gDad[i][:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      for i in xrange(brother.getHeight()):
         brother[i][:] = gMom[i][:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G2DBinaryStringXUniform_origins(size, details):

    """
    THis function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DBinaryStringXUniform(genome, **args):

    """
    The G2DBinaryString Uniform Crossover
    .. versionadded:: 0.6
      The *G2DBinaryStringXUniform* function
    """

    return_details = args.pop("return_details", False)

    from . import constants

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    sister = gMom.clone()
    brother = gDad.clone()
    sister.resetStats()
    brother.resetStats()

    h, w = gMom.getSize()

    swapped = []

    for i in xrange(h):
      for j in xrange(w):
         if utils.randomFlipCoin(constants.CDefG2DBinaryStringUniformProb):
            swapped.append((i,j))
            temp = sister.getItem(i, j)
            sister.setItem(i, j, brother.getItem(i, j))
            brother.setItem(i, j, temp)

    # Return the sister and brother
    if return_details: return sister, brother, swapped
    else: return sister, brother

# -----------------------------------------------------------------

def G2DBinaryStringXSingleVPoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DBinaryStringXSingleVPoint(genome, **args):

    """
    The crossover of G2DBinaryString, Single Vertical Point
    .. versionadded:: 0.6
      The *G2DBinaryStringXSingleVPoint* function
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    cut = prng.randint(1, gMom.getWidth())

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      for i in xrange(sister.getHeight()):
         sister[i][cut:] = gDad[i][cut:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      for i in xrange(brother.getHeight()):
         brother[i][cut:] = gMom[i][cut:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def G2DBinaryStringXSingleHPoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def G2DBinaryStringXSingleHPoint(genome, **args):

    """
    The crossover of G2DBinaryString, Single Horizontal Point
    .. versionadded:: 0.6
      The *G2DBinaryStringXSingleHPoint* function
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"]
    gDad = args["dad"]

    cut = prng.randint(1, gMom.getHeight())

    if args["count"] >= 1:
      sister = gMom.clone()
      sister.resetStats()
      for i in xrange(cut, sister.getHeight()):
         sister[i][:] = gDad[i][:]

    if args["count"] == 2:
      brother = gDad.clone()
      brother.resetStats()
      for i in xrange(brother.getHeight()):
         brother[i][:] = gMom[i][:]

    # Return the sister and brother
    if return_details: return sister, brother, cut
    else: return sister, brother

# -----------------------------------------------------------------

def GTreeCrossoverSinglePoint(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    return NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def GTreeCrossoverSinglePoint(genome, **args):

    """
    The crossover for GTree, Single Point
    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None
    gMom = args["mom"].clone()
    gDad = args["dad"].clone()

    gMom.resetStats()
    gDad.resetStats()

    node_mom_stack = []
    all_mom_nodes = []
    node_mom_tmp = None

    node_dad_stack = []
    all_dad_nodes = []
    node_dad_tmp = None

    node_mom_stack.append(gMom.getRoot())
    node_dad_stack.append(gDad.getRoot())

    while (len(node_mom_stack) > 0) and (len(node_dad_stack) > 0):
      node_mom_tmp = node_mom_stack.pop()
      node_dad_tmp = node_dad_stack.pop()

      if node_mom_tmp != gMom.getRoot():
         all_mom_nodes.append(node_mom_tmp)
         all_dad_nodes.append(node_dad_tmp)

      node_mom_stack.extend(node_mom_tmp.getChilds())
      node_dad_stack.extend(node_dad_tmp.getChilds())

    if len(all_mom_nodes) == 0 or len(all_dad_nodes) == 0:
      return (gMom, gDad)

    if len(all_dad_nodes) == 1:
      nodeDad = all_dad_nodes[0]
    else:
      nodeDad = prng.choice(all_dad_nodes)

    if len(all_mom_nodes) == 1:
      nodeMom = all_mom_nodes[0]
    else:
      nodeMom = prng.choice(all_mom_nodes)

    nodeMom_parent = nodeMom.getParent()
    nodeDad_parent = nodeDad.getParent()

    # Sister
    if args["count"] >= 1:
      sister = gMom
      nodeDad.setParent(nodeMom_parent)
      nodeMom_parent.replaceChild(nodeMom, nodeDad)
      sister.processNodes()

    # Brother
    if args["count"] == 2:
      brother = gDad
      nodeMom.setParent(nodeDad_parent)
      nodeDad_parent.replaceChild(nodeDad, nodeMom)
      brother.processNodes()

    # Return the sister and brother
    if return_details: return sister, brother, None # too complicated?
    else: return sister, brother

# -----------------------------------------------------------------

def GTreeCrossoverSinglePointStrict_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def GTreeCrossoverSinglePointStrict(genome, **args):

    """
    The crossover of Tree, Strict Single Point
    ..note:: This crossover method creates offspring with restriction of the
            *max_depth* parameter.
    Accepts the *max_attempt* parameter, *max_depth* (required), and
    the distr_leaft (>= 0.0 and <= 1.0), which represents the probability
    of leaf selection when findin random nodes for crossover.

    """

    return_details = args.pop("return_details", False)

    sister = None
    brother = None

    gMom = args["mom"].clone()
    gDad = args["dad"].clone()

    gMom.resetStats()
    gDad.resetStats()

    max_depth = gMom.getParam("max_depth", None)
    max_attempt = gMom.getParam("max_attempt", 10)
    distr_leaf = gMom.getParam("distr_leaf", None)

    if max_depth is None:
      utils.raiseException("You must specify the max_depth genome parameter !", ValueError)

    if max_depth < 0:
      utils.raiseException("The max_depth must be >= 1, if you want to use GTreeCrossoverSinglePointStrict crossover !", ValueError)

    momRandom = None
    dadRandom = None

    for i in xrange(max_attempt):

      if distr_leaf is None:
         dadRandom = gDad.getRandomNode()
         momRandom = gMom.getRandomNode()
      else:
         if utils.randomFlipCoin(distr_leaf):
            momRandom = gMom.getRandomNode(1)
         else:
            momRandom = gMom.getRandomNode(2)

         if utils.randomFlipCoin(distr_leaf):
            dadRandom = gDad.getRandomNode(1)
         else:
            dadRandom = gDad.getRandomNode(2)

      assert momRandom is not None
      assert dadRandom is not None

      # Optimize here
      mH = gMom.getNodeHeight(momRandom)
      dH = gDad.getNodeHeight(dadRandom)

      mD = gMom.getNodeDepth(momRandom)
      dD = gDad.getNodeDepth(dadRandom)

      # The depth of the crossover is greater than the max_depth
      if (dD + mH <= max_depth) and (mD + dH <= max_depth):
         break

    if i == (max_attempt - 1):
      assert gMom.getHeight() <= max_depth
      return (gMom, gDad)
    else:
      nodeMom, nodeDad = momRandom, dadRandom

    nodeMom_parent = nodeMom.getParent()
    nodeDad_parent = nodeDad.getParent()

    # Sister
    if args["count"] >= 1:
      sister = gMom
      nodeDad.setParent(nodeMom_parent)

      if nodeMom_parent is None:
         sister.setRoot(nodeDad)
      else:
         nodeMom_parent.replaceChild(nodeMom, nodeDad)
      sister.processNodes()
      assert sister.getHeight() <= max_depth

    # Brother
    if args["count"] == 2:
      brother = gDad
      nodeMom.setParent(nodeDad_parent)

      if nodeDad_parent is None:
         brother.setRoot(nodeMom)
      else:
         nodeDad_parent.replaceChild(nodeDad, nodeMom)
      brother.processNodes()
      assert brother.getHeight() <= max_depth

    # Return the sister and brother
    if return_details: return sister, brother, None # too complicated?
    else: return sister, brother

# -----------------------------------------------------------------

def GTreeGPCrossoverSinglePoint_origins(size, details):

    """
    This function ...
    :param size:
    :param details:
    :return:
    """

    raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------

def GTreeGPCrossoverSinglePoint(genome, **args):

    """
    The crossover of the GTreeGP, Single Point for Genetic Programming
    ..note:: This crossover method creates offspring with restriction of the
            *max_depth* parameter.
    Accepts the *max_attempt* parameter, *max_depth* (required).
    """

    return_details = args.pop("return_details", False)

    from . import constants

    sister = None
    brother = None

    gMom = args["mom"].clone()
    gDad = args["dad"].clone()

    gMom.resetStats()
    gDad.resetStats()

    max_depth = gMom.getParam("max_depth", None)
    max_attempt = gMom.getParam("max_attempt", 15)

    if max_depth is None:
      utils.raiseException("You must specify the max_depth genome parameter !", ValueError)

    if max_depth < 0:
      utils.raiseException("The max_depth must be >= 1, if you want to use GTreeCrossoverSinglePointStrict crossover !", ValueError)

    momRandom = None
    dadRandom = None

    for i in xrange(max_attempt):

      dadRandom = gDad.getRandomNode()

      if dadRandom.getType() == constants.nodeType["TERMINAL"]:
         momRandom = gMom.getRandomNode(1)
      elif dadRandom.getType() == constants.nodeType["NONTERMINAL"]:
         momRandom = gMom.getRandomNode(2)

      mD = gMom.getNodeDepth(momRandom)
      dD = gDad.getNodeDepth(dadRandom)

      # Two nodes are root
      if mD == 0 and dD == 0:
         continue

      mH = gMom.getNodeHeight(momRandom)
      if dD + mH > max_depth:
         continue

      dH = gDad.getNodeHeight(dadRandom)
      if mD + dH > max_depth:
         continue

      break

    if i == (max_attempt - 1):
      assert gMom.getHeight() <= max_depth
      return (gMom, gDad)
    else:
      nodeMom, nodeDad = momRandom, dadRandom

    nodeMom_parent = nodeMom.getParent()
    nodeDad_parent = nodeDad.getParent()

    # Sister
    if args["count"] >= 1:
      sister = gMom
      nodeDad.setParent(nodeMom_parent)

      if nodeMom_parent is None:
         sister.setRoot(nodeDad)
      else:
         nodeMom_parent.replaceChild(nodeMom, nodeDad)
      sister.processNodes()
      assert sister.getHeight() <= max_depth

    # Brother
    if args["count"] == 2:
      brother = gDad
      nodeMom.setParent(nodeDad_parent)

      if nodeDad_parent is None:
         brother.setRoot(nodeMom)
      else:
         nodeDad_parent.replaceChild(nodeDad, nodeMom)
      brother.processNodes()
      assert brother.getHeight() <= max_depth

    # Return the sister and brother
    if return_details: return sister, brother, None # too complicated?
    else: return sister, brother

# -----------------------------------------------------------------
