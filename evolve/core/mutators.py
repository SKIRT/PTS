#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.mutators In this module we have the genetic operators of mutation for each chromosome
#  representation.
#

# -----------------------------------------------------------------

# Import standard modules
import math

# Import other evolve modules
from . import utils
from . import tree

# Import the relevant PTS classes and modules
from ...core.tools.random import prng

# -----------------------------------------------------------------

def G1DBinaryStringMutatorSwap(genome, **args):

   """
   The 1D Binary String Swap Mutator
   """

   if args["pmut"] <= 0.0:
      return 0
   stringLength = len(genome)
   mutations = args["pmut"] * (stringLength)

   if mutations < 1.0:
      mutations = 0
      for it in xrange(stringLength):
         if utils.randomFlipCoin(args["pmut"]):
            utils.listSwapElement(genome, it, prng.randint(0, stringLength))
            mutations += 1

   else:
      for it in xrange(int(round(mutations))):
         utils.listSwapElement(genome, prng.randint(0, stringLength), prng.randint(0, stringLength))

   return int(mutations)

# -----------------------------------------------------------------

def G1DBinaryStringMutatorFlip(genome, **args):

   """
   The classical flip mutator for binary strings
   """

   if args["pmut"] <= 0.0:
      return 0
   stringLength = len(genome)
   mutations = args["pmut"] * (stringLength)

   if mutations < 1.0:
      mutations = 0
      for it in xrange(stringLength):
         if utils.randomFlipCoin(args["pmut"]):
            if genome[it] == 0:
               genome[it] = 1
            else:
               genome[it] = 0
            mutations += 1

   else:
      for it in xrange(int(round(mutations))):
         which = prng.randint(0, stringLength)
         if genome[which] == 0:
            genome[which] = 1
         else:
            genome[which] = 0

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorSwap(genome, **args):

   """ The mutator of G1DList, Swap Mutator
   .. note:: this mutator is :term:`Data Type Independent`
   """

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * listSize

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            utils.listSwapElement(genome, it, prng.randint(0, listSize))
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         utils.listSwapElement(genome, prng.randint(0, listSize), prng.randint(0, listSize))

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorSIM(genome, **args):

   """
   The mutator of G1DList, Simple Inversion Mutation
   .. note:: this mutator is :term:`Data Type Independent`
   """

   mutations = 0
   if args["pmut"] <= 0.0:
      return 0

   cuts = [prng.randint(0, len(genome) + 1), prng.randint(0, len(genome) + 1)] # HERE IT SHOULD BE INCLUSIVE

   if cuts[0] > cuts[1]:
      utils.listSwapElement(cuts, 0, 1)

   if (cuts[1] - cuts[0]) <= 0:
      cuts[1] = prng.randint(cuts[0], len(genome) + 1) # HERE IT SHOULD BE INCLUSIVE

   if utils.randomFlipCoin(args["pmut"]):
      part = genome[cuts[0]:cuts[1]]
      if len(part) == 0:
         return 0
      part.reverse()
      genome[cuts[0]:cuts[1]] = part
      mutations += 1

   return mutations

# -----------------------------------------------------------------

def G1DListMutatorIntegerRange(genome, **args):

   """
   Simple integer range mutator for G1DList
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * listSize

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            genome[it] = prng.randint(genome.getParam("rangemin", constants.CDefRangeMin), # HERE IT SHOULD BE INCLUSIVE
                                      genome.getParam("rangemax", constants.CDefRangeMax) + 1) # HERE IT SHOULD BE INCLUSIVE
            mutations += 1

   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         genome[which_gene] = prng.randint(genome.getParam("rangemin", constants.CDefRangeMin), # HERE IT SHOULD BE INCLUSIVE
                                           genome.getParam("rangemax", constants.CDefRangeMax) + 1) # HERE IT SHOULD BE INCLUSIVE

   return int(mutations)

# -----------------------------------------------------------------

def HeterogeneousListMutatorIntegerRange(genome, **args):

   """
   This function ...
   :param genome:
   :param args:
   :return:
   """

   if args["pmut"] <= 0.0:
       return 0

   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   if mutations < 1.0:

       mutations = 0
       for it in xrange(listSize):

           if utils.randomFlipCoin(args["pmut"]):

               genome[it] = prng.randint(int(math.ceil(genome.getParam("minima")[it])), int(math.floor(genome.getParam("maxima")[it])) + 1) # HERE IT SHOULD BE INCLUSIVE
               mutations += 1

   else:

       for it in xrange(int(round(mutations))):

            which_gene = prng.randint(0, listSize)
            genome[which_gene] = prng.randint(int(math.ceil(genome.getParam("minima")[which_gene])), int(math.floor(genome.getParam("maxima")[which_gene])) + 1) # HERE IT SHOULD BE INCLUSIVE

   # Return the number of mutations
   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorRealRange(genome, **args):

   """ Simple real range mutator for G1DList
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            genome[it] = prng.uniform(genome.getParam("rangemin", constants.CDefRangeMin),
                                      genome.getParam("rangemax", constants.CDefRangeMax))
            mutations += 1

   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         genome[which_gene] = prng.uniform(genome.getParam("rangemin", constants.CDefRangeMin),
                                           genome.getParam("rangemax", constants.CDefRangeMax))

   return int(mutations)

# -----------------------------------------------------------------

def HeterogeneousListMutatorRealRange(genome, **args):

    """
    Real range mutator for HeterogeneousList
    :param genome:
    :param args:
    :return:
    """

    if args["pmut"] <= 0.0:
        return 0

    listSize = len(genome)
    mutations = args["pmut"] * (listSize)

    if mutations < 1.0:

        mutations = 0
        for it in xrange(listSize):

            if utils.randomFlipCoin(args["pmut"]):

                genome[it] = prng.uniform(genome.getParam("minima")[it], genome.getParam("maxima")[it])
                mutations += 1

    else:

        for it in xrange(int(round(mutations))):

            which_gene = prng.randint(0, listSize)

            genome[which_gene] = prng.uniform(genome.getParam("minima")[which_gene], genome.getParam("maxima")[which_gene])

    return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorIntegerGaussianGradient(genome, **args):

   """ A gaussian mutator for G1DList of Integers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. The
   random distribution is set with mu=1.0 and std=0.0333

   Same as IntegerGaussian, except that this uses relative gradient rather than
   absolute gaussian. A value is randomly generated about gauss(mu=1, sigma=.0333)
   and multiplied by the gene to drift it up or down (depending on what side of
   1 the random value falls on) and cast to integer
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   mu = constants.CDefGaussianGradientMU
   sigma = constants.CDefGaussianGradientSIGMA

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            final_value = int(genome[it] * abs(prng.normal(mu, sigma)))

            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

            genome[it] = final_value
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         final_value = int(genome[which_gene] * abs(prng.normal(mu, sigma)))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome[which_gene] = final_value

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorIntegerGaussian(genome, **args):

   """ A gaussian mutator for G1DList of Integers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   mu = genome.getParam("gauss_mu")
   sigma = genome.getParam("gauss_sigma")

   if mu is None:
      mu = constants.CDefG1DListMutIntMU

   if sigma is None:
      sigma = constants.CDefG1DListMutIntSIGMA

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            final_value = genome[it] + int(prng.normal(mu, sigma))

            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

            genome[it] = final_value
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         final_value = genome[which_gene] + int(prng.normal(mu, sigma))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome[which_gene] = final_value

   return int(mutations)

# -----------------------------------------------------------------

def HeterogeneousListMutatorIntegerGaussian(genome, **args):

    """
    This function ...
    :param genome:
    :param args:
    :return:
    """

    if args["pmut"] <= 0.0:
        return 0
    listSize = len(genome)
    mutations = args["pmut"] * (listSize)

    if mutations < 1.0:
        mutations = 0

        # Let the fact whether we do a mutation for each genome depend on a random 'coin flip'
        for it in xrange(listSize):

            # Get the mu and the sigma
            mu = genome.getParam("centers")[it]
            sigma = genome.getParam("sigmas")[it]

            if utils.randomFlipCoin(args["pmut"]):

                final_value = genome[it] + int(prng.normal(mu, sigma))

                final_value = min(final_value, int(math.floor(genome.getParam("maxima")[it])))
                final_value = max(final_value, int(math.ceil(genome.getParam("minima")[it])))

                genome[it] = final_value
                mutations += 1
    else:

        # Do a specific number of mutations
        for it in xrange(int(round(mutations))):

            # Get the mu and the sigma
            mu = genome.getParam("centers")[it]
            sigma = genome.getParam("sigmas")[it]

            which_gene = prng.randint(0, listSize)
            final_value = genome[which_gene] + int(prng.normal(mu, sigma))

            final_value = min(final_value, int(math.floor(genome.getParam("maxima")[which_gene])))
            final_value = max(final_value, int(math.ceil(genome.getParam("minima")[which_gene])))

            genome[which_gene] = final_value

    # Return the number of mutations
    return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorRealGaussian(genome, **args):

   """
   The mutator of G1DList, Gaussian Mutator
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   mu = genome.getParam("gauss_mu")
   sigma = genome.getParam("gauss_sigma")

   if mu is None:
      mu = constants.CDefG1DListMutRealMU

   if sigma is None:
      sigma = constants.CDefG1DListMutRealSIGMA

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            final_value = genome[it] + prng.normal(mu, sigma)

            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

            genome[it] = final_value
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         final_value = genome[which_gene] + prng.normal(mu, sigma)

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome[which_gene] = final_value

   return int(mutations)

# -----------------------------------------------------------------

def HeterogeneousListMutatorRealGaussian(genome, **args):

   """
   Heteregeneous version of real gaussian list mutator
   """

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   #mu = genome.getParam("gauss_mu")
   #sigma = genome.getParam("gauss_sigma")

   #if mu is None:
   #   mu = constants.CDefG1DListMutRealMU

   #if sigma is None:
   #   sigma = constants.CDefG1DListMutRealSIGMA

   if mutations < 1.0:
      mutations = 0

      # Let the fact whether we do a mutation for each genome depend on a random 'coin flip'
      for it in xrange(listSize):

         # Get the mu and the sigma
         mu = genome.getParam("centers")[it]
         sigma = genome.getParam("sigmas")[it]

         if utils.randomFlipCoin(args["pmut"]):

            final_value = genome[it] + prng.normal(mu, sigma)

            final_value = min(final_value, genome.getParam("maxima")[it])
            final_value = max(final_value, genome.getParam("minima")[it])

            genome[it] = final_value
            mutations += 1
   else:

      # Do a specific number of mutations
      for it in xrange(int(round(mutations))):

         # Get the mu and the sigma
         mu = genome.getParam("centers")[it]
         sigma = genome.getParam("sigmas")[it]

         which_gene = prng.randint(0, listSize)
         final_value = genome[which_gene] + prng.normal(mu, sigma)

         final_value = min(final_value, genome.getParam("maxima")[which_gene])
         final_value = max(final_value, genome.getParam("minima")[which_gene])

         genome[which_gene] = final_value

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorRealGaussianGradient(genome, **args):

   """ The mutator of G1DList, Gaussian Gradient Mutator
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. The
   random distribution is set with mu=1.0 and std=0.0333

   The difference between this routine and the normal Gaussian Real is that the
   other function generates a gaussian value and adds it to the value. If the
   mu is 0, and the std is 1, a typical value could be 1.8 or -0.5. These small
   values are fine if your range is 0-10, but if your range is much larger, like
   0-100,000, a relative gradient makes sense.

   This routine generates a gaussian value with mu=1.0 and std=0.0333 and then
   the gene is multiplied by this value. This will cause the gene to drift
   no matter how large it is.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   mu = constants.CDefGaussianGradientMU
   sigma = constants.CDefGaussianGradientSIGMA

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            final_value = genome[it] * abs(prng.normal(mu, sigma))

            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

            genome[it] = final_value
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         final_value = genome[which_gene] * abs(prng.normal(mu, sigma))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome[which_gene] = final_value

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorIntegerBinary(genome, **args):

   """ The mutator of G1DList, the binary mutator
   This mutator will random change the 0 and 1 elements of the 1D List.

   """

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * (listSize)

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            if genome[it] == 0:
               genome[it] = 1
            elif genome[it] == 1:
               genome[it] = 0

            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         if genome[which_gene] == 0:
            genome[which_gene] = 1
         elif genome[which_gene] == 1:
            genome[which_gene] = 0

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorAllele(genome, **args):

   """ The mutator of G1DList, Allele Mutator
   To use this mutator, you must specify the *allele* genome parameter with the
   :class:`GAllele.GAlleles` instance.
   """

   if args["pmut"] <= 0.0:
      return 0
   listSize = len(genome)
   mutations = args["pmut"] * listSize

   allele = genome.getParam("allele", None)
   if allele is None:
      utils.raiseException("to use the G1DListMutatorAllele, you must specify the 'allele' parameter", TypeError)

   if mutations < 1.0:
      mutations = 0
      for it in xrange(listSize):
         if utils.randomFlipCoin(args["pmut"]):
            new_val = allele[it].getRandomAllele()
            genome[it] = new_val
            mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_gene = prng.randint(0, listSize)
         new_val = allele[which_gene].getRandomAllele()
         genome[which_gene] = new_val

   return int(mutations)

# -----------------------------------------------------------------

def G1DListMutatorAlleleGaussian(genome, **arguments):

    """An allele-based mutator based on G1DListMutatorRealGaussian.
    Accepts the parameter *gauss_mu* and the *gauss_sigma* which
    respectively represents the mean and the std. dev. of the random
    distribution.
    """

    from . import constants

    if arguments["pmut"] <= 0.0:
        return 0
    listSize = len(genome)
    mutations = arguments["pmut"] * listSize

    mu = genome.getParam("gauss_mu")
    sigma = genome.getParam("gauss_sigma")
    if mu is None:
        mu = constants.CDefG1DListMutRealMU
    if sigma is None:
        sigma = constants.CDefG1DListMutRealSIGMA

    allele = genome.getParam("allele", None)
    if allele is None:
        utils.raiseException("to use this mutator, you must specify the 'allele' parameter", TypeError)

    if mutations < 1.0:
        mutations = 0
        for it in xrange(listSize):
            if utils.randomFlipCoin(arguments["pmut"]):
                final_value = genome[it] + prng.normal(mu, sigma)
                assert len(allele[it].beginEnd) == 1, "only single ranges are supported"
                rangemin, rangemax = allele[it].beginEnd[0]
                final_value = min(final_value, rangemax)
                final_value = max(final_value, rangemin)
                genome[it] = final_value
                mutations += 1
    else:
        for it in xrange(int(round(mutations))):
            which_gene = prng.randint(0, listSize)
            final_value = genome[which_gene] + prng.normal(mu, sigma)
            assert len(allele[which_gene].beginEnd) == 1, "only single ranges are supported"
            rangemin, rangemax = allele[which_gene].beginEnd[0]
            final_value = min(final_value, rangemax)
            final_value = max(final_value, rangemin)
            genome[which_gene] = final_value
    return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorSwap(genome, **args):

   """ The mutator of G1DList, Swap Mutator
   .. note:: this mutator is :term:`Data Type Independent`
   """

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   if mutations < 1.0:
      mutations = 0
      for i in xrange(height):
         for j in xrange(width):
            if utils.randomFlipCoin(args["pmut"]):
               index_b = (prng.randint(0, height), prng.randint(0, width))
               utils.list2DSwapElement(genome.genomeList, (i, j), index_b)
               mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         index_a = (prng.randint(0, height), prng.randint(0, width))
         index_b = (prng.randint(0, height), prng.randint(0, width))
         utils.list2DSwapElement(genome.genomeList, index_a, index_b)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorIntegerRange(genome, **args):

   """ Simple integer range mutator for G2DList
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   range_min = genome.getParam("rangemin", constants.CDefRangeMin)
   range_max = genome.getParam("rangemax", constants.CDefRangeMax)

   if mutations < 1.0:
      mutations = 0
      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               random_int = prng.randint(range_min, range_max + 1) # HERE IT SHOULD BE INCLUSIVE
               genome.setItem(i, j, random_int)
               mutations += 1

   else:
      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())
         random_int = prng.randint(range_min, range_max + 1) # HERE IT SHOULD BE INCLUSIVE
         genome.setItem(which_y, which_x, random_int)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorIntegerGaussianGradient(genome, **args):

   """
   A gaussian mutator for G2DList of Integers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   This routine generates a gaussian value with mu=1.0 and std=0.0333 and then
   the gene is multiplied by this value. This will cause the gene to drift
   no matter how large it is.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   mu = constants.CDefGaussianGradientMU
   sigma = constants.CDefGaussianGradientSIGMA

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               final_value = int(genome[i][j] * abs(prng.normal(mu, sigma)))

               final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
               final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

               genome.setItem(i, j, final_value)
               mutations += 1
   else:

      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())

         final_value = int(genome[which_y][which_x] * abs(prng.normal(mu, sigma)))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome.setItem(which_y, which_x, final_value)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorIntegerGaussian(genome, **args):

   """
   A gaussian mutator for G2DList of Integers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   mu = genome.getParam("gauss_mu")
   sigma = genome.getParam("gauss_sigma")

   if mu is None:
      mu = constants.CDefG2DListMutIntMU

   if sigma is None:
      sigma = constants.CDefG2DListMutIntSIGMA

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               final_value = genome[i][j] + int(prng.normal(mu, sigma))

               final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
               final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

               genome.setItem(i, j, final_value)
               mutations += 1
   else:

      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())

         final_value = genome[which_y][which_x] + int(prng.normal(mu, sigma))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome.setItem(which_y, which_x, final_value)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorAllele(genome, **args):

   """ The mutator of G2DList, Allele Mutator
   To use this mutator, you must specify the *allele* genome parameter with the
   :class:`GAllele.GAlleles` instance.
   .. warning:: the :class:`GAllele.GAlleles` instance must have the homogeneous flag enabled
   """

   if args["pmut"] <= 0.0:
      return 0
   listSize = genome.getHeight() * genome.getWidth() - 1
   mutations = args["pmut"] * (listSize + 1)

   allele = genome.getParam("allele", None)
   if allele is None:
      utils.raiseException("to use the G2DListMutatorAllele, you must specify the 'allele' parameter", TypeError)

   if not allele.homogeneous:
      utils.raiseException("to use the G2DListMutatorAllele, the 'allele' must be homogeneous")

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               new_val = allele[0].getRandomAllele()
               genome.setItem(i, j, new_val)
               mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getHeight())
         which_y = prng.randint(0, genome.getWidth())

         new_val = allele[0].getRandomAllele()
         genome.setItem(which_x, which_y, new_val)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorRealGaussian(genome, **args):

   """ A gaussian mutator for G2DList of Real
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   mu = genome.getParam("gauss_mu")
   sigma = genome.getParam("gauss_sigma")

   if mu is None:
      mu = constants.CDefG2DListMutRealMU

   if sigma is None:
      sigma = constants.CDefG2DListMutRealSIGMA

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               final_value = genome[i][j] + prng.normal(mu, sigma)

               final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
               final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

               genome.setItem(i, j, final_value)
               mutations += 1
   else:

      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())

         final_value = genome[which_y][which_x] + prng.normal(mu, sigma)

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome.setItem(which_y, which_x, final_value)

   return int(mutations)

# -----------------------------------------------------------------

def G2DListMutatorRealGaussianGradient(genome, **args):

   """ A gaussian gradient mutator for G2DList of Real
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   The difference is that this multiplies the gene by gauss(1.0, 0.0333), allowing
   for a smooth gradient drift about the value.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   mu = constants.CDefGaussianGradientMU
   sigma = constants.CDefGaussianGradientSIGMA

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               final_value = genome[i][j] * abs(prng.normal(mu, sigma))

               final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
               final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

               genome.setItem(i, j, final_value)
               mutations += 1
   else:

      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())

         final_value = genome[which_y][which_x] * abs(prng.normal(mu, sigma))

         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))

         genome.setItem(which_y, which_x, final_value)

   return int(mutations)

# -----------------------------------------------------------------

def G2DBinaryStringMutatorSwap(genome, **args):

   """ The mutator of G2DBinaryString, Swap Mutator
   .. versionadded:: 0.6
      The *G2DBinaryStringMutatorSwap* function
   """

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   if mutations < 1.0:
      mutations = 0
      for i in xrange(height):
         for j in xrange(width):
            if utils.randomFlipCoin(args["pmut"]):
               index_b = (prng.randint(0, height), prng.randint(0, width))
               utils.list2DSwapElement(genome.genomeString, (i, j), index_b)
               mutations += 1
   else:
      for it in xrange(int(round(mutations))):
         index_a = (prng.randint(0, height), prng.randint(0, width))
         index_b = (prng.randint(0, height), prng.randint(0, width))
         utils.list2DSwapElement(genome.genomeString, index_a, index_b)

   return int(mutations)

# -----------------------------------------------------------------

def G2DBinaryStringMutatorFlip(genome, **args):

   """ A flip mutator for G2DBinaryString
   .. versionadded:: 0.6
      The *G2DBinaryStringMutatorFlip* function
   """

   if args["pmut"] <= 0.0:
      return 0
   height, width = genome.getSize()
   elements = height * width

   mutations = args["pmut"] * elements

   if mutations < 1.0:
      mutations = 0

      for i in xrange(genome.getHeight()):
         for j in xrange(genome.getWidth()):
            if utils.randomFlipCoin(args["pmut"]):
               if genome[i][j] == 0:
                  genome.setItem(i, j, 1)
               else:
                  genome.setItem(i, j, 0)
               mutations += 1
   else:

      for it in xrange(int(round(mutations))):
         which_x = prng.randint(0, genome.getWidth())
         which_y = prng.randint(0, genome.getHeight())

         if genome[which_y][which_x] == 0:
            genome.setItem(which_y, which_x, 1)
         else:
            genome.setItem(which_y, which_x, 0)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeMutatorSwap(genome, **args):

   """ The mutator of GTree, Swap Mutator
   .. versionadded:: 0.6
      The *GTreeMutatorSwap* function
   """

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            nodeOne = genome.getRandomNode()
            nodeTwo = genome.getRandomNode()
            nodeOne.swapNodeData(nodeTwo)
   else:
      for it in xrange(int(round(mutations))):
         nodeOne = genome.getRandomNode()
         nodeTwo = genome.getRandomNode()
         nodeOne.swapNodeData(nodeTwo)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeMutatorIntegerRange(genome, **args):

   """ The mutator of GTree, Integer Range Mutator
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   .. versionadded:: 0.6
      The *GTreeMutatorIntegerRange* function
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements

   range_min = genome.getParam("rangemin", constants.CDefRangeMin)
   range_max = genome.getParam("rangemax", constants.CDefRangeMax)

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            rand_node = genome.getRandomNode()
            random_int = prng.randint(range_min, range_max + 1) # HERE IT SHOULD BE INCLUSIVE
            rand_node.setData(random_int)

   else:
      for it in xrange(int(round(mutations))):
         rand_node = genome.getRandomNode()
         random_int = prng.randint(range_min, range_max + 1) # HERE IT SHOULD BE INCLUSIVE
         rand_node.setData(random_int)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeMutatorRealRange(genome, **args):

   """ The mutator of GTree, Real Range Mutator
   Accepts the *rangemin* and *rangemax* genome parameters, both optional.
   .. versionadded:: 0.6
      The *GTreeMutatorRealRange* function
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements

   range_min = genome.getParam("rangemin", constants.CDefRangeMin)
   range_max = genome.getParam("rangemax", constants.CDefRangeMax)

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            rand_node = genome.getRandomNode()
            random_real = prng.uniform(range_min, range_max)
            rand_node.setData(random_real)

   else:
      for it in xrange(int(round(mutations))):
         rand_node = genome.getRandomNode()
         random_real = prng.uniform(range_min, range_max)
         rand_node.setData(random_real)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeMutatorIntegerGaussian(genome, **args):

   """ A gaussian mutator for GTree of Integers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements

   mu = genome.getParam("gauss_mu", constants.CDefG1DListMutIntMU)
   sigma = genome.getParam("gauss_sigma", constants.CDefG1DListMutIntSIGMA)

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            rand_node = genome.getRandomNode()
            final_value = rand_node.getData() + int(prng.normal(mu, sigma))
            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))
            rand_node.setData(final_value)
   else:
      for it in xrange(int(round(mutations))):
         rand_node = genome.getRandomNode()
         final_value = rand_node.getData() + int(prng.normal(mu, sigma))
         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))
         rand_node.setData(final_value)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeMutatorRealGaussian(genome, **args):

   """ A gaussian mutator for GTree of Real numbers
   Accepts the *rangemin* and *rangemax* genome parameters, both optional. Also
   accepts the parameter *gauss_mu* and the *gauss_sigma* which respectively
   represents the mean and the std. dev. of the random distribution.
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements

   mu = genome.getParam("gauss_mu", constants.CDefG1DListMutRealMU)
   sigma = genome.getParam("gauss_sigma", constants.CDefG1DListMutRealSIGMA)

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            rand_node = genome.getRandomNode()
            final_value = rand_node.getData() + prng.normal(mu, sigma)
            final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
            final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))
            rand_node.setData(final_value)
   else:
      for it in xrange(int(round(mutations))):
         rand_node = genome.getRandomNode()
         final_value = rand_node.getData() + prng.normal(mu, sigma)
         final_value = min(final_value, genome.getParam("rangemax", constants.CDefRangeMax))
         final_value = max(final_value, genome.getParam("rangemin", constants.CDefRangeMin))
         rand_node.setData(final_value)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeGPMutatorOperation(genome, **args):

   """ The mutator of GTreeGP, Operation Mutator
   .. versionadded:: 0.6
      The *GTreeGPMutatorOperation* function
   """

   from . import constants

   if args["pmut"] <= 0.0:
      return 0
   elements = len(genome)
   mutations = args["pmut"] * elements
   ga_engine = args["ga_engine"]

   gp_terminals = ga_engine.getParam("gp_terminals")
   assert gp_terminals is not None

   gp_function_set = ga_engine.getParam("gp_function_set")
   assert gp_function_set is not None

   if mutations < 1.0:
      mutations = 0
      for i in xrange(len(genome)):
         if utils.randomFlipCoin(args["pmut"]):
            mutations += 1
            rand_node = genome.getRandomNode()
            assert rand_node is not None
            if rand_node.getType() == constants.nodeType["TERMINAL"]:
               term_operator = prng.choice(gp_terminals)
            else:
               op_len = gp_function_set[rand_node.getData()]
               fun_candidates = []
               for o, l in gp_function_set.items():
                  if l == op_len:
                     fun_candidates.append(o)

               if len(fun_candidates) <= 0:
                  continue

               term_operator = prng.choice(fun_candidates)
            rand_node.setData(term_operator)
   else:
      for it in xrange(int(round(mutations))):
         rand_node = genome.getRandomNode()
         assert rand_node is not None
         if rand_node.getType() == constants.nodeType["TERMINAL"]:
            term_operator = prng.choice(gp_terminals)
         else:
            op_len = gp_function_set[rand_node.getData()]
            fun_candidates = []
            for o, l in gp_function_set.items():
               if l == op_len:
                  fun_candidates.append(o)

            if len(fun_candidates) <= 0:
               continue

            term_operator = prng.choice(fun_candidates)
         rand_node.setData(term_operator)

   return int(mutations)

# -----------------------------------------------------------------

def GTreeGPMutatorSubtree(genome, **args):

   """ The mutator of GTreeGP, Subtree Mutator
   This mutator will recreate random subtree of the tree using the grow algorithm.
   .. versionadded:: 0.6
      The *GTreeGPMutatorSubtree* function
   """

   if args["pmut"] <= 0.0:
      return 0
   ga_engine = args["ga_engine"]
   max_depth = genome.getParam("max_depth", None)
   mutations = 0

   if max_depth is None:
      utils.raiseException("You must specify the max_depth genome parameter !", ValueError)

   if max_depth < 0:
      utils.raiseException("The max_depth must be >= 1, if you want to use GTreeGPMutatorSubtree crossover !", ValueError)

   branch_list = genome.nodes_branch
   elements = len(branch_list)

   for i in xrange(elements):

      node = branch_list[i]
      assert node is not None

      if utils.randomFlipCoin(args["pmut"]):
         depth = genome.getNodeDepth(node)
         mutations += 1

         root_subtree = tree.buildGTreeGPGrow(ga_engine, 0, max_depth - depth)
         node_parent = node.getParent()

         if node_parent is None:
            genome.setRoot(root_subtree)
            genome.processNodes()
            return mutations
         else:
            root_subtree.setParent(node_parent)
            node_parent.replaceChild(node, root_subtree)
         genome.processNodes()

   return int(mutations)

# -----------------------------------------------------------------
