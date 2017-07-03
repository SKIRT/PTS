#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.genomes.nucleobases

# -----------------------------------------------------------------

# Import other evolve modules
from ..core.genome import GenomeBase, G1DBase
from ..core import constants

# -----------------------------------------------------------------

class NucleoBases(G1DBase):

    """ 
    NucleoBases genome
    """

    __slots__ = ["nbases"]

    # -----------------------------------------------------------------

    def __init__(self, nbases):

        """
        The initializator of the NucleoBases genome representation
        """

        # Call the constructor of the base class
        super(NucleoBases, self).__init__(nbases)

        # Set nbases
        self.nbases = nbases

        # Set function slots
        self.initializator.set(constants.CDefG1DBinaryStringInit)
        self.mutator.set(constants.CDefG1DBinaryStringMutator)
        self.crossover.set(constants.CDefG1DBinaryStringCrossover)

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        Return a string representation of the genome
        """

        ret = GenomeBase.__repr__(self)
        ret += "- G1DBinaryString\n"
        ret += "\tNumber of bases:\t %s\n" % (self.getListSize(),)
        ret += "\tBases:\t\t" + "".join(self.genomeList) + "\n\n"
        return ret

# -----------------------------------------------------------------
