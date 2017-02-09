#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.step Contains the StepOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from textwrap import wrap

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ..genomes.list1d import G1DList
from ..genomes.list2d import G2DList
from ..genomes.binarystring1d import G1DBinaryString
from ..genomes.binarystring2d import G2DBinaryString
from ...core.tools.logging import log
from ..core.initializators import G1DListInitializatorReal, G1DListInitializatorInteger
from ..core.mutators import G1DListMutatorIntegerRange, G1DListMutatorIntegerGaussian, G1DListMutatorIntegerBinary, G1DListMutatorRealGaussian, G1DListMutatorRealRange
from ..core.engine import GeneticEngine, RawScoreCriteria
from ..core import constants
from ...core.basics.range import RealRange, IntegerRange
from ...core.tools import formatting as fmt
from ...core.tools import stringify

# -----------------------------------------------------------------

class StepOptimizer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(Optimizer, self).__init__(config)

        # The evaluating function
        self.evaluator = None

        # The initializator function
        self.initializator = None

        # The mutator function
        self.mutator = None

        # The crossover function
        self.crossover = None

        # The callback function
        self.callback = None

        # The database adapter
        self.adapter = None

        # Kwargs for the various functions
        self.evaluator_kwargs = None
        self.initializator_kwargs = None
        self.mutator_kwargs = None
        self.crossover_kwargs = None
        self.callback_kwargs = None

        # The starting genome
        self.genome = None

        # The genetic algorithm engine
        self.engine = None

        # The best individual
        self.best = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Initialize
        self.initialize()

        # 3. Evolve
        self.evolve()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Optimizer, self).setup(**kwargs)

        # Get the genome
        if "genome" in kwargs: self.genome = kwargs.pop("genome").clone()

        # Get evaluator
        if "evaluator" in kwargs: self.evaluator = kwargs.pop("evaluator")

        # Get initializator
        if "initializator" in kwargs: self.initializator = kwargs.pop("initializator")

        # Get mutator
        if "mutator" in kwargs: self.mutator = kwargs.pop("mutator")

        # Get crossover function
        if "crossover" in kwargs: self.crossover = kwargs.pop("crossover")

        # Get the callback function
        if "callback" in kwargs: self.callback = kwargs.pop("callback")

        # Get the database adapter function
        if "adapter" in kwargs: self.adapter = kwargs.pop("adapter")

        # Get kwargs for the different functions
        if "evaluator_kwargs" in kwargs: self.evaluator_kwargs = kwargs.pop("evaluator_kwargs")
        if "initializator_kwargs" in kwargs: self.initializator_kwargs = kwargs.pop("initializator_kwargs")
        if "mutator_kwargs" in kwargs: self.mutator_kwargs = kwargs.pop("mutator_kwargs")
        if "crossover_kwargs" in kwargs: self.crossover_kwargs = kwargs.pop("crossover_kwargs")
        if "callback_kwargs" in kwargs: self.callback_kwargs = kwargs.pop("callback_kwargs")

    # -----------------------------------------------------------------



# -----------------------------------------------------------------
