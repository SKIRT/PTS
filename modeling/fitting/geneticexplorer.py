#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.geneticexplorer Contains the GeneticParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .explorer import ParameterExplorer
from .genetic import GeneticAlgorithm
from ...core.tools.logging import log

# -----------------------------------------------------------------

class GeneticParameterExplorer(ParameterExplorer):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(GeneticParameterExplorer, self).__init__(config)

        # The genetic algorithm
        self.ga = None

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parameter ranges ...")

        # Get the current values in the ski file prepared by InputInitializer
        young_luminosity_guess, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
        ionizing_luminosity_guess, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
        dust_mass_guess = self.ski.get_dust_component_mass(0)

        # Inform the user
        log.info("Determining every possible combination of parameter values ...")

        # Determine the ranges
        #fuv_young_range = self.young_luminosity_range(young_luminosity_guess)
        #fuv_ionizing_range = self.ionizing_luminosity_range(ionizing_luminosity_guess)
        #dust_mass_range = self.dust_mass_range(dust_mass_guess)

        # Create the genetic algorithm
        self.ga = GeneticAlgorithm()





    # -----------------------------------------------------------------

# -----------------------------------------------------------------
