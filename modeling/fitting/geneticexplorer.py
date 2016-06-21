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

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ParameterExplorer instance
        explorer = cls(arguments.config)

        # Set the modeling path
        explorer.config.path = arguments.path

        # Set the remote host IDs
        if arguments.remote is not None: explorer.config.remotes = arguments.remote

        # Set options for the young stellar population
        if arguments.young_nvalues is not None: explorer.config.young_stars.nvalues = arguments.young_nvalues
        if arguments.young_range is not None:
            explorer.config.young_stars.rel_min = arguments.young_range[0]
            explorer.config.young_stars.rel_max = arguments.young_range[1]
        if arguments.young_log: explorer.config.young_stars.scale = "log"
        else: explorer.config.young_stars.scale = "linear"

        # Set options for the ionizing stellar population
        if arguments.ionizing_nvalues is not None: explorer.config.ionizing_stars.nvalues = arguments.ionizing_nvalues
        if arguments.ionizing_range is not None:
            explorer.config.ionizing_stars.rel_min = arguments.ionizing_range[0]
            explorer.config.ionizing_stars.rel_max = arguments.ionizing_range[1]
        if arguments.ionizing_log: explorer.config.ionizing_stars = "log"
        else: explorer.config.ionizing_stars.scale = "linear"

        # Set options for the dust component
        if arguments.dust_nvalues is not None: explorer.config.dust.nvalues = arguments.dust_nvalues
        if arguments.dust_range is not None:
            explorer.config.dust.rel_min = arguments.dust_range[0]
            explorer.config.dust.rel_max = arguments.dust_range[1]
        if arguments.dust_log: explorer.config.dust.scale = "log"
        else: explorer.config.dust.scale = "linear"

        # Make visualisations
        explorer.config.visualise = arguments.visualise

        # Return the new instance
        return explorer

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
