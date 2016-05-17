#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialparameterexplorer Contains the InitialParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .parameterexploration import ParameterExplorer
from ...core.tools.logging import log

# -----------------------------------------------------------------

class InitialParameterExplorer(ParameterExplorer):
    
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
        super(InitialParameterExplorer, self).__init__(config)

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
        if arguments.remotes is not None: explorer.config.remotes = arguments.remotes

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

        # Return the new instance
        return explorer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the current parameter table
        self.load_table()

        # 3. Load the ski file
        self.load_ski()

        # 4. Set the parameter combinations
        self.set_parameters()

        # 5. Launch the simulations for different parameter values
        self.simulate()

        # 6. Writing
        self.write()

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

        # Loop over the different values of the young stellar luminosity
        for young_luminosity in self.young_luminosity_range(young_luminosity_guess):

            # Loop over the different values of the ionizing stellar luminosity
            for ionizing_luminosity in self.ionizing_luminosity_range(ionizing_luminosity_guess):

                # Loop over the different values of the dust mass
                for dust_mass in self.dust_mass_range(dust_mass_guess):

                    # Add the parameter values to the dictionary
                    self.parameters["FUV young"].append(young_luminosity)
                    self.parameters["FUV ionizing"].append(ionizing_luminosity)
                    self.parameters["Dust mass"].append(dust_mass)

    # -----------------------------------------------------------------

    def young_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the FUV luminosity of the young stars ...")

        # Set the range of the FUV luminosity of the young stellar population
        min_value = self.config.young_stars.rel_min * luminosity
        max_value = self.config.young_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.young_stars.scale == "linear":
            young_luminosities = np.linspace(min_value, max_value, num=self.config.young_stars.nvalues, endpoint=True)
        elif self.config.young_stars.scale == "logarithmic":
            young_luminosities = np.logspace(min_value, max_value, num=self.config.young_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the young stellar luminosity values")

        # Return the range of FUV luminosities of the young stellar component
        return young_luminosities

    # -----------------------------------------------------------------

    def ionizing_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the FUV luminosity of the ionizing stars ...")

        # Determine the minimum and maximum FUV luminosity of the ionizing stellar population
        min_value = self.config.ionizing_stars.rel_min * luminosity
        max_value = self.config.ionizing_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.ionizing_stars.scale == "linear":
            ionizing_luminosities = np.linspace(min_value, max_value, num=self.config.ionizing_stars.nvalues, endpoint=True)
        elif self.config.ionizing_stars.scale == "log":
            ionizing_luminosities = np.logspace(min_value, max_value, num=self.config.ionizing_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the ionizing stellar luminosity values")

        # Return the range of FUV luminosities of the ionizing stellar component
        return ionizing_luminosities

    # -----------------------------------------------------------------

    def dust_mass_range(self, mass):

        """
        This function ...
        :param mass:
        :return:
        """

        # Inform the user
        log.info("Setting the range for the dust mass ...")

        # Set the dust mass range
        min_value = self.config.dust.rel_min * mass
        max_value = self.config.dust.rel_max * mass

        # Create a linear or logarithmic range of dust masses
        if self.config.dust.scale == "linear": dust_masses = np.linspace(min_value, max_value, num=self.config.dust.nvalues, endpoint=True)
        elif self.config.dust.scale == "log": dust_masses = np.logspace(min_value, max_value, num=self.config.dust.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the dust mass values")

        # Return the range of dust masses
        return dust_masses

# -----------------------------------------------------------------
