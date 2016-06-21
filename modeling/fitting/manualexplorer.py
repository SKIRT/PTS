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
from .explorer import ParameterExplorer
from ...core.tools.logging import log
from ...magic.animation.scatter import ScatterAnimation
from ...core.tools import filesystem as fs
from ...core.tools import time

# -----------------------------------------------------------------

class ManualParameterExplorer(ParameterExplorer):
    
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
        super(ManualParameterExplorer, self).__init__(config)

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
        fuv_young_range = self.young_luminosity_range(young_luminosity_guess)
        fuv_ionizing_range = self.ionizing_luminosity_range(ionizing_luminosity_guess)
        dust_mass_range = self.dust_mass_range(dust_mass_guess)

        # Make animation
        if self.config.visualise:
            animation = ScatterAnimation([fuv_young_range[0], fuv_young_range[-1]], [fuv_ionizing_range[0], fuv_ionizing_range[-1]], [dust_mass_range[0].to("Msun").value, dust_mass_range[-1].to("Msun").value])
            animation.x_label = "FUV young"
            animation.y_label = "FUV ionizing"
            animation.z_label = "Dust mass"
            animation.density = False
        else: animation = None

        # Loop over the different values of the young stellar luminosity
        for young_luminosity in fuv_young_range:

            # Loop over the different values of the ionizing stellar luminosity
            for ionizing_luminosity in fuv_ionizing_range:

                # Loop over the different values of the dust mass
                for dust_mass in dust_mass_range:

                    # Add the parameter values to the dictionary
                    self.parameters["FUV young"].append(young_luminosity)
                    self.parameters["FUV ionizing"].append(ionizing_luminosity)
                    self.parameters["Dust mass"].append(dust_mass)

                    # Add a point (and thus a frame) to the animation
                    if animation is not None: animation.add_point(young_luminosity, ionizing_luminosity, dust_mass.to("Msun").value)

        # Save the animation
        if animation is not None:
            path = fs.join(self.visualisation_path, time.unique_name("initialparameterexploration") + ".gif")
            animation.save(path)

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
