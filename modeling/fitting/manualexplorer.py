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
import math
import numpy as np

# Import the relevant PTS classes and modules
from .explorer import ParameterExplorer
from ...core.tools.logging import log

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
        log.info("Setting the parameter values for the different models ...")

        # Loop over the different values of the young stellar luminosity
        for young_luminosity in self.get_fuv_young_values():

            # Loop over the different values of the ionizing stellar luminosity
            for ionizing_luminosity in self.get_fuv_ionizing_values():

                # Loop over the different values of the dust mass
                for dust_mass in self.get_dust_mass_values():

                    # Add the parameter values to the dictionary
                    self.parameters["FUV young"].append(young_luminosity)
                    self.parameters["FUV ionizing"].append(ionizing_luminosity)
                    self.parameters["Dust mass"].append(dust_mass)

                    # Update the animation
                    if self.config.visualise: self.scatter_animation.add_point(young_luminosity, ionizing_luminosity, dust_mass)

    # -----------------------------------------------------------------

    def get_fuv_young_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the values for the FUV luminosity of the young stars ...")

        # Get min and max value
        min_value = self.ranges["FUV young"][0]
        max_value = self.ranges["FUV young"][1]

        # Get number of values
        nvalues = int(math.ceil(self.config.simulations**(1./3.)))

        # Create a linear or logarithmic range of values
        if self.config.young_log: values = np.logspace(min_value, max_value, num=nvalues, endpoint=True)
        else: values = np.linspace(min_value, max_value, num=nvalues, endpoint=True)

        # Return the FUV luminosities of the young stellar component
        return values

    # -----------------------------------------------------------------

    def get_fuv_ionizing_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the values for the FUV luminosity of the ionizing stars ...")

        # Get min and max value
        min_value = self.ranges["FUV ionizing"][0]
        max_value = self.ranges["FUV ionizing"][1]

        # Get number of values
        nvalues = int(math.ceil(self.config.simulations ** (1. / 3.)))

        # Create a linear or logarithmic range of luminosities
        if self.config.ionizing_log: values = np.logspace(min_value, max_value, num=nvalues, endpoint=True)
        else: values = np.linspace(min_value, max_value, num=nvalues, endpoint=True)

        # Return the FUV luminosities of the ionizing stellar component
        return values

    # -----------------------------------------------------------------

    def get_dust_mass_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the values for the dust mass ...")

        # Get min and max value
        min_value = self.ranges["Dust mass"][0]
        max_value = self.ranges["Dust mass"][1]

        # Get number of values
        nvalues = int(math.ceil(self.config.simulations ** (1. / 3.)))

        # Create a linear or logarithmic range of dust masses
        if self.config.dust_log: values = np.logspace(min_value, max_value, num=nvalues, endpoint=True)
        else: values = np.linspace(min_value, max_value, num=nvalues, endpoint=True)

        # Return the dust masses
        return values

# -----------------------------------------------------------------
