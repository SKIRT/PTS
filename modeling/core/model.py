#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.model Contains the Model class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..core.mappings import Mappings
from ...core.basics.filter import Filter

# -----------------------------------------------------------------

class Model(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        self.simulation_name = None
        self.chi_squared = None
        self.parameter_values = None

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values.keys() if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @property
    def sfr(self):

        """
        This function derives the SFR (in Msun / year) from the FUV luminosity of the model and the intrinsic MAPPINGS SED
        :return:
        """

        # Get the relevant parameters
        metallicity = self.parameter_values["metallicity"]
        compactness = self.parameter_values["sfr_compactness"]
        pressure = self.parameter_values["sfr_pressure"]
        covering_factor = self.parameter_values["sfr_covering"]

        # Get the FUV luminosity of the ionizing stella
        fuv_luminosity = self.parameter_values["fuv_ionizing"]

        # Get the FUV pivot wavelength
        fuv_wavelength = Filter.from_string("GALEX FUV").pivot

        # Get the SFR
        sfr = Mappings.sfr_for_luminosity(metallicity, compactness, pressure, covering_factor, fuv_luminosity, fuv_wavelength)

        # Return the SFR
        return sfr

    # -----------------------------------------------------------------

    @property
    def sfr_dust_mass(self):

        """
        This function ...
        :return:
        """

        # Get the relevant parameters
        metallicity = self.parameter_values["metallicity"]
        compactness = self.parameter_values["sfr_compactness"]
        pressure = self.parameter_values["sfr_pressure"]
        covering_factor = self.parameter_values["sfr_covering"]

        # Get the FUV luminosity of the ionizing stella
        fuv_luminosity = self.parameter_values["fuv_ionizing"]

        # Get the FUV pivot wavelength
        fuv_wavelength = Filter.from_string("GALEX FUV").pivot

        # Get the dust mass
        dust_mass = Mappings.dust_mass_for_luminosity(metallicity, compactness, pressure, covering_factor, fuv_luminosity, fuv_wavelength)

        # Return the dust mass
        return dust_mass

    # -----------------------------------------------------------------

    @property
    def mappings(self):

        """
        This function ...
        :return:
        """

        # Get the relevant parameters
        metallicity = self.parameter_values["metallicity"]
        compactness = self.parameter_values["sfr_compactness"]
        pressure = self.parameter_values["sfr_pressure"]
        covering_factor = self.parameter_values["sfr_covering"]

        # Create the MAPPINGS template
        mappings = Mappings(metallicity, compactness, pressure, covering_factor, self.sfr.to("Msun / yr").value)

        # Return the mappings template
        return mappings

# -----------------------------------------------------------------
