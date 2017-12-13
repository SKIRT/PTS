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
from ..core.mappings import Mappings
from ...core.filter.filter import parse_filter
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class Model(object):
    
    """
    This class...
    """

    def __init__(self, simulation_name=None, chi_squared=None, parameter_values=None):

        """
        The constructor ...
        :param simulation_name:
        :param chi_squared:
        :param parameter_values:
        :return:
        """

        # The attributes
        self.simulation_name = simulation_name
        self.chi_squared = chi_squared
        self.parameter_values = parameter_values

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values.keys() if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def i1_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter("IRAC I1")

    # -----------------------------------------------------------------

    @property
    def i1_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.i1_filter.wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):

        """
        This function ...
        :return:
        """

        return parse_filter("GALEX FUV")

    # -----------------------------------------------------------------

    @property
    def fuv_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.fuv_filter.wavelength

    # -----------------------------------------------------------------

    @property
    def metallicity(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values["metallicity"] if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @property
    def sfr_compactness(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values["sfr_compactness"] if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @property
    def sfr_pressure(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values["sfr_pressure"] if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @property
    def sfr_covering_factor(self):

        """
        This function ...
        :return:
        """

        return self.parameter_values["sfr_covering"] if self.parameter_values is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr(self):

        """
        This function derives the SFR (in Msun / year) from the FUV luminosity of the model and the intrinsic MAPPINGS SED
        :return:
        """

        # Get the FUV luminosity of the ionizing stellar component
        fuv_luminosity = self.parameter_values["fuv_ionizing"]

        # Get the SFR
        return Mappings.sfr_for_luminosity(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, fuv_luminosity, self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings(self):

        """
        This function ...
        :return:
        """

        # Create the MAPPINGS template and return it
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor, self.sfr)

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings(self):

        """
        This function ...
        :return:
        """

        # Create the MAPPINGS template
        return Mappings(self.metallicity, self.sfr_compactness, self.sfr_pressure, self.sfr_covering_factor)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.mappings.dust_mass

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_mass(self):

        """
        This function ...
        :return:
        """

        return self.mappings.stellar_mass

# -----------------------------------------------------------------
