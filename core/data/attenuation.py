#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.attenuation Contains the AttenuationCurve classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import interpolate

# Import astronomical modules
from astropy.units import Unit, spectral
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables, introspection, arrays
from ...core.tools import filesystem as fs
from ..basics.curve import Curve

# -----------------------------------------------------------------

attenuation_data_path = fs.join(introspection.pts_dat_dir("core"), "attenuation")

# -----------------------------------------------------------------

class AttenuationCurve(Curve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Units
        x_unit = "micron"

        # Names
        x_name = "Wavelength"
        y_name = "Attenuation"

        # Descriptions
        x_description = "Wavelength"
        y_description = "Attenuation"

        kwargs["x_unit"] = x_unit
        kwargs["y_unit"] = None
        kwargs["x_name"] = x_name
        kwargs["y_name"] = y_name
        kwargs["x_description"] = x_description
        kwargs["y_description"] = y_description

        # If data is passed
        if "wavelengths" in kwargs and "attenuations" in kwargs:

            wavelengths = kwargs.pop("wavelengths")
            attenuations = kwargs.pop("attenuations")

        else: wavelengths = attenuations = None

        # Call the constructor of the base class
        super(AttenuationCurve, self).__init__(*args, **kwargs)

        # Add the data
        if wavelengths is not None:
            for index in range(len(wavelengths)): self.add_row([wavelengths[index], attenuations[index]])

    # -----------------------------------------------------------------

    @classmethod
    def from_seds(cls, total, transparent):

        """
        This function ...
        :param total:
        :param transparent:
        :return:
        """

        # Get the wavelengths
        wavelengths = total.wavelengths(unit="micron", add_unit=False)

        # Get the total and transparent fluxes
        total_fluxes = total.photometry(asarray=True)
        transparent_fluxes = transparent.photometry(asarray=True)

        # Calculate the attenuations
        attenuations = -2.5 * np.log10(total_fluxes / transparent_fluxes)

        # Create a new AttenuationCurve instance
        return cls(wavelengths, attenuations)

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self["Wavelength"], unit=unit, array_unit=self.column_unit("Wavelength"))
        else: return arrays.array_as_list(self["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Wavelength"))

    # -----------------------------------------------------------------

    def attenuations(self, asarray=False):

        """
        This function ...
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["Attenuation"])
        else: return arrays.array_as_list(self["Attenuation"])

    # -----------------------------------------------------------------

    def attenuation_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        interpolated = interpolate.interp1d(self.wavelengths(unit="micron", asarray=True), self.attenuations(asarray=True), kind='linear')
        return interpolated(wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def normalize_at(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        attenuation_wavelength = self.attenuation_at(wavelength)
        self["Attenuation"] /= attenuation_wavelength * value

# -----------------------------------------------------------------


# -----------------------------------------------------------------
