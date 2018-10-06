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

# Import the relevant PTS classes and modules
from ..tools import tables, introspection, arrays
from ..tools import filesystem as fs
from ..basics.curve import Curve
from ..units.parsing import parse_unit as u

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
        :param total: OBSERVED
        :param transparent: INTRINSIC
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

    def closest_wavelength_index(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        return arrays.find_closest_index(self["Wavelength"], wavelength, array_unit=self["Wavelength"].unit)

    # -----------------------------------------------------------------

    def closest_wavelength(self, wavelength, return_index=False):

        """
        This function ...
        :param wavelength:
        :param return_index:
        :return:
        """

        index = self.closest_wavelength_index(wavelength)
        wavelength = self["Wavelength"][index] * self["Wavelength"].unit
        if return_index: return wavelength, index
        else: return wavelength

    # -----------------------------------------------------------------

    def normalize_at(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        attenuation_wavelength = self.attenuation_at(wavelength)
        self["Attenuation"] *= (value / attenuation_wavelength)

    # -----------------------------------------------------------------

    @property
    def peak_attenuation(self):
        return np.max(self["Attenuation"])

    # -----------------------------------------------------------------

    def normalize_at_peak(self, value=1.):

        """
        This function ...
        :param value:
        :return:
        """

        self["Attenuation"] *= (value / self.peak_attenuation)

# -----------------------------------------------------------------

def load_mappings_attenuation_data():

    """
    This function ...
    :return:
    """

    # The filepath
    path = fs.join(attenuation_data_path, "AttenuationLawMAPPINGS.dat")

    # wl in micron from long to short wl.
    # ABS attenuations (see header of data file)
    wavelengths, abs_attenuations = np.loadtxt(path, unpack=True)

    # CREATE A TABLE SO WE CAN EASILY SORT THE COLUMNS FOR INCREASING WAVELENGTH
    names = ["Wavelength", "ABS attenuation"]
    # Create the table
    abs_table = tables.new([wavelengths, abs_attenuations], names)
    abs_table["Wavelength"].unit = u("micron")
    # Sort the table on wavelength
    abs_table.sort("Wavelength")

    # Return
    return abs_table

# -----------------------------------------------------------------

class MappingsAttenuationCurve(AttenuationCurve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Check whether constructor is called by Astropy
        if "from_astropy" in kwargs: from_astropy = kwargs.pop("from_astropy")
        elif "wavelength" in kwargs: from_astropy = False
        else: from_astropy = True

        # Not called by Astropy
        if not from_astropy:

            # Get attenuation at a
            attenuation = kwargs.pop("attenuation", None)
            wavelength = kwargs.pop("wavelength", None)
            if attenuation is not None and wavelength is None: raise ValueError("If attenuation value is specified, wavelength must be specified")
            if wavelength is not None and attenuation is None: raise ValueError("If wavelength is specified, attenuation value must be specified")

            # Load the data
            abs_table = load_mappings_attenuation_data()
            wavelengths = np.array(list(abs_table["Wavelength"]))
            abs_attenuations = np.array(list(abs_table["ABS attenuation"]))

            #print(wavelengths)
            #print(abs_attenuations)
            #print(len(wavelengths), len(abs_attenuations))

            # Find the ABS attenuation at the specified wavelength
            if attenuation is not None and wavelength is not None:

                # Normalize
                interpolated = interpolate.interp1d(wavelengths, abs_attenuations, kind='linear')
                abs_wavelength = interpolated(wavelength.to("micron").value)

                # 'Real' attenuations
                attenuations = abs_attenuations / abs_wavelength * attenuation

            # Just use the ABS attenuations
            else: attenuations = abs_attenuations

            #print(len(wavelengths), wavelengths)
            #print(len(attenuations), attenuations)

            # Set arguments
            kwargs["wavelengths"] = wavelengths
            kwargs["attenuations"] = attenuations

        # Call the constructor of the base class
        super(MappingsAttenuationCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------
