#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.extinction Contains the GalacticExtinction class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astroquery.irsa_dust import IrsaDust

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.tools import tables
from ...core.tools import types
from ...core.basics.curve import FilterCurve
from ...core.filter.filter import parse_filter
from ..tools import wavelengths

# -----------------------------------------------------------------

# The path to the PTS user/extinction directory
extinction_path = fs.join(introspection.pts_user_dir, "extinction")

# -----------------------------------------------------------------

irsa_names = {"SDSS u": "SDSS u",
              "SDSS g": "SDSS g",
              "SDSS r": "SDSS r",
              "SDSS i": "SDSS i",
              "SDSS z": "SDSS z",
              "2MASS J": "2MASS J",
              "2MASS H": "2MASS H",
              "2MASS Ks": "2MASS Ks",
              "IRAC I1": "IRAC-1",
              "IRAC I2": "IRAC-2",
              "IRAC I3": "IRAC-3",
              "IRAC I4": "IRAC-4",
              "WISE W1": "WISE-1",
              "WISE W2": "WISE-2"}

# -----------------------------------------------------------------

# Determine the path to the Extinction data directory
extinction_data_path = fs.join(introspection.pts_dat_dir("magic"), "Extinction")

# Determine thbe path to the RV=3.1 curve file
rv31_path = fs.join(extinction_data_path, "al_av3_1.dat")

# -----------------------------------------------------------------

def calculate_factor(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # lees de RV=3.1 curve in van Sinopsis
    #file_RV3_1 = '/Users/idelooze/OtherProjects/SWIFT/SINOPSISv1.6.1/data/extinction_curves/al_av3_1.dat'
    #readcol, file_RV3_1, wave_RV3_1, ext_RV3_1, format = 'F,F'

    # interpol:
    filtFUV2 = INTERPOL(filtFUV, waveFUV, wave_RV3_1)

    #j = where(filtFUV2 LT 0.)
    #filtFUV2[j] = 0.

    sumFUV = 0
    sum2FUV = 0

    nwave = n_elements(wave_RV3_1)

    #for i=0, nwave-2 do begin

        #FUV
        sumFUV = sumFUV + (ext_RV3_1[i]*filtFUV2[i]*(wave_RV3_1[i+1]-wave_RV3_1[i]))
        sum2FUV = sum2FUV + (filtFUV2[i]*(wave_RV3_1[i+1]-wave_RV3_1[i]))

    #print,sumFUV,sum2FUV

    #print, 'FUV'
    #print, sumFUV / sum2FUV


# -----------------------------------------------------------------

class IRSAGalacticExtinction(object):

    """
    This class ...
    """

    def __init__(self, coordinate_or_galaxy):

        """
        This function ...
        :param coordinate_or_galaxy:
        :return:
        """

        # Check whether the user/extinction directory exists
        if not fs.is_directory(extinction_path): fs.create_directory(extinction_path)

        # Get query
        if types.is_string_type(coordinate_or_galaxy): query = coordinate_or_galaxy
        else: query = coordinate_or_galaxy.to_astropy()

        # Determine the path to the local extinction table
        path = fs.join(extinction_path, str(query))

        # Check if the local file exists
        if not fs.is_file(path):

            # Get the extinction table from IRSA
            self.table = IrsaDust.get_extinction_table(query)

            # Save the table
            tables.write(self.table, path)

        # Load the table
        else: self.table = tables.from_file(path)

    # -----------------------------------------------------------------

    def extinction_for_filter_name(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # GALEX bands
        if "GALEX" in filter_name or "UVOT" in filter_name:

            # Get the A(V) / E(B-V) ratio
            v_band_index = tables.find_index(self.table, "CTIO V")
            av_ebv_ratio = self.table["A_over_E_B_V_SandF"][v_band_index]

            # Get the attenuation of the V band A(V)
            attenuation_v = self.table["A_SandF"][v_band_index]

            # Determine the factor
            if "NUV" in filter_name: factor = 8.0
            elif "FUV" in filter_name: factor = 7.9
            elif "W2" in filter_name: factor = 8.81867
            elif "M2" in filter_name: factor = 9.28435
            elif "W1" in filter_name: factor = 6.59213
            else: raise ValueError("Unsure which GALEX or Swift UVOT band this is: " + filter_name)

            # Calculate the attenuation
            attenuation = factor * attenuation_v / av_ebv_ratio

        # Fill in the Ha attenuation manually
        elif "Halpha" in filter_name or "656_1" in filter_name: attenuation = 0.174

        # Other bands for which attenuation is listed by IRSA
        elif filter_name in irsa_names:

            irsa_name = irsa_names[filter_name]

            # Find the index of the corresponding table entry
            index = tables.find_index(self.table, irsa_name)

            # Get the attenuation
            attenuation = self.table["A_SandF"][index]

        # All other bands: set attenuation to zero
        #else: #attenuation = 0.0
        # NO: RAISE AND ERROR
        else: raise ValueError("Cannot determine the extinction for the '"+  filter_name + "' filter")

        # Return the galactic attenuation coefficient
        return attenuation

    # -----------------------------------------------------------------

    def extinction_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # No extinction for this filter
        if not wavelengths.wavelength_in_regimes(fltr.wavelength, ["FUV-NIR"]): return 0.0

        filter_name = str(fltr)
        return self.extinction_for_filter_name(filter_name)

    # -----------------------------------------------------------------

    def extinction_for_filters(self, filters):

        """
        This function ...
        :param filters:
        :return:
        """

        extinctions = []

        # Loop over the filters
        for fltr in filters:

            # Parse the filter
            if types.is_string_type(fltr): fltr = parse_filter(fltr)

            # Get the extinction
            extinction = self.extinction_for_filter(fltr)

            # Add to the list
            extinctions.append(extinction)

        # Return the list
        return extinctions

    # -----------------------------------------------------------------

    def extinction_curve(self, filters, ignore_errors=False):

        """
        This function ...
        :param filters:
        :param ignore_errors:
        :return:
        """

        kwargs = dict()
        kwargs["y_name"] = "Extinction"
        kwargs["y_description"] = "Extinction value"
        kwargs["y_unit"] = None

        curve = FilterCurve(**kwargs)

        # Loop over the filter
        for fltr in filters:

            # Parse the filter
            if types.is_string_type(fltr): fltr = parse_filter(fltr)

            # Get the extinction
            try:
                extinction = self.extinction_for_filter(fltr)
                # Add extinction point
                curve.add_point(fltr, extinction)

            except ValueError:
                if ignore_errors: pass
                else: raise ValueError("Could not determine the extinction for the " + str(fltr) + "filter")

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    def correct_sed_for_extinction(self, sed):

        """
        This function ...
        :return:
        """

        # Loop over all entries
        for i in range(len(sed)):

            # Get the filter name
            instrument = sed["Instrument"][i]
            band = sed["Band"][i]
            filter_name = instrument + " " + band

            # Get the extinction
            extinction = self.extinction_for_filter_name(filter_name)

            # Correct the flux value (and error) for galactic extinction
            correction_factor = 10 ** (0.4 * extinction)
            sed["Flux"][i] *= correction_factor
            sed["Error-"][i] *= correction_factor
            sed["Error+"][i] *= correction_factor

# -----------------------------------------------------------------
