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
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.tools import tables

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

class GalacticExtinction(object):

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
        if isinstance(coordinate_or_galaxy, basestring): query = coordinate_or_galaxy
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
            else: raise ValueError("Unsure which GALEX or Swift UVOT band this is")

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
        else: attenuation = 0.0

        # Return the galactic attenuation coefficient
        return attenuation

    # -----------------------------------------------------------------

    def extinction_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        filter_name = str(fltr)
        return self.extinction_for_filter_name(filter_name)

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
