#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.services.attenuation Contains the GalacticAttenuation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from scipy.interpolate import interp1d
import numpy as np

# Import astronomical modules
from astroquery.irsa_dust import IrsaDust

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.tools import tables
from ...core.tools import types
from ...core.basics.curve import FilterCurve, Curve
from ...core.filter.filter import parse_filter
from ..tools import wavelengths
from ...core.basics.log import log
from ...core.tools import sequences
from ...core.filter.broad import BroadBandFilter
from ...core.filter.narrow import NarrowBandFilter
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# The path to the PTS user/attenuation directory
attenuation_path = fs.join(introspection.pts_user_dir, "attenuation")

# The path to the RFILTER table
rfilter_table_path = fs.join(attenuation_path, "rfilter.dat")

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

# Determine the path to the RV=3.1 curve file
rv31_path = fs.join(extinction_data_path, "al_av3_1.dat")

# -----------------------------------------------------------------

class RFilterTable(FilterCurve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set kwargs
        kwargs["y_name"] = "R"
        kwargs["y_description"] = "R value"
        kwargs["y_unit"] = None

        # Call the constructor of the base class
        super(RFilterTable, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

def calculate_r(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # Calcualte R
    if isinstance(fltr, BroadBandFilter): return calculate_r_broad(fltr)
    elif isinstance(fltr, NarrowBandFilter): return calculate_r_narrow(fltr)
    else: raise ValueError("Filter object not recognized")

# -----------------------------------------------------------------

def calculate_r_narrow(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    log.debug("Calculating R for the " + str(fltr) + " filter ...")

    rv31_curve = Curve.from_data_file(rv31_path, x_name="Wavelength", y_name="R")
    # To micron from Angstrom
    wavelengths = np.array(rv31_curve.x_data) * 0.0001
    rvs = rv31_curve.y_data

    # Return
    rvs_interp = interp1d(wavelengths, rvs)
    return rvs_interp(fltr.wavelength.to("micron").value)

# -----------------------------------------------------------------

def calculate_r_broad(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # lees de RV=3.1 curve in van Sinopsis
    #file_RV3_1 = '/Users/idelooze/OtherProjects/SWIFT/SINOPSISv1.6.1/data/extinction_curves/al_av3_1.dat'
    #readcol, file_RV3_1, wave_RV3_1, ext_RV3_1, format = 'F,F'

    # interpol:
    #filtFUV2 = INTERPOL(filtFUV, waveFUV, wave_RV3_1)

    log.debug("Calculating R for the " + str(fltr) + " filter ...")

    #print("filter wavelengths", fltr.wavelengths)

    #transmissions = fltr.transmissions
    transmissionfun = interp1d(fltr.wavelengths, fltr.transmissions)

    minwavelength = min(fltr.wavelengths)
    maxwavelength = max(fltr.wavelengths)

    rv31_curve = Curve.from_data_file(rv31_path, x_name="Wavelength", y_name="R")
    # To micron from Angstrom
    wavelengths = np.array(rv31_curve.x_data) * 0.0001
    rvs = rv31_curve.y_data

    mask = np.array([wavelength > maxwavelength or wavelength < minwavelength for wavelength in wavelengths])
    inverse_mask = np.logical_not(mask)
    wavelengths = wavelengths[inverse_mask]
    rvs = rvs[inverse_mask]

    #print("wavelengths", wavelengths)

    # Calculate interpolated transmissions
    interpolated = transmissionfun(wavelengths)
    for index in range(len(interpolated)):
        if interpolated[index] < 0: interpolated[index] = 0.0

    #j = where(filtFUV2 LT 0.)
    #filtFUV2[j] = 0.

    sumFUV = 0.
    sum2FUV = 0.

    #nwave = n_elements(wave_RV3_1)
    nwave = len(wavelengths)

    #for i=0, nwave-2 do begin

    for i in range(nwave-1):

        wave_delta = wavelengths[i+1] - wavelengths[i]

        #FUV
        sumFUV += rvs[i] * interpolated[i] * wave_delta
        sum2FUV += interpolated[i] * wave_delta

    #print("SUMS", sumFUV, sum2FUV)

    result = sumFUV / sum2FUV

    index = sequences.find_closest_index(wavelengths, fltr.wavelength.to("micron").value)
    lazy_result = rvs[index]

    #print("RESULT", result, lazy_result)

    #print, 'FUV'
    #print, sumFUV / sum2FUV

    # Return the result
    return result

# -----------------------------------------------------------------

class GalacticAttenuation(object):

    """
    This class ...
    """

    def __init__(self, coordinate_or_galaxy):

        """
        This function ...
        :param coordinate_or_galaxy:
        :return:
        """

        # Check whether the user/attenuation directory exists
        if not fs.is_directory(attenuation_path): fs.create_directory(attenuation_path)

        # Get query
        if types.is_string_type(coordinate_or_galaxy): query = coordinate_or_galaxy
        else: query = coordinate_or_galaxy.to_astropy()

        # Determine the path to the local extinction table
        path = fs.join(attenuation_path, str(query))

        # Check if the local file exists
        if not fs.is_file(path):

            # Get the extinction table from IRSA
            self.table = IrsaDust.get_extinction_table(query)

            # Save the table
            tables.write(self.table, path)

        # Load the table
        else: self.table = tables.from_file(path)

        # Check whether the RFILTER table is present, load
        if fs.is_file(rfilter_table_path): self.rfilter = RFilterTable.from_file(rfilter_table_path)
        else: # or create new
            self.rfilter = RFilterTable()
            self.rfilter.path = rfilter_table_path
            self.rfilter.save()

    # -----------------------------------------------------------------

    def has_rfilter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # CHeck in table
        return self.rfilter.has_filter(fltr)

    # -----------------------------------------------------------------

    def set_rfilter(self, fltr, rfilter):

        """
        This function ...
        :param fltr:
        :param rfilter:
        :return:
        """

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Add
        self.rfilter.add_point(fltr, rfilter)

        # Save the table
        self.rfilter.save()

    # -----------------------------------------------------------------

    def get_rfilter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.rfilter.value_for_filter(fltr)

    # -----------------------------------------------------------------

    def calculate_rfilter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        r = calculate_r(fltr)
        self.set_rfilter(fltr, r)
        return r

    # -----------------------------------------------------------------

    def extinction_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # No extinction for this filter
        #if not wavelengths.wavelength_in_regimes(fltr.wavelength, ["FUV-NIR"]): return 0.0
        if fltr.wavelength not in wavelengths.extinction_wavelength_range: return 0.0

        #return self.extinction_for_filter_name(filter_name)

        if self.has_rfilter(fltr): r = self.get_rfilter(fltr)
        else: r = self.calculate_rfilter(fltr)

        # reddening ratio Av / E(B-V) = 3.1

        # Calculate the attenuation
        attenuation = r * self.attenuation_v #/ self.av_ebv_ratio

        # Check with precalculated and with IRSA
        #if self.has_precalculated_extinction(fltr):
        #    print("pre", attenuation, self.precalculated_extinction(fltr), attenuation/self.precalculated_extinction(fltr))
        #if self.has_irsa_extinction(fltr):
        #    print("irsa", attenuation, self.irsa_extinction(fltr), attenuation/self.irsa_extinction(fltr))

        # Return the attenuation
        return attenuation

    # -----------------------------------------------------------------

    @lazyproperty
    def v_band_index(self):
        return tables.find_index(self.table, "CTIO V")

    # -----------------------------------------------------------------

    @lazyproperty
    def av_ebv_ratio(self):
        return self.table["A_over_E_B_V_SandF"][self.v_band_index]

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_v(self):
        return self.table["A_SandF"][self.v_band_index]

    # -----------------------------------------------------------------

    def has_precalculated_extinction(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        filter_name = str(fltr)

        if "GALEX" in filter_name or "UVOT" in filter_name:
            if "NUV" in filter_name or "FUV" in filter_name or "W2" in filter_name or "M2" in filter_name \
                or "W1" in filter_name: return True
            else: return False
        #elif "Halpha" in filter_name or "Ha" in filter_name or "656_1" in filter_name: return True # NO: only for a particular galaxy??
        else: return False

    # -----------------------------------------------------------------

    def has_irsa_extinction(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return str(fltr) in irsa_names

    # -----------------------------------------------------------------

    def precalculated_extinction(self, fltr):

        """
        :param fltr:
        :return:
        """

        filter_name = str(fltr)

        # GALEX bands
        if "GALEX" in filter_name or "UVOT" in filter_name:

            # Get the A(V) / E(B-V) ratio
            #v_band_index = tables.find_index(self.table, "CTIO V")
            #av_ebv_ratio = self.table["A_over_E_B_V_SandF"][v_band_index]

            # Get the attenuation of the V band A(V)
            #attenuation_v = self.table["A_SandF"][v_band_index]

            # Determine the factor
            if "NUV" in filter_name: factor = 8.0
            elif "FUV" in filter_name: factor = 7.9
            elif "W2" in filter_name: factor = 8.81867
            elif "M2" in filter_name: factor = 9.28435
            elif "W1" in filter_name: factor = 6.59213
            else: raise ValueError("Unsure which GALEX or Swift UVOT band this is: " + filter_name)

            # Calculate the attenuation
            attenuation = factor * self.attenuation_v / self.av_ebv_ratio

        # Fill in the Ha attenuation manually
        #elif "Halpha" in filter_name or "Ha" in filter_name or "656_1" in filter_name: attenuation = 0.174

        # No precalculated extinction
        else: raise ValueError("Don't have precalculated extinction for this filter")

        # Return the attenuation
        return attenuation

    # -----------------------------------------------------------------

    def irsa_extinction(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        filter_name = str(fltr)

        # Other bands for which attenuation is listed by IRSA
        #elif filter_name in irsa_names:

        irsa_name = irsa_names[filter_name]

        # Find the index of the corresponding table entry
        index = tables.find_index(self.table, irsa_name)

        # Get the attenuation
        attenuation = self.table["A_SandF"][index]

        # All other bands: set attenuation to zero
        #else: #attenuation = 0.0
        # NO: RAISE AND ERROR
        #else: raise ValueError("Cannot determine the extinction for the '"+  filter_name + "' filter")

        # Return the galactic attenuation coefficient
        return attenuation

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

            except ValueError as e:
                if ignore_errors:
                    log.warning("Could not determine the extinction for the " + str(fltr) + " filter:")
                    print(e)
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
            extinction = self.extinction_for_filter(filter_name)

            # Correct the flux value (and error) for galactic extinction
            correction_factor = 10 ** (0.4 * extinction)
            sed["Flux"][i] *= correction_factor
            sed["Error-"][i] *= correction_factor
            sed["Error+"][i] *= correction_factor

# -----------------------------------------------------------------
