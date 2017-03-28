#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.filter.broad Contains the BroadBandFilter class.
#
# An instance of the BroadBandFilter class in this module represents a particular wavelength band, including its response or
# transmission curve, and allows integrating a given spectrum over the band. A filter instance can be constructed from
# one of the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import sys
import types
import numpy as np
from scipy.interpolate import interp1d
from lxml import etree
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..basics.map import Map
from ..tools import strings
from ..tools.strings import str_from_real_or_integer
from .filter import Filter

# -----------------------------------------------------------------

identifiers = OrderedDict()
identifiers["GALEX.FUV"] = Map(instruments=["GALEX"], bands=["FUV"])
identifiers["GALEX.NUV"] = Map(instruments=["GALEX"], bands=["NUV"])
identifiers["SDSS.u"] = Map(instruments=["SDSS"], bands=["u"])
identifiers["SDSS.g"] = Map(instruments=["SDSS"], bands=["g"])
identifiers["SDSS.r"] = Map(instruments=["SDSS"], bands=["r"])
identifiers["SDSS.i"] = Map(instruments=["SDSS"], bands=["i"])
identifiers["SDSS.z"] = Map(instruments=["SDSS"], bands=["z"])
identifiers["2MASS.H"] = Map(instruments=["2MASS"], bands=["H"])
identifiers["2MASS.J"] = Map(instruments=["2MASS"], bands=["J"])
identifiers["2MASS.K"] = Map(instruments=["2MASS"], bands=["K", "Ks"])
identifiers["UKIDSS.H"] = Map(observatories=["UKIRT"], instruments=["WFCAM"], bands=["H"], surveys=["UKIDSS"])
identifiers["UKIDSS.J"] = Map(observatories=["UKIRT"], instruments=["WFCAM"], bands=["J"], surveys=["UKIDSS"])
identifiers["UKIDSS.K"] = Map(observatories=["UKIRT"], instruments=["WFCAM"], bands=["K"], surveys=["UKIDSS"])
identifiers["UKIDSS.Y"] = Map(observatories=["UKIRT"], instruments=["WFCAM"], bands=["Y"], surveys=["UKIDSS"])
identifiers["UKIDSS.Z"] = Map(observatories=["UKIRT"], instruments=["WFCAM"], bands=["Z"], surveys=["UKIDSS"])
identifiers["IRAC.I1"] = Map(observatories=["Spitzer"], instruments=["IRAC"], bands=["I1"], channel=1, wavelength="3.6 micron")
identifiers["IRAC.I2"] = Map(observatories=["Spitzer"], instruments=["IRAC"], bands=["I2"], channel=2, wavelength="4.5 micron")
identifiers["IRAC.I3"] = Map(observatories=["Spitzer"], instruments=["IRAC"], bands=["I3"], channel=3, wavelength="5.8 micron")
identifiers["IRAC.I4"] = Map(observatories=["Spitzer"], instruments=["IRAC"], bands=["I4"], channel=4, wavelength="8.0 micron")
identifiers["WISE.W1"] = Map(instruments=["WISE"], bands=["W1"], channel=1, wavelength="3.4 micron")
identifiers["WISE.W2"] = Map(instruments=["WISE"], bands=["W2"], channel=2, wavelength="4.6 micron")
identifiers["WISE.W3"] = Map(instruments=["WISE"], bands=["W3"], channel=3, wavelength="12 micron")
identifiers["WISE.W4"] = Map(instruments=["WISE"], bands=["W4"], channel=4, wavelength="22 micron")
identifiers["MIPS.24mu"] = Map(observatories=["Spitzer"], instruments=["MIPS"], wavelength="24 micron")
identifiers["MIPS.70mu"] = Map(observatories=["Spitzer"], instruments=["MIPS"], wavelength="70 micron")
identifiers["MIPS.160mu"] = Map(observatories=["Spitzer"], instruments=["MIPS"], wavelength="160 micron")
identifiers["Pacs.blue"] = Map(observatories=["Herschel"], instruments=["Pacs"], wavelength="70 micron", bands=["blue"])
identifiers["Pacs.green"] = Map(observatories=["Herschel"], instruments=["Pacs"], wavelength="100 micron", bands=["green"])
identifiers["Pacs.red"] = Map(observatories=["Herschel"], instruments=["Pacs"], wavelength="160 micron", bands=["red"])
identifiers["SPIRE.PSW"] = Map(observatories=["Herschel"], instruments=["SPIRE"], bands=["PSW", "PSW_ext"], wavelength="250 micron")
identifiers["SPIRE.PMW"] = Map(observatories=["Herschel"], instruments=["SPIRE"], bands=["PMW", "PMW_ext"], wavelength="350 micron")
identifiers["SPIRE.PLW"] = Map(observatories=["Herschel"], instruments=["SPIRE"], bands=["PLW", "PLW_ext"], wavelength="500 micron")
identifiers["Planck_350"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="350 micron", frequency="857 GHz", bands=["857"])
identifiers["Planck_550"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="550 micron", frequency="545 GHz", bands=["545"])
identifiers["Planck_850"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="850 micron", frequency="353 GHz", bands=["353"])
identifiers["Planck_1380"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="1380 micron", frequency="217 GHz", bands=["217"])
identifiers["Planck_2100"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="2100 micron", frequency="143 GHz", bands=["143"])
identifiers["Planck_3000"] = Map(observatories=["Planck"], instruments=["HFI"], wavelength="3000 micron", frequency="100 GHz", bands=["100"])
identifiers["Planck_4260"] = Map(observatories=["Planck"], instruments=["LFI"], wavelength="4260 micron", frequency="70 GHz", bands=["070"])
identifiers["Planck_6810"] = Map(observatories=["Planck"], instruments=["LFI"], wavelength="6810 micron", frequency="44 GHz", bands=["044"])
identifiers["Planck_10600"] = Map(observatories=["Planck"], instruments=["LFI"], wavelength="10600 micron", frequency="30 GHz", bands=["030"])
identifiers["Johnson.U"] = Map(system="Johnson", bands=["U"])
identifiers["Johnson.B"] = Map(system="Johnson", bands=["B"])
identifiers["Johnson.V"] = Map(system="Johnson", bands=["V"])
identifiers["Johnson.R"] = Map(system="Johnson", bands=["R"])
identifiers["Johnson.I"] = Map(system="Johnson", bands=["I"])
identifiers["IRAS.12mu"] = Map(instruments=["IRAS"], wavelength="12 micron")
identifiers["IRAS.25mu"] = Map(instruments=["IRAS"], wavelength="25 micron")
identifiers["IRAS.60mu"] = Map(instruments=["IRAS"], wavelength="60 micron")
identifiers["IRAS.100mu"] = Map(instruments=["IRAS"], wavelength="100 micron")
identifiers["Swift.UVOT.Uband"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["U"])
identifiers["Swift.UVOT.Bband"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["B"])
identifiers["Swift.UVOT.Vband"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["V"])
identifiers["UVOT.UVW2"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["UVW2", "W2"])
identifiers["UVOT.UVM2"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["UVM2", "M2"])
identifiers["UVOT.UVW1"] = Map(observatories=["Swift"], instruments=["UVOT"], bands=["UVW1", "W1"])
identifiers["ALMA.3"] = Map(observatories=["ALMA", "APEX"], channel=3)
identifiers["ALMA.4"] = Map(observatories=["ALMA", "APEX"], channel=4)
identifiers["ALMA.5"] = Map(observatories=["ALMA", "APEX"], channel=5)
identifiers["ALMA.6"] = Map(observatories=["ALMA", "APEX"], channel=6)
identifiers["ALMA.7"] = Map(observatories=["ALMA", "APEX"], channel=7)
identifiers["ALMA.8"] = Map(observatories=["ALMA", "APEX"], channel=8)
identifiers["ALMA.9"] = Map(observatories=["ALMA", "APEX"], channel=9)
identifiers["ALMA.10"] = Map(observatories=["ALMA", "APEX"], channel=10)
identifiers["OIG.U"] = Map(observatories=["TNG"], instruments=["Oig"], bands=["U"])
identifiers["OIG.V"] = Map(observatories=["TNG"], instruments=["Oig"], bands=["V"])
identifiers["OIG.B"] = Map(observatories=["TNG"], instruments=["Oig"], bands=["B"])
identifiers["OIG.R"] = Map(observatories=["TNG"], instruments=["Oig"], bands=["R"])
identifiers["NICS.H"] = Map(observatories=["TNG"], instruments=["Nics"], bands=["H"])
identifiers["NICS.J"] = Map(observatories=["TNG"], instruments=["Nics"], bands=["J"])
#identifiers["NICS.Js"] = Map(observatories=["TNG"], instruments=["Nics"], bands=["Js"]) we removed this because the spec was ambigious ...
identifiers["NICS.K"] = Map(observatories=["TNG"], instruments=["Nics"], bands=["K"])

# -----------------------------------------------------------------

# Filtername , Frequency, filename
# Planck_350    857 GHz   857   HFI
# Planck_550    545 GHz   545   HFI
# Planck_850    352 GHz   353   HFI
# Planck_1380   222 GHz   217   HFI
# Planck_2100   142 GHz   143   HFI
# Planck_3000   100 GHz   100   HFI
# Planck_4260   70 GHz    070   LFI
# Planck_6810   44 GHz    044   LFI
# Planck_10600  28 GHz    030   LFI

planck_info = dict()
planck_info["Planck_350"] = ("857", "HFI")
planck_info["Planck_550"] = ("545", "HFI")
planck_info["Planck_850"] = ("353", "HFI")
planck_info["Planck_1380"] = ("217", "HFI")
planck_info["Planck_2100"] = ("143", "HFI")
planck_info["Planck_3000"] = ("100", "HFI")
planck_info["Planck_4260"] = ("070", "LFI")
planck_info["Planck_6810"] = ("044", "LFI")
planck_info["Planck_10600"] = ("030", "LFI")

# -----------------------------------------------------------------

# ALMA BAND   Frequency range (GHz)      Wavelength range (mm)   Wavelength range (micron)
# 1             31 - 45
# 2             67 - 90
# 3             84 - 116                    3.57 - 2.59             2590 - 3570
# 4             125 - 163                   2.40 - 1.84             1840 - 2400
# 5             162 - 211                   1.84 - 1.42             1420 - 1840
# 6             211 - 275                   1.42 - 1.09             1090 - 1420
# 7             275 - 373                   1.09 - 0.80             800 - 1090
# 8             385 - 500                   0.78 - 0.60             600 - 780
# 9             602 - 720                   0.50 - 0.42             420 - 500
# 10            787 - 950                   0.38 - 0.32             320 - 380

alma_ranges = dict()
alma_ranges[3] = (2590, 3570)
alma_ranges[4] = (1840, 2400)
alma_ranges[5] = (1420, 1840)
alma_ranges[6] = (1090, 1420)
alma_ranges[7] = (800, 1090)
alma_ranges[8] = (600, 780)
alma_ranges[9] = (420, 500)
alma_ranges[10] = (320, 380)

# -----------------------------------------------------------------

def get_filter_description(spec):

    """
    This function ...
    :param spec:
    :return:
    """

    identifier = identifiers[spec]
    description = "the "

    if "bands" in identifier: description += identifier.bands[0] + " "
    elif "channel" in identifier: description += strings.num_to_ith(identifier.channel) + " "
    description += "band "

    if "wavelength" in identifier and "frequency" in identifier:
        description += "at " + identifier.wavelength + " or " + identifier.frequency + " "
    elif "wavelength" in identifier:
        description += "at " + identifier.wavelength + " "
    elif "frequency" in identifier:
        description += "at " + identifier.frequency + " "

    if "instruments" in identifier:
        description += "of the " + identifier.instruments[0] + " instrument "

    if "system" in identifier:
        description += "of the " + identifier.system + " system "

    if "observatories" in identifier:
        description += "on the " + identifier.observatories[0] + " observatory "

    # Return the description
    return description.strip()

# -----------------------------------------------------------------

def get_filter_descriptions():

    """
    This function ...
    :return:
    """

    descriptions = dict()

    # Loop over the filters
    for spec in identifiers: descriptions[spec] = get_filter_description(spec)

    # Return the dictionary
    return descriptions

# -----------------------------------------------------------------

def is_sdss_or_johnson(identifier):

    """
    This function ...
    :param identifier:
    :return:
    """

    if "instruments" in identifier:

        if "SDSS" in identifier.instruments: return True

    if "system" in identifier:

        if identifier.system == "Johnson": return True
        else: return False

    return False

# -----------------------------------------------------------------

def is_sdss_2mass_or_johnson(identifier):

    """
    This function ...
    :param identifier:
    :return:
    """

    if "instruments" in identifier:

        if "SDSS" in identifier.instruments: return True
        if "2MASS" in identifier.instruments: return True

    if "system" in identifier:

        if identifier.system == "Johnson": return True
        else: return False

    return False

# -----------------------------------------------------------------

def generate_aliases(identifier):

    """
    This function ...
    :param identifier:
    :return:
    """

    # Just the band
    if "bands" in identifier:
        for band in strings.case_combinations_list(identifier.bands, also_one_letter=False):

            #if len(band) == 1 and not is_sdss_2mass_or_johnson(identifier): continue
            if len(band) == 1 and not is_sdss_2mass_or_johnson(identifier): continue

            yield band
            yield "the " + band + " band"
            yield "the " + band + "-band"
            yield band + " band"

    # Combinations of system and band
    if "system" in identifier and "bands" in identifier:

        for band in strings.case_combinations_list(identifier.bands, also_one_letter=False):
            for string in strings.generate_from_two_parts(identifier.system, band, also_reverse=True): yield string

    # Combinations of instrument with band
    if "instruments" in identifier and "bands" in identifier:

        for instrument in strings.case_combinations_list(identifier.instruments):
            for band in strings.case_combinations_list(identifier.bands, also_one_letter=False):

                for string in strings.generate_from_two_parts(instrument, band, also_reverse=True): yield string
                for string in strings.generate_from_two_parts("the " + instrument, band + "-band"): yield string
                for string in strings.generate_from_two_parts("the " + instrument, band + " band"): yield string

    # Combinations of instrument with channel
    if "channel" in identifier and "instruments" in identifier:
        for instrument in strings.case_combinations_list(identifier.instruments):
            for string in strings.generate_from_two_parts(instrument, str_from_real_or_integer(identifier.channel)): yield string

    # Combinations of observatory with channel
    if "channel" in identifier and "observatories" in identifier:
        for instrument in strings.case_combinations_list(identifier.observatories):
            for string in strings.generate_from_two_parts(instrument, str_from_real_or_integer(identifier.channel)): yield string

    # Combinations of observatory with band
    if "observatories" in identifier and "bands" in identifier:

        for observatory in strings.case_combinations_list(identifier.observatories):
            for band in strings.case_combinations_list(identifier.bands, also_one_letter=False):

                for string in strings.generate_from_two_parts(observatory, band, also_reverse=True): yield string
                for string in strings.generate_from_two_parts("the " + observatory, band + "-band"): yield string
                for string in strings.generate_from_two_parts("the " + observatory, band + "-band"): yield string

    # Combinations of instrument and wavelength
    if "wavelength" in identifier and "instruments" in identifier:

        from ..tools import parsing
        wavelength = parsing.quantity(identifier.wavelength)

        for instrument in strings.case_combinations_list(identifier.instruments):
            for wavelength_string in strings.quantity_combinations(wavelength):
                for string in strings.generate_from_two_parts(instrument, wavelength_string): yield string
                for string in strings.generate_from_two_parts("the " + instrument, wavelength_string + " band"): yield string

    # Combinations of observatory with wavelength
    if "observatories" in identifier and "wavelength" in identifier:

        from ..tools import parsing
        wavelength = parsing.quantity(identifier.wavelength)

        for observatory in identifier.observatories:
            for wavelength_string in strings.quantity_combinations(wavelength):
                for string in strings.generate_from_two_parts(observatory, wavelength_string): yield string
                for string in strings.generate_from_two_parts("the " + observatory, wavelength_string + " band"): yield string

    # Combinations of instrument and frequency
    if "frequency" in identifier and "instruments" in identifier:

        from ..tools import parsing
        frequency = parsing.quantity(identifier.frequency)

        for instrument in identifier.instruments:
            for frequency_string in strings.quantity_combinations(frequency):
                for string in strings.generate_from_two_parts(instrument, frequency_string): yield string
                for string in strings.generate_from_two_parts("the " + instrument, frequency_string + " band"): yield string

    # Combinations of observatory and frequency
    if "observatories" in identifier and "frequency" in identifier:

        from ..tools import parsing
        frequency = parsing.quantity(identifier.frequency)

        for observatory in identifier.observatories:
            for frequency_string in strings.quantity_combinations(frequency):
                for string in strings.generate_from_two_parts(observatory, frequency_string): yield string
                for string in strings.generate_from_two_parts("the " + observatory, frequency_string + " band"): yield string

    # Combinations of survey and band
    if "surveys" in identifier and "bands" in identifier:

        for survey in strings.case_combinations_list(identifier.surveys):
            for band in strings.case_combinations_list(identifier.bands, also_one_letter=False):

                for string in strings.generate_from_two_parts(survey, band, also_reverse=True): yield string
                for string in strings.generate_from_two_parts("the " + survey, band + "-band"): yield string
                for string in strings.generate_from_two_parts("the " + survey, band + " band"): yield string

    # Combinations of survey and channel
    if "channel" in identifier and "surveys" in identifier:
        for survey in strings.case_combinations_list(identifier.surveys):
            for string in strings.generate_from_two_parts(survey, str_from_real_or_integer(identifier.channel)): yield string

    # Combinations of survey and wavelength
    if "wavelength" in identifier and "surveys" in identifier:

        from ..tools import parsing
        wavelength = parsing.quantity(identifier.wavelength)

        for survey in strings.case_combinations_list(identifier.surveys):
            for wavelength_string in strings.quantity_combinations(wavelength):
                for string in strings.generate_from_two_parts(survey, wavelength_string): yield string
                for string in strings.generate_from_two_parts("the " + survey, wavelength_string + " band"): yield string

    # Combinations of survey and frequency
    if "frequency" in identifier and "surveys" in identifier:

        from ..tools import parsing
        frequency = parsing.quantity(identifier.frequency)

        for survey in identifier.surveys:
            for frequency_string in strings.quantity_combinations(frequency):
                for string in strings.generate_from_two_parts(survey, frequency_string): yield string
                for string in strings.generate_from_two_parts("the " + survey, frequency_string + " band"): yield string

# -----------------------------------------------------------------

def generate_all_aliases():

    """
    This function ...
    :return:
    """

    for spec in identifiers:
        for alias in generate_aliases(identifiers[spec]):
            yield spec, alias

# -----------------------------------------------------------------

## An instance of the BroadBandFilter class represents a particular wavelength bandpass, including its response or
# transmission curve and some basic properties such as its mean and pivot wavelengths. The class provides a function to
# integrate a given spectrum over the band. A filter instance can be constructed by name from one of
# the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.
#
# The precise formalae involved in the integration over the filter and the calculation of the pivot wavelength
# depend on whether the instrument counts photons (photon counter) or measures energy (bolometer).
#
class BroadBandFilter(Filter):

    # ---------- Constructing -------------------------------------

    ## The constructor constructs a new Filter instance in one of two ways, depending on the type of the argument.
    #
    # If \em filterspec is a tuple with two numbers, a bolometer-type filter is constructed with a uniform
    # transmission curve in the indicated range. The (min,max) wavelength values must be expressed in micron.
    #
    # If \em filterspec is a string, the constructor locates a "$PYTHONPATH/dat/filter" directory, and searches
    # in that directory for a VOTable resource file with a name that contains the specified string. If there is
    # exactly one such file, its contents is loaded into the new Filter instance. If there is no such file, or
    # if the specified string matches multiple file names, the constructor raises an error.
    #
    # The VOTable resource files do not properly specify the type of filter (photon counter or bolometer).
    # The constructor classifies filters based on the instrument name in the filter identifier.
    # While this simple heuristic works for now, it is not guaranteed to be correct for future cases.
    #
    # The tables below list the filters available at the time of writing. More filter definitions can be
    # downloaded from the filter profile service web site http://svo2.cab.inta-csic.es/theory/fps.
    #
    # For ALMA, the filter profiles have been downloaded from http://www.obs.u-bordeaux1.fr/radio/SBontemps/atmo/alma.html.
    #
    # A list of all predefined filters can be generated with the PTS command 'pts filters'.
    #
    def __init__(self, filterspec, name=None):

        # Check aliases if the filterspec is not exactly equal to predefined specs
        if isinstance(filterspec, types.StringTypes):
            if filterspec not in identifiers:
                for spec, alias in generate_all_aliases():
                    if filterspec == alias:
                        filterspec = spec
                        break
                # Break not encountered
                else: raise ValueError("Could not recognize the filter: " + filterspec)

        # Planck filters have to be handled seperately
        if isinstance(filterspec, types.StringTypes) and "planck" in filterspec.lower():

            # Load properties
            wavelengths, transmissions, instrument, frequency_string = load_planck(filterspec)

            # Determine ID and description
            filter_id = "Planck/" + instrument + "." + frequency_string
            description = instrument + " " + frequency_string

            true_filter = True

            # Initialize
            self._initialize2(wavelengths, transmissions)

        # ALMA filters have to be handled seperately
        elif isinstance(filterspec, types.StringTypes) and "alma" in filterspec.lower():

            # Load
            wavelengths, transmissions = load_alma(filterspec)

            # Determine ID and description
            identifier = identifiers[filterspec]
            filter_id = "ALMA/ALMA." + str(identifier.channel)
            description = "ALMA " + str(identifier.channel)

            true_filter = True

            # Initialize
            self._initialize2(wavelengths, transmissions)

        # string --> load from appropriate resource file
        elif isinstance(filterspec, types.StringTypes):

            # Load
            min_wavelength, max_wavelength, center_wavelength, mean_wavelength, eff_wavelength, filter_id, description, \
            fwhm, eff_width, photon_counter, wavelengths, transmissions = load_svo(filterspec)

            # Set properties
            self._WavelengthMin = min_wavelength
            self._WavelengthMax = max_wavelength
            self._WavelengthCen = center_wavelength
            self._WavelengthMean = mean_wavelength
            self._WavelengthEff = eff_wavelength
            self._EffWidth = eff_width
            self._FWHM = fwhm
            self._Wavelengths = wavelengths
            self._Transmission = transmissions
            self._PhotonCounter = photon_counter

            true_filter = True

            # Initialize
            self._initialize1()

        # range --> construct ad hoc uniform bolometer
        else:

            if name is None: name = "Uniform"
            self._WavelengthMin, self._WavelengthMax = map(float, filterspec)
            self._WavelengthCen = 0.5 * (self._WavelengthMin + self._WavelengthMax)
            self._WavelengthMean = self._WavelengthCen
            self._WavelengthEff = self._WavelengthCen
            self._Wavelengths = np.array((self._WavelengthMin,self._WavelengthMax))
            self._Transmission = np.ones((2,))
            self._PhotonCounter = False
            self._IntegratedTransmission = self._WavelengthMax - self._WavelengthMin
            self._WavelengthPivot = np.sqrt(self._WavelengthMin * self._WavelengthMax)
            self._EffWidth = self._WavelengthMax - self._WavelengthMin
            self._FWHM = None

            true_filter = False

            # Determine filter ID and description
            filter_id = name + "_[{},{}]".format(self._WavelengthMin, self._WavelengthMax)
            description = name + " filter in range [{},{}]".format(self._WavelengthMin, self._WavelengthMax)

        self.true_filter = true_filter
        self.spec = filterspec

        # Call the constructor of the base class
        super(BroadBandFilter, self).__init__(filter_id, description)

    # -----------------------------------------------------------------

    @classmethod
    def from_instrument_and_band(cls, instrument, band):

        """
        This function ...
        :param instrument:
        :param band:
        :return:
        """

        return cls(instrument + "." + band)

    # -----------------------------------------------------------------

    @classmethod
    def sdss(cls, band):

        """
        This function ...
        :param band:
        :return:
        """

        return cls("SDSS " + band)

    # -----------------------------------------------------------------

    def _initialize1(self):

        """
        This function ...
        :return:
        """

        # calculate the pivot wavelength
        if self._PhotonCounter:
            integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission * self._Wavelengths)
            integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission / self._Wavelengths)
        else:
            integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission)
            integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission / (self._Wavelengths ** 2))
        self._IntegratedTransmission = integral1
        self._WavelengthPivot = np.sqrt(integral1 / integral2)

    # -----------------------------------------------------------------

    def _initialize2(self, wavelengths, transmissions):

        """
        This function ...
        :param wavelengths:
        :param transmissions:
        :return:
        """

        self._WavelengthMin = np.min(wavelengths)
        self._WavelengthMax = np.max(wavelengths)
        self._WavelengthCen = 0.5 * (self._WavelengthMin + self._WavelengthMax)

        self._Wavelengths = wavelengths
        self._Transmission = transmissions

        # Calculate integrals
        integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission)
        integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission / (self._Wavelengths ** 2))

        integral3 = np.trapz(x=self._Wavelengths, y=self._Transmission * self._Wavelengths)

        # Mean wavelength, ID and description
        self._WavelengthMean = integral3 / integral1
        self._WavelengthEff = None

        # Not photon counter
        self._PhotonCounter = False

        # Calculate integrated transmission and pivot wavelength
        self._IntegratedTransmission = integral1
        self._WavelengthPivot = np.sqrt(integral1 / integral2)

        self._EffWidth = None
        self._FWHM = None

    # ---------- Special functions --------------------------------------

    def __eq__(self, other):

        """
        This function compares two filter instances
        :param other:
        :return:
        """

        return str(other) == str(self)

    # ---------- Retrieving information -------------------------------

    @property
    def aliases(self):
        # returns the aliases
        return list(generate_aliases(identifiers[self.name]))

    @property
    def skirt_description(self): # returns the name as defined in the SKIRT LuminosityStellarCompNormalization class

        if self.name == "GALEX.FUV": return "FUV"
        elif self.name == "GALEX.NUV": return "NUV"
        elif self.name == "Johnson.U": return "U"
        elif self.name == "Johnson.B": return "B"
        elif self.name == "Johnson.V": return "V"
        elif self.name == "Johnson.R": return "R"
        elif self.name == "Johnson.I": return "I"
        elif self.name == "2MASS.J": return "J"
        elif self.name == "2MASS.H": return "H"
        elif self.name == "2MASS.Ks": return "K"
        elif self.name == "SDSS.u": return "SDSSu"
        elif self.name == "SDSS.g": return "SDSSg"
        elif self.name == "SDSS.r": return "SDSSr"
        elif self.name == "SDSS.i": return "SDSSi"
        elif self.name == "SDSS.z": return "SDSSz"
        elif self.name == "IRAC.I1": return "IRAC1"
        elif self.name == "IRAC.I2": return "IRAC2"
        elif self.name == "WISE.W1": return "WISE1"
        elif self.name == "WISE.W2": return "WISE2"
        else: raise ValueError("The band " +  self.name + " is not defined in SKIRT")

    @property
    def aniano_name(self): # Returns the name as appearing in the Aniano kernel and psf FITS files

        if self.name == "GALEX.FUV": return "GALEX_FUV"
        elif self.name == "GALEX.NUV": return "GALEX_NUV"
        elif self.name == "IRAC.I1": return "IRAC_3.6"
        elif self.name == "IRAC.I2": return "IRAC_4.5"
        elif self.name == "IRAC.I3": return "IRAC_5.8"
        elif self.name == "IRAC.I4": return "IRAC_8.0"
        elif self.name == "WISE.W1": return "WISE_FRAME_3.4"
        elif self.name == "WISE.W2": return "WISE_FRAME_4.6"
        elif self.name == "WISE.W3": return "WISE_FRAME_11.6"
        elif self.name == "WISE.W4": return "WISE_FRAME_22.1"
        elif self.name == "MIPS.24mu": return "MIPS_24"
        elif self.name == "MIPS.70mu": return "MIPS_70"
        elif self.name == "MIPS.160mu": return "MIPS_160"
        elif self.name == "Pacs.blue": return "PACS_70"
        elif self.name == "Pacs.green": return "PACS_100"
        elif self.name == "Pacs.red": return "PACS_160"
        elif self.name == "SPIRE.PSW_ext": return "SPIRE_250"
        elif self.name == "SPIRE.PMW_ext": return "SPIRE_350"
        elif self.name == "SPIRE.PLW_ext": return "SPIRE_500"
        else: raise ValueError("The band " + self.name + " is not defined for the Aniano set of kernels")

    ## This function returns the filter in string format ('instrument' 'band')
    def __str__(self):
        if self.true_filter: return self.instrument + " " + self.band
        else: return self.name

    @property
    def name(self):

        if self.true_filter: return self._FilterID.split("/")[1]
        else: return self._FilterID.split("_")[0]

    @property
    def observatory(self):

        if self.true_filter: return self._FilterID.split("/")[0]
        else: return None

    @property
    def instrument(self):

        if self.true_filter: return self._FilterID.split("/")[1].split(".")[0]
        else: return None

    @property
    def band(self):

        if self.true_filter: return self._FilterID.split("/")[1].split(".")[1].replace("_ext", "")
        elif self._FilterID.startswith("Uniform"): return None
        else: return self._FilterID.split("_")[0]

    ## This function returns the mean wavelength for the filter, in micron.
    def meanwavelength(self):
        return self._WavelengthMean

    ## This function returns the mean wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def mean(self):
        from ..basics.unit import parse_unit as u
        return self.meanwavelength() * u("micron")

    ## This function returns the effective wavelength for the filter, in micron.
    def effectivewavelength(self):
        return self._WavelengthEff

    ## This function returns the effective wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def effective(self):
        from ..basics.unit import parse_unit as u
        return self.effectivewavelength() * u("micron") if self._WavelengthEff is not None else None

    ## This function returns the minimum wavelength for the filter, in micron.
    def minwavelength(self):
        return self._WavelengthMin

    ## This function returns the minimum wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def min(self):
        from ..basics.unit import parse_unit as u
        return self.minwavelength() * u("micron")

    ## This function returns the maximum wavelength for the filter, in micron.
    def maxwavelength(self):
        return self._WavelengthMax

    ## This function returns the maximum wavelength for the filter as a quantity. Astropy is imported inside this
    #  function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def max(self):
        from ..basics.unit import parse_unit as u
        return self.maxwavelength() * u("micron")

    ## This function returns the center wavelength for the filter, in micron. The center wavelength is
    # defined as the wavelength halfway between the two points for which filter response or transmission
    # (depending on the filter type) is half maximum.
    def centerwavelength(self):
        return self._WavelengthCen

    ## This function returns the center wavelength as a quantity. Astropy is imported inside this function to
    #  support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def center(self):
        from ..basics.unit import parse_unit as u
        return self.centerwavelength() * u("micron")

    ## This function returns the pivot wavelength for the filter, in micron. The pivot wavelength is defined
    # as the wavelength that connects the filter-averaged wavelength and frequency-style fluxes through
    # \f$\left<F_\nu\right> = \left<F_\lambda\right>\lambda_\text{pivot}^2/c\f$. The value depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int\lambda R(\lambda) \,\mathrm{d}\lambda }
    #     {  \int R(\lambda) \,\mathrm{d}\lambda/\lambda } }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda/\lambda^2 } }. \f]
    def pivotwavelength(self):
        return self._WavelengthPivot

    ## This function returns the pivot wavelength as a quantity. Astropy is imported inside this function to
    #  support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def pivot(self):
        from ..basics.unit import parse_unit as u
        return self.pivotwavelength() * u("micron")

    ## This function returns the effective bandwith, in micron.
    def effective_bandwidth(self):
        return self._EffWidth

    @property
    def bandwidth(self):
        from ..basics.unit import parse_unit as u
        return self.effective_bandwidth() * u("micron") if self._EffWidth is not None else None

    @property
    def fwhm(self):
        from ..basics.unit import parse_unit as u
        return self._FWHM * u("micron") if self._FWHM is not None else None

    @property
    def fwhm_range(self):
        if self.fwhm is None: return None
        from ..basics.range import QuantityRange
        return QuantityRange(self.mean - self.fwhm, self.mean + self.fwhm)

    # ---------- Integrating --------------------------------------

    ## This function calculates and returns the filter-averaged value \f$\left<F_\lambda\right>\f$ for a given
    # spectral energy distribution \f$F_\lambda(\lambda)\f$. The calculation depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int\lambda F_\lambda(\lambda)R(\lambda) \,\mathrm{d}\lambda }
    #     { \int\lambda R(\lambda) \,\mathrm{d}\lambda }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int F_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda }. \f]
    #
    # The quantities \f$F_\lambda(\lambda)\f$ must be expressed per unit of wavelength (and \em not, for example,
    # per unit of frequency). The resulting \f$\left<F_\lambda\right>\f$ has the same units as the input distribition.
    # \f$F_\lambda(\lambda)\f$ can be expressed in any units (as long as it is per unit of wavelength) and it can
    # represent various quantities; for example a flux density, a surface density, or a luminosity density.
    #
    # The function accepts two arguments:
    # - \em wavelengths: a numpy array specifying the wavelengths \f$\lambda_\ell\f$, in micron, in increasing order,
    #   on which the spectral energy distribution is sampled. The integration is performed on a wavelength grid that
    #   combines the grid points given here with the grid points on which the filter response or transmission curve
    #   is defined.
    # - \em densities: a numpy array specifying the spectral energy distribution(s) \f$F_\lambda(\lambda_\ell)\f$
    #   per unit of wavelength. This can be an array with the same length as \em wavelengths, or a multi-dimensional
    #   array where the last dimension has the same length as \em wavelengths.
    #   The returned result will have the shape of \em densities minus the last (or only) dimension.
    def convolve(self, wavelengths, densities):

        # define short names for the involved wavelength grids
        wa = wavelengths
        wb = self._Wavelengths

        # create a combined wavelength grid, restricted to the overlapping interval
        w1 = wa[ (wa>=wb[0]) & (wa<=wb[-1]) ]
        w2 = wb[ (wb>=wa[0]) & (wb<=wa[-1]) ]
        w = np.unique(np.hstack((w1,w2)))
        if len(w) < 2: return 0

        # log-log interpolate SED and transmission on the combined wavelength grid
        # (use scipy interpolation function for SED because np.interp does not support broadcasting)
        F = np.exp(interp1d(np.log(wa), _log(densities), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._Transmission), left=0., right=0.))

        # perform the integration
        if self._PhotonCounter: return np.trapz(x=w, y=w*F*T) / self._IntegratedTransmission
        else: return np.trapz(x=w, y=F*T) / self._IntegratedTransmission

    ## This function calculates and returns the integrated value for a given spectral energy distribution over the
    #  filter's wavelength range,
    def integrate(self, wavelengths, densities):

        # define short names for the involved wavelength grids
        wa = wavelengths
        wb = self._Wavelengths

        # create a combined wavelength grid, restricted to the overlapping interval
        w1 = wa[(wa >= wb[0]) & (wa <= wb[-1])]
        w2 = wb[(wb >= wa[0]) & (wb <= wa[-1])]
        w = np.unique(np.hstack((w1, w2)))
        if len(w) < 2: return 0

        # log-log interpolate SED and transmission on the combined wavelength grid
        # (use scipy interpolation function for SED because np.interp does not support broadcasting)
        F = np.exp(interp1d(np.log(wa), _log(densities), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._Transmission), left=0., right=0.))

        # perform the integration
        if self._PhotonCounter: return np.trapz(x=w, y=w * F * T)
        else: return np.trapz(x=w, y=F * T)

## This private helper function returns the natural logarithm for positive values, and a large negative number
# (but not infinity) for zero or negative values.
def _log(X):
    zeromask = X<=0
    logX = np.empty(X.shape)
    logX[zeromask] = -750.  # the smallest (in magnitude) negative value x for which np.exp(x) returns zero
    logX[~zeromask] = np.log(X[~zeromask])
    return logX

# -----------------------------------------------------------------

def load_planck(filterspec):

    """
    This function ...
    :param filterspec:
    :return:
    """

    # Import Astropy stuff
    from astropy.units import spectral
    from ..basics.unit import parse_unit as u

    this_path = os.path.dirname(os.path.abspath(__file__))
    core_path = os.path.dirname(this_path)
    planck_transmissions_path = os.path.join(core_path, "dat", "filters", "Planck")

    # Get info
    frequency_string = planck_info[filterspec][0]
    instrument = planck_info[filterspec][1]

    filename = instrument + "_BANDPASS_F" + frequency_string + ".txt"
    filepath = os.path.join(planck_transmissions_path, filename)

    wavenumbers, transmissions = np.loadtxt(filepath, unpack=True, skiprows=2, usecols=(0, 1))

    # Remove zero from the wavenumbers
    if wavenumbers[0] == 0:
        wavenumbers = wavenumbers[1:]
        transmissions = transmissions[1:]

    # HFI
    if instrument == "HFI":
        wavenumbers = wavenumbers * u("1/cm")
        wavelengths = (1.0 / wavenumbers).to("micron").value

    # LFI
    else:
        frequencies = wavenumbers * u("GHz")
        wavelengths = frequencies.to("micron", equivalencies=spectral()).value

    # REVERSE
    wavelengths = np.flipud(wavelengths)
    transmissions = np.flipud(transmissions)

    # Success
    return wavelengths, transmissions, instrument, frequency_string

# -----------------------------------------------------------------

def load_alma(filterspec, pwv=0.2):

    """
    This function ...
    :param filterspec:
    :param pwv: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 or 1 mm
    :return:
    """

    identifier = identifiers[filterspec]

    possible_pwvs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1]
    assert pwv in possible_pwvs, "Pwv must be one of " + str(possible_pwvs)

    # Import Astropy stuff
    from astropy.units import spectral
    from ..basics.unit import parse_unit as u

    this_path = os.path.dirname(os.path.abspath(__file__))
    core_path = os.path.dirname(this_path)
    alma_transmissions_path = os.path.join(core_path, "dat", "filters", "ALMA")

    # data in two-column format. First column: Frequency in GHz, second column: Transmission (to get opacity take -ln(transmission)). The frequency resolution (df) is indicated for each model.
    #filename = "alma-" + str(min_frequency) + "-" + str(max_frequency) + "-" + pwv_string + ".dat" # example: alma-535-750-02.dat"

    filename = "alma-0-2000-" + str(pwv).replace(".", "") + ".dat"
    filepath = os.path.join(alma_transmissions_path, filename)

    # Load the data
    frequencies, transmissions = np.loadtxt(filepath, unpack=True, usecols=(0, 1))

    # Get array of wavelengths
    frequencies = frequencies * u("GHz")
    wavelengths = frequencies.to("micron", equivalencies=spectral()).value

    # REVERSE
    wavelengths = np.flipud(wavelengths)
    transmissions = np.flipud(transmissions)

    # Cut out the relevant part
    min_wavelength, max_wavelength = alma_ranges[identifier.channel]
    mask = (min_wavelength < wavelengths) * (wavelengths < max_wavelength)
    wavelengths = wavelengths[mask]
    transmissions = transmissions[mask]

    # Make sure that the transmission curves are not floating, but are attached to the x axis
    transmissions[0] = 0.0
    transmissions[-1] = 0.0

    # Return
    return wavelengths, transmissions

# -----------------------------------------------------------------

def load_svo(filterspec):

    """
    This function ...
    :param filterspec:
    :return:
    """

    for pythondir in sys.path:
        filterdir = os.path.join(pythondir, "pts", "core", "dat", "filters", "SVO")
        if os.path.isdir(filterdir):

            filterfiles = filter(lambda fn: fn.endswith(".xml") and filterspec in fn, os.listdir(filterdir))
            if len(filterfiles) > 1: raise ValueError("filter spec " + filterspec + " is ambiguous")
            if len(filterfiles) < 1: raise ValueError("no filter found with spec " + filterspec)

            # load the XML tree
            tree = etree.parse(open(os.path.join(filterdir, filterfiles[0]), 'r'))

            # verify the wavelength unit to be Angstrom
            unit = tree.xpath("//RESOURCE/PARAM[@name='WavelengthUnit'][1]/@value")[0]
            if unit != 'Angstrom': raise ValueError("VOTable uses unsupported unit: " + unit)

            # load some basic properties (converting from Angstrom to micron)
            min_wavelength = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMin'][1]/@value")[0])
            max_wavelength = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMax'][1]/@value")[0])
            center_wavelength = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthCen'][1]/@value")[0])
            mean_wavelength = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMean'][1]/@value")[0])
            eff_wavelength = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthEff'][1]/@value")[0])
            filterid = tree.xpath("//RESOURCE/PARAM[@name='filterID'][1]/@value")[0]
            description = tree.xpath("//RESOURCE/PARAM[@name='Description'][1]/@value")[0]
            description = description.replace("&#956;m", "micron")
            fwhm = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='FWHM'][1]/@value")[0])
            eff_width = 1e-4 * float(tree.xpath("//RESOURCE/PARAM[@name='WidthEff'][1]/@value")[0])

            # load the transmission table (converting wavelengths from Angstrom to micron)
            values = np.array(tree.xpath("//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD/text()"), dtype=float)
            if len(values) < 4: raise ValueError("transmission table not found in filter definition")
            wavelengths, transmissions = np.reshape(values, (-1, 2)).T
            wavelengths *= 1e-4

            # determine the filter type (there seems to be no better heuristic than using the instrument name)
            photon_counter = not any(["/" + x in filterid.lower() for x in ("mips", "pacs", "spire")])

            # Return the properties
            return min_wavelength, max_wavelength, center_wavelength, mean_wavelength, eff_wavelength, filterid, \
                   description, fwhm, eff_width, photon_counter, wavelengths, transmissions

    # Not recognized
    raise ValueError("Could not detect the filter")

# -----------------------------------------------------------------
