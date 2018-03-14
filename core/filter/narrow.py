#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.filter.narrow Contains the NarrowBandFilter class.

# -----------------------------------------------------------------

# Import standard modules
from collections import defaultdict

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from .filter import Filter
from ..basics.map import Map
from ..tools import strings, types

# -----------------------------------------------------------------

# Identifiers
identifiers = dict()
identifiers["Halpha"] = Map(bands=["Ha", "H alpha", "H-alpha", "H-a", "Halpha", "H_alpha", "Mosaic Halpha", "H_alph", "H-alph"], wavelength="656.28 nm")
identifiers["ALMA"] = Map(observatories=["ALMA", "APEX"], frequency_range="84 GHz > 950 GHz")

# -----------------------------------------------------------------

def get_filters_for_regime(regime, subregime=None, categorize=False):

    """
    This function ...
    :param regime:
    :param subregime:
    :param categorize:
    :return:
    """

    from ...magic.tools.wavelengths import wavelength_range_for_regime
    wavelength_range = wavelength_range_for_regime(regime, subregime=subregime)
    return get_filters(wavelength_range.min, wavelength_range.max, categorize=categorize)

# -----------------------------------------------------------------

def get_filters(min_wavelength=None, max_wavelength=None, categorize=False):

    """
    This function ...
    :param min_wavelength:
    :param max_wavelength:
    :param categorize:
    :return:
    """

    raise NotImplementedError("Not yet implemented")

# -----------------------------------------------------------------

def categorize_filters(wavelength_range=None):

    """
    This function ...
    :param wavelength_range:
    :return:
    """

    narrow = defaultdict(list)

    # Categorize
    for spec in identifiers: narrow[spec].append(spec)

    # Filter out those outside of the wavelength range
    if wavelength_range is not None:

        # Loop over the filters
        for spec in narrow.keys():

            if defines_wavelength(spec):
                fltr = NarrowBandFilter(spec)
                wavelength = fltr.wavelength
                wavelength_range_filter = None
            else:
                wavelength = None
                wavelength_range_filter = wavelength_range_for_spec(spec)

            # Check wavelength
            if wavelength is not None and wavelength not in wavelength_range: narrow.pop(spec)
            if wavelength_range_filter is not None and (wavelength_range_filter.min not in wavelength_range or wavelength_range_filter.max not in wavelength_range): narrow.pop(spec)

    # Return
    return narrow

# -----------------------------------------------------------------

def categorized_filters_sorted_labels(narrow):

    """
    Ths function ...
    :param narrow:
    :return:
    """

    return sorted(narrow.keys(), key=lambda x: identifiers.keys().index(narrow[x][0]))

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

            yield band
            yield "the " + band + " band"
            yield "the " + band + "-band"
            yield band + " band"

    if "observatories" in identifier:

        if "bands" in identifier:

            for observatory in strings.case_combinations_list(identifier.observatories):
                for band in strings.case_combinations_list(identifier.bands):
                    for string in strings.generate_from_two_parts(observatory, band): yield string

        else:

            for observatory in strings.case_combinations_list(identifier.observatories):
                yield observatory

    if "instruments" in identifier:

        if "bands" in identifier:

            for instrument in strings.case_combinations_list(identifier.instruments):
                for band in strings.case_combinations_list(identifier.bands):
                    for string in strings.generate_from_two_parts(instrument, band): yield string

        else:

            for instrument in strings.case_combinations_list(identifier.instruments):
                yield instrument

# -----------------------------------------------------------------

def generate_aliases_no_ranges():

    """
    This function ...
    :return:
    """

    for spec in identifiers:
        if "wavelength_range" in identifiers[spec] or "frequency_range" in identifiers[spec]: continue
        for alias in generate_aliases(identifiers[spec]):
            yield spec, alias

# -----------------------------------------------------------------

def generate_aliases_ranges():

    """
    This function ...
    """

    for spec in identifiers:
        if "wavelength_range" not in identifiers[spec] and "frequency_range" not in identifiers[spec]: continue
        for alias in generate_aliases(identifiers[spec]):
            yield spec, alias

# -----------------------------------------------------------------

def defines_wavelength(spec):

    """
    This function ...
    :param spec:
    :return:
    """

    identifier = identifiers[spec]
    return "wavelength" in identifier

# -----------------------------------------------------------------

def wavelength_range_for_spec(spec):

    """
    This function ...
    :param spec:
    :return:
    """

    # Import
    from ..tools import parsing

    identifier = identifiers[spec]

    # Parse the wavelength range
    if "wavelength_range" in identifier: wavelength_range = parsing.quantity_range(identifier.wavelength_range)
    elif "frequency_range" in identifier: wavelength_range = parsing.quantity_range(identifier.frequency_range).to("micron", equivalencies=spectral())
    else: raise ValueError("Wavelength range or frequency range is not defined for filter spec " + spec)

    # Return the wavelength range
    return wavelength_range

# -----------------------------------------------------------------

def wavelengths_for_spec(spec):

    """
    This function ...
    :param spec:
    :return:
    """

    # Import statements
    from ..units.parsing import parse_quantity

    wavelengths = dict()

    identifier = identifiers[spec]
    if defines_wavelength(spec): wavelengths[spec] = parse_quantity(identifier.wavelength)
    else:
        # Get the wavelength range
        wavelength_range = wavelength_range_for_spec(spec)

        # Add the minimum and maximum wavelength
        min_wavelength = wavelength_range.min
        max_wavelength = wavelength_range.max
        wavelengths[spec + " min"] = min_wavelength
        wavelengths[spec + " max"] = max_wavelength

    # Return the wavelengths
    return wavelengths

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

    if "observatories" in identifier:
        description += "on the " + identifier.observatories[0] + " observatory "

    # Return the description
    return description.strip()

# -----------------------------------------------------------------

class NarrowBandFilter(Filter):

    """
    An instance of the NarrowBandFilter class represents a narrow band filter around a certain wavelength.
    """

    # -----------------------------------------------------------------

    cached = {}

    # -----------------------------------------------------------------

    def __new__(cls, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        if len(args) == 0: return super(NarrowBandFilter, cls).__new__(cls)

        filterspec = args[0]
        if filterspec in cls.cached: return cls.cached[filterspec]
        else:
            fltr = super(NarrowBandFilter, cls).__new__(cls)
            fltr.__init__(*args, **kwargs)
            if fltr.true_filter: cls.cached[filterspec] = fltr
            return fltr

    # -----------------------------------------------------------------

    def __init__(self, filterspec, name=None):

        """
        This function ...
        :param filterspec:
        :param name:
        """

        from astropy.units import Quantity
        from ..units.parsing import parse_quantity
        from ..units.stringify import represent_unit

        observatory = None
        instrument = None

        # String is passed
        if types.is_string_type(filterspec):

            # Find exact match
            if filterspec not in identifiers:
                for spec, alias in generate_aliases_no_ranges():
                    if filterspec == alias:
                        filterspec = spec
                        break

                # Break not encountered
                else:

                    for spec, alias in generate_aliases_ranges():
                        if alias in filterspec:
                            identifier = identifiers[spec]
                            filterspec = parse_quantity(filterspec.split(alias)[1], physical_type="length") # get wavelength
                            observatory = identifier.observatories[0] if "observatories" in identifier else None
                            instrument = identifier.instruments[0] if "instruments" in identifier else None
                            name = alias
                            break

                    # Break not encountered
                    else:  raise ValueError("Could not recognize the filter: " + filterspec)

        # String is converted to a valid filterspec
        if types.is_string_type(filterspec):

            identifier = identifiers[filterspec]
            wavelength = identifier.wavelength

            self._observatory = identifier.observatories[0] if "observatories" in identifier else None
            self._instrument = identifier.instruments[0] if "instruments" in identifier else None
            self._band = identifier.bands[0] if "bands" in identifier else None

            self.wavelength = parse_quantity(wavelength)
            self._name = filterspec

            filter_id = self._name
            description = get_filter_description(filterspec)

            true_filter = True

        # String is converted to a wavelength and a name or was a wavelength to begin with
        else:

            self.wavelength = filterspec
            if not isinstance(self.wavelength, Quantity): raise ValueError("Narrow band filter specification should be identifier (string) or wavelength (quantity)")
            self._name = name

            wavelength_as_string = strings.str_from_real_or_integer(self.wavelength.value) + " " + represent_unit(self.wavelength.unit)

            self._observatory = observatory
            self._instrument = instrument

            if self._observatory is None and self._instrument is not None: self._observatory = self._instrument
            if self._instrument is None and self._observatory is not None: self._instrument = self._observatory

            self._band = wavelength_as_string

            filter_id = self._name
            description = "Filter at a wavelength of " + wavelength_as_string

            true_filter = False

        self.true_filter = true_filter
        self.spec = filterspec

        # Call the constructor of the base class
        super(NarrowBandFilter, self).__init__(filter_id, description)

    # -----------------------------------------------------------------

    @property
    def is_halpha(self):

        """
        This property ...
        :return:
        """

        return self._name == "Halpha"

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.instrument is not None: return self.instrument + " " + self.band
        elif self.observatory is not None: return self.observatory + " " + self.band
        #elif self.name is not None: return self.name + " " + self.band
        else: return self.band

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This fucntion produces a string representation of this instance
        :return: 
        """

        prefix = "<NarrowBandFilter "
        suffix = " >"
        return prefix + str(self) + suffix

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self._name

    # -----------------------------------------------------------------

    @property
    def observatory(self):

        """
        This property ...
        :return:
        """

        return self._observatory

    # -----------------------------------------------------------------

    @property
    def instrument(self):

        """
        This property ...
        :return:
        """

        return self._instrument

    # -----------------------------------------------------------------

    @property
    def band(self):

        """
        This property ...
        :return:
        """

        return self._band

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This property ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def effective(self):

        """
        This property ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def min(self):

        """
        THis function ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def max(self):

        """
        This fucntion ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def pivot(self):

        """
        This function ...
        :return:
        """

        return self.wavelength

    # -----------------------------------------------------------------

    @property
    def effective(self):

        """
        This function ...
        :return:
        """

        return self.wavelength

# -----------------------------------------------------------------
