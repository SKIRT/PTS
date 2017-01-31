#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.filter.narrow Contains the NarrowBandFilter class.

# -----------------------------------------------------------------

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from .filter import Filter
from ..basics.map import Map
from ..tools import strings

# -----------------------------------------------------------------

# Identifiers
identifiers = dict()
identifiers["Halpha"] = Map(bands=["Ha", "H alpha", "H-alpha", "H-a", "Halpha", "H_alpha", "Mosaic Halpha", "H_alph", "H-alph"], wavelength="656.28 nm")
identifiers["ALMA"] = Map(observatories=["ALMA", "APEX"], frequency_range="84 GHz > 950 GHz")

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
    from ..basics.quantity import parse_quantity

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

    def __init__(self, filterspec, name=None):

        """
        This function ...
        :param filterspec:
        :param name:
        """

        from astropy.units import Quantity
        from ..basics.quantity import parse_quantity
        from ..basics.unit import represent_unit

        observatory = None
        instrument = None

        # String is passed
        if isinstance(filterspec, basestring):

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
                            filterspec = parse_quantity(filterspec.split(alias)[1]) # get wavelength
                            observatory = identifier.observatories[0] if "observatories" in identifier else None
                            instrument = identifier.instruments[0] if "instruments" in identifier else None
                            name = alias
                            break

                    # Break not encountered
                    else:  raise ValueError("Could not recognize the filter: " + filterspec)

        # String is converted to a valid filterspec
        if isinstance(filterspec, basestring):

            identifier = identifiers[filterspec]
            wavelength = identifier.wavelength

            self._observatory = identifier.observatories[0] if "observatories" in identifier else None
            self._instrument = identifier.instruments[0] if "instruments" in identifier else None
            self._band = identifier.bands[0] if "bands" in identifier else None

            self.wavelength = parse_quantity(wavelength)
            self._name = filterspec

            filter_id = self._name
            description = get_filter_description(filterspec)

        # String is converted to a wavelength and a name or was a wavelength to begin with
        else:

            self.wavelength = filterspec
            if not isinstance(self.wavelength, Quantity): raise ValueError("Narrow band filter specification should be identifier (string) or wavelength (quantity)")
            self._name = name

            wavelength_as_string = strings.str_from_real_or_integer(self.wavelength.value) + " " + represent_unit(self.wavelength.unit)

            self._observatory = observatory
            self._instrument = instrument
            self._band = wavelength_as_string

            filter_id = self._name
            description = "Filter at a wavelength of " + wavelength_as_string

        # Call the constructor of the base class
        super(NarrowBandFilter, self).__init__(filter_id, description)

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
