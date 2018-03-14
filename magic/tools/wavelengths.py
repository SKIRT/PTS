#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.wavelengths Provides ...

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.map import Map
from ...core.basics.range import QuantityRange
from ...core.units.parsing import parse_unit as u
from ...core.basics.range import QuantityRange

# -----------------------------------------------------------------

black_body_wavelength_range = QuantityRange(50., 2000., "micron")

# -----------------------------------------------------------------

extinction_wavelength_range = QuantityRange(0.005, 10., "micron")

# -----------------------------------------------------------------

# Define the wavelength ranges
spectrum_wavelengths = OrderedDict([(("UV", "EUV"), (0.01, 0.121)),
                                   (("UV", "Lyman-alpha"), (0.121, 0.122)),
                                   (("UV", "FUV"), (0.122, 0.2)),
                                   (("UV", "MUV"), (0.2, 0.3)),
                                   (("UV", "NUV"), (0.3, 0.39)),
                                   (("Optical", "Violet"), (0.39, 0.45)),
                                   (("Optical", "Blue"), (0.45, 0.495)),
                                   (("Optical", "Green"), (0.495, 0.570)),
                                   (("Optical", "Yellow"), (0.57, 0.59)),
                                   (("Optical", "Orange"), (0.59, 0.62)),
                                   (("Optical", "Red"), (0.62, 0.75)),
                                   #(("Optical/IR", "Red/NIR"), (0.75, 1.0)),
                                   #(("IR", "NIR"), (1.0, 5.0)),
                                   (("IR", "NIR"), (0.75, 5.0)),
                                   (("IR", "MIR"), (5.0, 25.0)),
                                   (("IR", "FIR"), (25.0, 350.0)),
                                   (("Submm", "Submm"), (350.0, 1000.)),
                                   (("Radio", "Microwave"), (1000., 1e6))]) #, # 1 mm to 1 m
                                   # (("Radio", "VHF"), (1e6, 1e7)), # 1 m to 10 m
                                   # (("Radio", "HF"), (1e7, 1e8)), # 10 m to 100 m
                                   # (("Radio", "MF"), (1e8, 1e9)), # 100 m to 1 km
                                   # (("Radio", "LF"), (1e9, 1e10)), # 1 km to 10 km
                                   # (("Radio", "VLF"), (1e10, 1e11)), # 10 km to 100 km
                                   # (("Radio", "ULF"), (1e11, 1e12)),  # 100 km to 1,000 km
                                   # (("Radio", "SLF"), (1e11, 1e12)),  # 1,000 km to 10,000 km
                                   # (("Radio", "ELF"), (1e12, 1e13))]) # 1,000 km to 100,000 km

# -----------------------------------------------------------------

# Define names of the physical regimes
ionizing = "ionizing"
young = "young"
evolved = "evolved"
mix = "mix"
aromatic = "aromatic"
thermal = "thermal"

# -----------------------------------------------------------------

# The names of the physical regimes
physical_regimes = [ionizing, young, evolved, mix, aromatic, thermal]

# Define the ranges of the subgrids
physical_ranges = OrderedDict()

# FIRST: USED FOR CREATING WAVELENGTH GRIDS
#physical_ranges[sf] = QuantityRange(0.02, 0.085, unit="micron")
#physical_ranges[stellar] = QuantityRange(0.085, 3., unit="micron")
#physical_ranges[aromatic] = QuantityRange(3., 27., unit="micron")
#physical_ranges[thermal] = QuantityRange(27., 1000., unit="micron")
#physical_ranges[microwave] = QuantityRange(1000., 2000, unit="micron")

# NEW
physical_ranges[ionizing] = QuantityRange(0.01, 0.085, unit="micron")
physical_ranges[young] = QuantityRange(0.085, 0.39, unit="micron")
#physical_ranges[sf] = QuantityRange(0.01, 0.39, unit="micron")
physical_ranges[evolved] = QuantityRange(0.39, 1., unit="micron")
physical_ranges[mix] = QuantityRange(1., 7., unit="micron")
physical_ranges[aromatic] = QuantityRange(7., 27., unit="micron")
#physical_ranges[thermal] = QuantityRange(27., 1000., unit="micron")
physical_ranges[thermal] = QuantityRange(27., 2000., unit="micron")
#physical_ranges[microwave] = QuantityRange(1000., 2000, unit="micron")

# -----------------------------------------------------------------

def all_regimes(lower=True, divisions=True, subdivisions=True, as_dict=False):

    """
    This function ...
    :param lower:
    :param divisions:
    :param subdivisions:
    :param as_dict:
    :return:
    """

    from ...core.basics.containers import DefaultOrderedDict
    if as_dict: regimes = DefaultOrderedDict(list)
    else: regimes = set()

    # Loop over the dictionary
    for key in spectrum_wavelengths:

        # Get division
        if lower: division = key[0].lower()
        else: division = key[0]

        # Get subdivision
        if lower: subdivision = key[1].lower()
        else: subdivision = key[1]

        # Add
        if as_dict: regimes[division].append(subdivision)
        else:

            # Add both division and subdivision
            if divisions: regimes.add(division)
            if subdivisions: regimes.add(subdivision)

    # Return
    if as_dict: return regimes
    else: return list(regimes)

# -----------------------------------------------------------------

def regimes_in_range(wavelength_range, lower=True, divisions=True, subdivisions=True, as_dict=False):

    """
    This function ...
    :param wavelength_range:
    :param lower:
    :param divisions:
    :param subdivisions:
    :param as_dict:
    :return:
    """

    from ...core.basics.containers import DefaultOrderedDict
    if as_dict: regimes = DefaultOrderedDict(list)
    else: regimes = set()

    # Loop over the dictionary
    for key in spectrum_wavelengths:

        # In range?
        key_range = wavelength_range_for_regime(key[0], key[1])
        both_below = key_range.max < wavelength_range.min
        both_above = key_range.min > wavelength_range.max
        if both_below or both_above: continue

        # Get division
        if lower: division = key[0].lower()
        else: division = key[0]

        # Get subdivision
        if lower: subdivision = key[1].lower()
        else: subdivision = key[1]

        # Add
        if as_dict: regimes[division].append(subdivision)
        else:

            # Add both division and subdivision
            if divisions: regimes.add(division)
            if subdivisions: regimes.add(subdivision)

    # Return
    if as_dict: return regimes
    else: return list(regimes)

# -----------------------------------------------------------------

def physical_regimes_in_range(wavelength_range):

    """
    This function ...
    :param wavelength_range:
    :return:
    """

    # List for the names of the regimes
    regimes = []

    # Loop over the physical regimes
    for name in physical_regimes:

        # In range
        regime_range = physical_ranges[name]
        both_below = regime_range.max < wavelength_range.min
        both_above = regime_range.min > wavelength_range.max
        if both_below or both_above: continue

        # Add
        regimes.append(name)

    # Return the regime names
    return regimes

# -----------------------------------------------------------------

ranges = Map()
for key in spectrum_wavelengths:

    division = key[0].lower()
    subdivision = key[1].lower()

    if division not in ranges: ranges[division] = Map()
    ranges[division][subdivision] = Map()
    ranges[division][subdivision].min = spectrum_wavelengths[key][0] * u("micron")
    ranges[division][subdivision].max = spectrum_wavelengths[key][1] * u("micron")

# -----------------------------------------------------------------

def name_in_spectrum(wavelength):
    
    """
    This function ...
    :param wavelength:
    """

    # Get the value of the wavelength in micron
    wavelength_in_micron = wavelength.to("micron").value

    # Loop over all possible names
    for key in spectrum_wavelengths:

        # Get the range for this name and compare it to the given wavelength
        range = spectrum_wavelengths[key]
        if range[0] <= wavelength_in_micron <= range[1]: return key

    # Nothing?
    return None

# -----------------------------------------------------------------

def regime_for_wavelength(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return name_in_spectrum(wavelength)

# -----------------------------------------------------------------

def physical_regime_for_wavelength(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    # Loop over all possible physical regimes
    for name in physical_regimes:

        # Get the range
        regime_range = physical_ranges[name]

        # Check the range
        if wavelength in regime_range: return name

    # Nothing?
    return None

# -----------------------------------------------------------------

def regime_for_filter(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # Parse the filter
    from ...core.tools import types
    from ...core.filter.filter import parse_filter
    if types.is_string_type(fltr): fltr = parse_filter(fltr)

    # Return the regime
    return name_in_spectrum(fltr.wavelength)

# -----------------------------------------------------------------

def physical_regime_for_filter(fltr):

    """
    This function ...
    :param fltr:
    :return:
    """

    # Parse the filter
    from ...core.tools import types
    from ...core.filter.filter import parse_filter
    if types.is_string_type(fltr): fltr = parse_filter(fltr)

    # Return the physical regime
    return physical_regime_for_wavelength(fltr.wavelength)

# -----------------------------------------------------------------

def is_ionizing(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == ionizing

# -----------------------------------------------------------------

def is_young(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == young

# -----------------------------------------------------------------

def is_evolved(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == evolved

# -----------------------------------------------------------------

def is_mix(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == mix

# -----------------------------------------------------------------

def is_aromatic(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == aromatic

# -----------------------------------------------------------------

def is_thermal(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    return physical_regime_for_wavelength(wavelength) == thermal

# -----------------------------------------------------------------

def is_optical(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[0] == "Optical"

# -----------------------------------------------------------------

def is_uv(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[0] == "UV"

# -----------------------------------------------------------------

def is_ir(wavelength):

    """
    Thisj function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[0] == "IR"

# -----------------------------------------------------------------

def is_nir(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[1] == "NIR"

# -----------------------------------------------------------------

def is_mir(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[1] == "MIR"

# -----------------------------------------------------------------

def is_fir(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[1] == "FIR"

# -----------------------------------------------------------------

def is_submm(wavelength):

    """
    Thisn function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[0] == "Submm"

# -----------------------------------------------------------------

def is_fir_or_submm(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[1] == "FIR" or regime[0] == "Submm"

# -----------------------------------------------------------------

def is_radio(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    regime = regime_for_wavelength(wavelength)
    return regime[0] == "Radio"

# -----------------------------------------------------------------

def wavelength_range_for_regime(regime, subregime=None):

    """
    This function ...
    :param regime:
    :param subregime:
    :return:
    """

    # If subregime is defined
    if subregime is not None:

        if (regime, subregime) not in spectrum_wavelengths: raise ValueError("Invalid regime: (" + regime + ", " + subregime + ")")
        else:
            min_value = spectrum_wavelengths[(regime, subregime)][0]
            max_value = spectrum_wavelengths[(regime, subregime)][1]
            return QuantityRange(min_value, max_value, "micron")

    else:

        if regime in spectrum_wavelengths:
            min_value = spectrum_wavelengths[regime][0]
            max_value = spectrum_wavelengths[regime][1]
            return QuantityRange(min_value, max_value, "micron")

        min_value = None
        max_value = None

        # Search in the division names
        for key in spectrum_wavelengths:

            if key[0] == regime:

                min_value_key = spectrum_wavelengths[key][0]
                max_value_key = spectrum_wavelengths[key][1]

                #return QuantityRange(min_value, max_value, "micron")

                # Adapt the min and max
                if min_value is None or min_value_key < min_value: min_value = min_value_key
                if max_value is None or max_value_key > max_value: max_value = max_value_key

        if min_value is not None:

            return QuantityRange(min_value, max_value, "micron")

        # Break is not encountered (or return)
        #else:
        else: # no match found yet

            for key in spectrum_wavelengths:

                if key[1] == regime:

                    min_value = spectrum_wavelengths[key][0]
                    max_value = spectrum_wavelengths[key][1]
                    return QuantityRange(min_value, max_value, "micron")

            # Break not encountered
            else: raise ValueError("Invalid regime: " + regime)

# -----------------------------------------------------------------

def find_wavelength_range(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    keys = find_keys(regime)

    min_wavelength = None
    max_wavelength = None

    # Loop over the matching subregimes
    for regime, subregime in keys:

        # Get the wavelength range
        wavelength_range = wavelength_range_for_regime(regime, subregime)

        # Adjust
        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    #print(min_wavelength, max_wavelength)

    if min_wavelength is None: raise ValueError("Minimum wavelength is undefined")
    if max_wavelength is None: raise ValueError("Maximum wavelength is undefined")

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def find_single_wavelength_range(regime, return_string=False, only_subregime_string=False):

    """
    This function ...
    :param regime:
    :param return_string:
    :param only_subregime_string:
    :return:
    """

    regime, subregime = find_single_key(regime)

    if subregime is None:
        wavelength_range = wavelength_range_for_regime(regime)
        string = regime
    else:
        wavelength_range = wavelength_range_for_regime(regime, subregime)
        if only_subregime_string: string = subregime
        else: string = regime + "/" + subregime

    # Return the wavelength range
    if return_string: return wavelength_range, string
    else: return wavelength_range

# -----------------------------------------------------------------

def find_keys(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    keys = []

    # Search in the division names
    for key in spectrum_wavelengths:

        if key == regime: keys.append(key)
        if key[0].lower() == regime.lower(): keys.append(key)
        if key[1].lower() == regime.lower(): keys.append(key)

    # Return the keys
    return keys

# -----------------------------------------------------------------

def find_key_indices(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    indices = []

    for index, key in enumerate(spectrum_wavelengths.keys()):

        if key == regime: indices.append(index)
        if key[0].lower() == regime.lower(): indices.append(index)
        if key[1].lower() == regime.lower(): indices.append(index)

    return indices

# -----------------------------------------------------------------

def find_single_key(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    results = []

    for key in spectrum_wavelengths:

        if key == regime: return results.append(key)
        if key[0].lower() == regime.lower(): return key[0], None #results.append((key[0], None)) # here it's OK!
        if key[1].lower() == regime.lower(): results.append((key[0], key[1]))

    if len(results) == 0: raise ValueError("Invalid regime: " + regime)
    elif len(results) > 1: raise ValueError("Multiple matching keys for '" + str(regime) + "': " + ", ".join(str(r) for r in results))
    else: return results[0]

# -----------------------------------------------------------------

def find_first_key(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    for key in spectrum_wavelengths:

        if key == regime: return key
        if key[0].lower() == regime.lower(): return key
        if key[1].lower() == regime.lower(): return key

    raise ValueError("Invalid regime: " + regime)

# -----------------------------------------------------------------

def find_first_key_index(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    for index, key in enumerate(spectrum_wavelengths.keys()):

        if key == regime: return index
        if key[0].lower() == regime.lower(): return index
        if key[1].lower() == regime.lower(): return index

    raise ValueError("Invalid regime: " + regime)

# -----------------------------------------------------------------

def find_last_key(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    for key in reversed(spectrum_wavelengths.keys()):

        if key == regime: return key
        if key[0].lower() == regime.lower(): return key
        if key[1].lower() == regime.lower(): return key

    raise ValueError("Invalid regime: " + regime)

# -----------------------------------------------------------------

def find_last_key_index(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    for index, key in reversed(enumerate(spectrum_wavelengths.keys())):

        if key == regime: return key
        if key[0].lower() == regime.lower(): return key
        if key[1].lower() == regime.lower(): return key

    raise ValueError("Invalid regime: " + regime)

# -----------------------------------------------------------------

def wavelength_range_for_regimes(*regimes):

    """
    This function ...
    :param regimes:
    :return:
    """

    min_wavelength = None
    max_wavelength = None

    for regime in regimes:

        if "-" in regime:

            regime_a, regime_b = regime.split("-")
            regime_range = regimes_between_and_including(regime_a, regime_b)
            #print("regime range", regime_range)
            wavelength_range = wavelength_range_for_regimes(*regime_range)

        else:

            wavelength_range = wavelength_range_for_regime(regime)

        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def regimes_between(regime_a, regime_b):

    """
    This function ...
    :param regime_a:
    :param regime_b:
    :return:
    """

    triggered = True

    regimes = []

    for key in spectrum_wavelengths:

        # Already triggered, keep adding until regime_b
        if triggered:

            # Check if end
            if key[0] == regime_b or key[1] == regime_b: triggered = False

            # Not end
            else: regimes.append(key)

        # Not triggered yet, wait for regime_a
        else:

            # Check if begin
            if key[0] == regime_a or key[0] == regime_a: triggered = True
            else: pass

    # Return the regimes
    return regimes

# -----------------------------------------------------------------

def regimes_between_and_including(regime_a, regime_b):

    """
    Thisf unction ...
    :param regime_a:
    :param regime_b:
    :return:
    """

    return [regime_a] + regimes_between(regime_a, regime_b) + [regime_b]

# -----------------------------------------------------------------

def regimes_after(regime):

    """
    Thisf unction ...
    :param regime:
    :return:
    """

    regimes = []

    # Loop over the regimes in reversed order
    for key in reversed(spectrum_wavelengths):

        # Check if end
        if key[0] == regime or key[1] == regime: break

        # Not end
        else: regimes.append(key)

    # Return the regimes
    return list(reversed(regimes))

# -----------------------------------------------------------------

def regimes_after_and_including(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    return [find_first_key(regime)] + regimes_after(regime)

# -----------------------------------------------------------------

def wavelength_range_after(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    min_wavelength = None
    max_wavelength = None

    # Loop over the regimes after
    for regime_i in regimes_after(regime):

        # Get the wavelength range
        wavelength_range = wavelength_range_for_regime(regime_i)

        # Adapt min and max wavelength
        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def wavelength_range_after_and_including(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    min_wavelength = None
    max_wavelength = None

    # Loop over the regimes after and including
    for regime_i in regimes_after_and_including(regime):

        # Get the wavelength range
        wavelength_range = wavelength_range_for_regime(regime_i)

        # Adapt min and max wavelength
        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def wavelength_range_before(regime):

    """
    Thisnfunction ...
    :param regime:
    :return:
    """

    min_wavelength = None
    max_wavelength = None

    # Loop over the regimes before
    for regime_i in regimes_before(regime):

        # Get the wavelength range
        wavelength_range = wavelength_range_for_regime(regime_i)

        # Adapt min and max wavelength
        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def wavelength_range_before_and_including(regime):

    """
    Thisf unction ...
    :param regime:
    :return:
    """

    min_wavelength = None
    max_wavelength = None

    # Loop over the regimes before and including
    for regime_i in regimes_before_and_including(regime):

        # Get the wavelength range
        wavelength_range = wavelength_range_for_regime(regime_i)

        # Adapt min and max wavelength
        if min_wavelength is None or wavelength_range.min < min_wavelength: min_wavelength = wavelength_range.min
        if max_wavelength is None or wavelength_range.max > max_wavelength: max_wavelength = wavelength_range.max

    # Return the wavelength range
    return QuantityRange(min_wavelength, max_wavelength)

# -----------------------------------------------------------------

def regimes_before(regime):

    """
    Thisj function ...
    :param regime:
    :return:
    """

    regimes = []

    # Loop over the regimes
    for key in spectrum_wavelengths:

        # Check if end
        if key[0] == regime or key[1] == regime: break

        # Not end
        else: regimes.append(key)

    # Return the regimes
    return regimes

# -----------------------------------------------------------------

def regimes_before_and_including(regime):

    """
    This function ...
    :param regime:
    :return:
    """

    return regimes_before(regime) + [find_first_key(regime)]

# -----------------------------------------------------------------

def wavelength_in_regime(wavelength, regime):

    """
    This function ...
    :param wavelength:
    :param regime:
    :return:
    """

    wavelength_range = wavelength_range_for_regime(regime)
    return wavelength in wavelength_range

# -----------------------------------------------------------------

def wavelength_in_regimes(wavelength, regimes):

    """
    This function ...
    :param wavelength:
    :param regimes:
    :return:
    """

    wavelength_range = wavelength_range_for_regimes(*regimes)
    return wavelength in wavelength_range

# -----------------------------------------------------------------
