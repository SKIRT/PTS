#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.regimes Show wavelength regimes

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.tools.wavelengths import all_regimes, regimes_in_range, wavelength_range_for_regime
from pts.magic.tools.wavelengths import physical_regimes, physical_ranges, physical_regimes_in_range
from pts.core.tools import formatting as fmt
from pts.core.filter.broad import get_filters as get_broad_band_filters
from pts.core.tools.stringify import tostr, stringify_list_fancy

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()

# Wavelength range
definition.add_positional_optional("wavelength_range", "quantity_range", "wavelength range")

# Flags
definition.add_flag("filters", "show filters in regimes", False)
definition.add_flag("physical", "use physical regimes (star formation, stellar emission, aromatic features, dust thermal emission, microwave)", False)

# Parse the command line arguments
config = parse_arguments("regimes", definition)

# -----------------------------------------------------------------

# Physical
if config.physical:

    # Get the list of regimes
    if config.wavelength_range is not None: regimes = physical_regimes_in_range(config.wavelength_range)
    else: regimes = physical_regimes

    print("")
    for name in regimes:

        print(fmt.bold + fmt.blue + name + fmt.reset)
        print("")

        # Get wavelength range
        wavelength_range = physical_ranges[name]

        # Show minimum wavelength
        if config.wavelength_range is not None and wavelength_range.min < config.wavelength_range.min: print("   - minimum wavelength: " + fmt.red + tostr(wavelength_range.min) + fmt.reset)
        else: print("   - minimum wavelength: " + tostr(wavelength_range.min))

        # Show maximum wavelength
        if config.wavelength_range is not None and wavelength_range.max > config.wavelength_range.max: print("   - maximum wavelength: " + fmt.red + tostr(wavelength_range.max) + fmt.reset)
        else: print("   - maximum wavelength: " + tostr(wavelength_range.max))

        if config.filters:

            filters = get_broad_band_filters(wavelength_range.min, wavelength_range.max)

            if config.wavelength_range is not None:
                colour = "red"
                indices_outside_range = [index for index in range(len(filters)) if filters[index].wavelength not in config.wavelength_range]
            else: colour = indices_outside_range = None

            if len(filters) > 0:
                print("   - filters:")
                print("")
                print(stringify_list_fancy(filters, lines_prefix="       ", colour=colour, colour_indices=indices_outside_range)[1])

        print("")

# -----------------------------------------------------------------

# Not physical
else:

    # Get the list of regimes
    if config.wavelength_range is not None: regimes = regimes_in_range(config.wavelength_range, lower=False, divisions=False, as_dict=True)
    else: regimes = all_regimes(lower=False, divisions=False, as_dict=True)

    print("")
    for division in regimes:

        print(fmt.bold + fmt.blue + division + fmt.reset)
        print("")

        # Loop over the subdivision
        for subdivision in regimes[division]:

            print("  " + fmt.bold + fmt.green + subdivision + fmt.reset)
            print("")

            # Get wavelength range
            wavelength_range = wavelength_range_for_regime(division, subdivision)

            # Show minimum wavelength
            if config.wavelength_range is not None and wavelength_range.min < config.wavelength_range.min: print("   - minimum wavelength: " + fmt.red + tostr(wavelength_range.min) + fmt.reset)
            else: print("   - minimum wavelength: " + tostr(wavelength_range.min))

            # Show maximum wavelength
            if config.wavelength_range is not None and wavelength_range.max > config.wavelength_range.max: print("   - maximum wavelength: " + fmt.red + tostr(wavelength_range.max) + fmt.reset)
            else: print("   - maximum wavelength: " + tostr(wavelength_range.max))

            if config.filters:

                filters = get_broad_band_filters(wavelength_range.min, wavelength_range.max)

                if config.wavelength_range is not None:
                    colour = "red"
                    indices_outside_range = [index for index in range(len(filters)) if filters[index].wavelength not in config.wavelength_range]
                else: colour = indices_outside_range = None

                if len(filters) > 0:
                    print("   - filters:")
                    print("")
                    print(stringify_list_fancy(filters, lines_prefix="       ", colour=colour, colour_indices=indices_outside_range)[1])

            print("")

# -----------------------------------------------------------------
