#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.normalizations Show the normalizations of stellar components in a ski file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.skifile import LabeledSkiFile
from pts.core.filter.filter import Filter
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Ski path
definition.add_required("ski", "file_path", "path of the ski file")

# Optional
definition.add_optional("flux", "photometric_unit", "also show flux in a particular unit")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("normalizations", definition, "Show the normalizations of stellar components in a ski file")

# -----------------------------------------------------------------

# Load the ski file
ski = LabeledSkiFile(config.ski)

# -----------------------------------------------------------------

# Get stellar component IDs
stellar_component_ids = ski.get_stellar_component_ids()

# -----------------------------------------------------------------

wavelengths = defaultdict(dict)
filters = defaultdict(dict)

# -----------------------------------------------------------------

# Loop over the stellar components
for component_id in stellar_component_ids:

    # Get normalization filter/wavelength and luminosity
    filter_or_wavelength = ski.get_stellar_component_normalization_wavelength_or_filter(component_id)
    luminosity = ski.get_stellar_component_luminosity(component_id, return_wavelength=False)

    # Add to dict
    if isinstance(filter_or_wavelength, Filter): filters[filter_or_wavelength][component_id] = luminosity
    else: wavelengths[filter_or_wavelength][component_id] = luminosity

# -----------------------------------------------------------------

nwavelengths = len(wavelengths)
nfilters = len(filters)
has_wavelengths = nwavelengths > 0
has_filters = nfilters > 0

# -----------------------------------------------------------------

# Get distance (returns None if couldn't define)
distance = ski.get_distance(return_none=True)

# -----------------------------------------------------------------

if has_wavelengths:

    print("")
    print(fmt.underlined + fmt.blue + "Wavelength normalizations:" + fmt.reset)
    print("")

    # Loop over the wavelengths
    for wavelength in wavelengths:

        print(fmt.bold + fmt.green + tostr(wavelength) + fmt.reset + ":")
        print("")

        total_luminosity = 0.

        # Loop over the components
        for component_id in wavelengths[wavelength]:

            luminosity = wavelengths[wavelength][component_id]
            if config.flux is not None: flux = luminosity.to(config.flux, distance=distance, wavelength=wavelength)
            else: flux = None
            if flux is not None: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity) + " (" + tostr(flux) + ")")
            else: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity))
            total_luminosity += luminosity

        if config.flux is not None: total_flux = total_luminosity.to(config.flux, distance=distance, wavelength=wavelength)
        else: total_flux = None

        print("")
        if total_flux is not None: print(" - TOTAL: " + tostr(total_luminosity) + " (" + tostr(total_flux) + ")")
        else: print(" - TOTAL: " + tostr(total_luminosity))
        print("")

# -----------------------------------------------------------------

if has_filters:

    print("")
    print(fmt.underlined + fmt.blue + "Filter normalizations:" + fmt.reset)
    print("")

    # Loop over the filters
    for fltr in filters:

        print(fmt.bold + fmt.green + str(fltr) + fmt.reset + ":")
        print("")

        # Get wavelength
        wavelength = fltr.wavelength

        total_luminosity = 0.

        # Loop over the components
        for component_id in filters[fltr]:

            luminosity = filters[fltr][component_id]
            if config.flux is not None: flux = luminosity.to(config.flux, distance=distance, wavelength=wavelength)
            else: flux = None
            if flux is not None: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity) + " (" + tostr(flux) + ")")
            else: print(" - " + fmt.bold + component_id + fmt.reset + ": " + tostr(luminosity))
            total_luminosity += luminosity

        if config.flux is not None: total_flux = total_luminosity.to(config.flux, distance=distance, wavelength=wavelength)
        else: total_flux = None

        print("")
        if total_flux is not None: print(" - TOTAL: " + tostr(total_luminosity) + " (" + tostr(total_flux) + ")")
        else: print(" - TOTAL: " + tostr(total_luminosity))
        print("")

# -----------------------------------------------------------------
