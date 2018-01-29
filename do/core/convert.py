#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.convert Convert a quantity from one unit to another.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments, prompt_yn
from pts.core.tools.stringify import tostr
from pts.core.units.parsing import is_photometric_unit, is_photometric_quantity, parse_unit, parse_quantity
from pts.core.units.parsing import possible_physical_types, possible_density_flags, possible_brightness_flags

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition(write_config=False)

# Quantity to convert and unit to convert to
definition.add_required("quantity", "string", "quantity to convert to another unit")
definition.add_required("unit", "string", "unit to convert the quantity to")

# Extra information
definition.add_optional("distance", "length_quantity", "distance")
definition.add_optional("wavelength", "length_quantity", "wavelength")
definition.add_optional("frequency", "frequency_quantity", "frequency")
definition.add_optional("pixelscale", "pixelscale", "pixelscale")
definition.add_optional("solid_angle", "solid_angle", "solid angle")
definition.add_optional("filter", "filter", "filter")

# Create the configuration
config = parse_arguments("convert", definition, "Convert a quantity from one unit to another")

# -----------------------------------------------------------------

# Check quantity
if is_photometric_quantity(config.quantity):

    #physical_types = possible_physical_types(config.quantity)
    #print(physical_types)

    # Density?
    density_flags = possible_density_flags(config.quantity)
    if len(density_flags) == 1: density = density_flags[0]
    elif len(density_flags) == 2: density = prompt_yn("density", "is this quantity (" + config.quantity + ") a spectral density?")
    else: raise RuntimeError("Something went wrong")

    # Brightness?
    brightness_flags = possible_brightness_flags(config.quantity)
    if len(brightness_flags) == 1: brightness = brightness_flags[0]
    elif len(brightness_flags) == 2: brightness = prompt_yn("brightness", "is this quantity (" + config.quantity + ") a surface brightness?")
    else: raise RuntimeError("Something went wrong")

    # Parse
    quantity = parse_quantity(config.quantity, density=density, brightness=brightness, density_strict=True, brightness_strict=True)

# Not a photometric quantity
else: quantity = parse_quantity(config.quantity)

# -----------------------------------------------------------------

# Check unit
if is_photometric_unit(config.unit):

    #physical_types = possible_physical_types(config.unit)
    #print(physical_types)

    # Density?
    density_flags = possible_density_flags(config.unit)
    if len(density_flags) == 1: density = density_flags[0]
    elif len(density_flags) == 2: density = prompt_yn("density", "is this unit (" + config.unit + ") a spectral density?")
    else: raise RuntimeError("Something went wrong")

    # Brightness?
    brightness_flags = possible_brightness_flags(config.unit)
    if len(brightness_flags) == 1: brightness = brightness_flags[0]
    elif len(brightness_flags) == 2: brightness = prompt_yn("brightness", "is this unit (" + config.unit + ") a surface brightness?")
    else: raise RuntimeError("Something went wrong")

    # Parse
    unit = parse_unit(config.unit, density=density, brightness=brightness, density_strict=True, brightness_strict=True)

# Not a photometric unit
else: unit = parse_unit(config.unit)

# -----------------------------------------------------------------

# Set conversion info
conversion_info = dict()
if config.distance is not None: conversion_info["distance"] = config.distance
if config.wavelength is not None: conversion_info["wavelength"] = config.wavelength
if config.frequency is not None: conversion_info["frequency"] = config.frequency
if config.pixelscale is not None: conversion_info["pixelscale"] = config.pixelscale
if config.solid_angle is not None: conversion_info["solid_angle"] = config.solid_angle
if config.filter is not None: conversion_info["fltr"] = config.filter

# -----------------------------------------------------------------

# Convert
converted = quantity.to(unit, **conversion_info)

# -----------------------------------------------------------------

# Show
print(tostr(converted))

# -----------------------------------------------------------------
