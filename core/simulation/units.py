#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.units Working with SKIRT output units.
#
# An instance of the SkirtUnits class in this module provides support for working with SKIRT input/output units.
# The constructor arguments specify the name of a SKIRT unit system (SI, stellar, or extragalactic units)
# and the flux style (neutral, wavelength or frequency) to set the default units for physical quantities.
# The instance then offers functions to convert a physical quantity from its default unit to some specified unit
# (or between two specified units).

# -----------------------------------------------------------------

# Import standard modules
import types
import numpy as np

# -----------------------------------------------------------------
#  SkirtUnits class
# -----------------------------------------------------------------

## An instance of the SkirtUnits class represents a particular SKIRT input/output unit system, specified at
# construction through the name of the unit system (SI, stellar, or extragalactic units) and the flux style
# (neutral, wavelength or frequency). Based on this information, the object knows the output units used by SKIRT
# for a series of supported physical quantities. This allows converting a value extracted from SKIRT output to
# some specified unit without having to know the actual output unit.
#
# The SkirtUnits class supports the following physical quantities and unit specifiers:
#
#| Physical Quantity | Flux Style | Units
#|-------------------|------------|-------
#| length, distance, wavelength | | A, nm, micron, mm, cm, m, km, AU, pc, kpc, Mpc
#| volume | | m3, AU3, pc3
#| mass | | g, kg, Msun
#| luminosity | | W, Lsun
#| luminositydensity | wavelength | W/m, W/micron, Lsun/micron
#| luminositydensity | frequency | W/Hz, erg/s/Hz, Lsun/Hz
#| fluxdensity | neutral | W/m2
#| fluxdensity | wavelength | W/m3, W/m2/micron
#| fluxdensity | frequency | W/m2/Hz, Jy, mJy, MJy, erg/s/cm2/Hz
#| surfacebrightness | neutral | W/m2/sr, W/m2/arcsec2
#| surfacebrightness | wavelength | W/m3/sr, W/m2/micron/sr, W/m2/micron/sr, W/m2/micron/arcsec2
#| surfacebrightness | frequency | W/m2/Hz/sr, W/m2/Hz/arcsec2, Jy/sr, Jy/arcsec2, MJy/sr, MJy/arcsec2
#
# Flux style 'neutral' indicates \f$\lambda F_\lambda = \nu F_\nu\f$; 'wavelength' indicates
# \f$F_\lambda\f$; and 'frequency' indicates \f$F_\nu\f$.

class SkirtUnits:

    ## The constructor accepts the name of the SKIRT unit system ('SI', 'stellar', or 'extragalactic')
    # and the flux style ('neutral', 'wavelength' or 'frequency') to be represented by this instance.
    # The specified strings are case-insensitive, and any portion beyond the recognized names is ignored.
    # Based on this information, it initializes the default SKIRT units for a series of supported
    # physical quantities.
    #
    def __init__(self, unitsystem, fluxstyle):
        unitsystem = unitsystem.lower()
        fluxstyle = fluxstyle.lower()

        if unitsystem.startswith('si'):
            self._defaultunit = { 'length': 'm', 'distance': 'm', 'wavelength': 'm',
                                  'volume': 'm3', 'mass': 'kg',
                                  'luminosity': 'W', 'luminositydensity': 'W/m' }
            if fluxstyle.startswith('neutral'):
                self._defaultunit.update(fluxdensity='W/m2', surfacebrightness='W/m2/sr')
            elif fluxstyle.startswith('wavelength'):
                self._defaultunit.update(fluxdensity='W/m3', surfacebrightness='W/m3/sr')
            elif fluxstyle.startswith('frequency'):
                self._defaultunit.update(fluxdensity='W/m2/Hz', surfacebrightness='W/m2/Hz/sr')
            else:
                raise ValueError("Unsupported flux style: " + fluxstyle)

        elif unitsystem.startswith('stellar'):
            self._defaultunit = { 'length': 'AU', 'distance': 'pc', 'wavelength': 'micron',
                                  'volume': 'AU3', 'mass': 'Msun',
                                  'luminosity': 'Lsun', 'luminositydensity': 'Lsun/micron' }
            if fluxstyle.startswith('neutral'):
                self._defaultunit.update(fluxdensity='W/m2', surfacebrightness='W/m2/arcsec2')
            elif fluxstyle.startswith('wavelength'):
                self._defaultunit.update(fluxdensity='W/m2/micron', surfacebrightness='W/m2/micron/arcsec2')
            elif fluxstyle.startswith('frequency'):
                self._defaultunit.update(fluxdensity='Jy', surfacebrightness='MJy/sr')
            else:
                raise ValueError("Unsupported flux style: " + fluxstyle)

        elif unitsystem.startswith('extragalactic'):
            self._defaultunit = { 'length': 'pc', 'distance': 'Mpc', 'wavelength': 'micron',
                                  'volume': 'pc3', 'mass': 'Msun',
                                  'luminosity': 'Lsun', 'luminositydensity': 'Lsun/micron' }
            if fluxstyle.startswith('neutral'):
                self._defaultunit.update(fluxdensity='W/m2', surfacebrightness='W/m2/arcsec2')
            elif fluxstyle.startswith('wavelength'):
                self._defaultunit.update(fluxdensity='W/m2/micron', surfacebrightness='W/m2/micron/arcsec2')
            elif fluxstyle.startswith('frequency'):
                self._defaultunit.update(fluxdensity='Jy', surfacebrightness='MJy/sr')
            else:
                raise ValueError("Unsupported flux style: " + fluxstyle)

        else:
            raise ValueError("Unsupported unit system: " + unitsystem)


    ## This function performs unit conversion for a specified value (or for a numpy array of values). The first
    # argument specifies the value to be converted. This can be a number, a numpy array of numbers (in which case
    # the conversion is performed for each array element), or a string representing a number optionally followed
    # by a unit specifier (e.g. "0.76 micron"). The second argument specifies the target unit. In addition the
    # function accepts the following optional arguments:
    #  - \em from_unit: specifies the units in which the incoming value is expressed.
    #  - \em quantity: specifies the physical quantity of the incoming value; this is used to determine the appropriate
    #    default unit in case the \em from_unit argument is missing (and the value is not a string including units).
    #  - \em wavelength: provides the wavelength, in micron, for converting between flux styles; this argument can be
    #    omitted for any other conversions; if \em value is a numpy array, \em wavelength can be an array of the same
    #    length (one wavelength per flux value), or it can be a single number (in which case all fluxes are considered
    #    to be at the same wavelength).
    #
    # The unit of the incoming value is determined using three mechanisms in the following order:
    #  - if the value is a string with two segments, the second segement determines the unit.
    #  - otherwise, if \em from_unit is specified (and not None), its value determines the unit.
    #  - otherwise, the default SKIRT unit corresponding to the specified \em quantity is used.
    #
    def convert(self, value, to_unit, from_unit=None, quantity=None, wavelength=None):

        # if the value is a string, it may include a unit specifier that overrides the from_unit argument
        if isinstance(value, types.StringTypes):
            segments = value.split()
            if len(segments) == 2:
                value = float(segments[0])
                from_unit = segments[1]
            elif len(segments) == 1:
                value = float(segments[0])
            else:
                raise ValueError("Invalid value/unit string")

        # if the from_unit has not been specified, use the default for the specified quantity
        if from_unit==None:
            from_unit = self._defaultunit[quantity]

        # skip the conversion if the units are identical
        if from_unit == to_unit:
            return value

        # perform straightforward conversion between units of the same physical quantity
        from_quantity = _quantity[from_unit]
        to_quantity = _quantity[to_unit]
        if from_quantity == to_quantity:
            return value * (_conversion[from_unit]/_conversion[to_unit])

        # perform conversion between styles of flux density or surface brightness
        if ('fluxdensity' in from_quantity and 'fluxdensity' in to_quantity) or  \
           ('luminositydensity' in from_quantity and 'luminositydensity' in to_quantity) or  \
           ('surfacebrightness' in from_quantity and 'surfacebrightness' in to_quantity):

            # convert to/from SI units within the respective flux styles
            flux = value * (_conversion[from_unit]/_conversion[to_unit])

            # convert between flux styles (convert specified wavelength from micron to m)
            wave = wavelength * 1e-6
            if 'wavelength' in from_quantity: flux *= wave
            elif 'frequency' in from_quantity: flux *= _c/wave
            if 'wavelength' in to_quantity: flux *= 1./wave
            elif 'frequency' in to_quantity: flux *= wave/_c
            return flux

        else:
            raise ValueError("Can't convert from " + from_unit + " to " + to_unit)

    ## This function returns the absolute AB magnitude corresponding to a given flux density and distance
    # from the source. The units in which these values are expressed can be explicitly specified. If not,
    # the default units for respectively flux density and distance are used instead. If the flux density
    # is expressed per unit of frequency, the \em wavelength argument may be omitted. Otherwise, the
    # wavelength is used to convert between flux styles.
    #
    # Given a flux density \f$F_\nu\f$, measured in ergs per second per square cm per Hz, the corresponding
    # AB magnitude is defined as \f$\text{AB}=-2.5\log_{10} F_\nu -48.60\f$. The resulting apparent magnitude
    # is converted to the absolute magnitude using the standard formula \f$M=m-5\log_{10}d^{(pc)}+5\f$.
    def absolutemagnitude(self, fluxdensity, distance, fluxdensity_unit=None, distance_unit=None, wavelength=None):
        fluxdensity = self.convert(fluxdensity, to_unit='erg/s/cm2/Hz', from_unit=fluxdensity_unit,
                                   quantity='fluxdensity', wavelength=wavelength)
        distance = self.convert(distance, to_unit='pc', from_unit=distance_unit, quantity='distance')
        apparent = -2.5*np.log10(fluxdensity) - 48.60
        absolute = apparent - 5*np.log10(distance) + 5
        return absolute

    ## This function returns the luminosity density corresponding to a given flux density and distance
    # from the source. The units in which these values are expressed can be explicitly specified. If not,
    # the default units for respectively flux density and distance are used instead. The units for the
    # returned luminosity must be specified (there is no default). If both the flux density and the
    # luminosity density are expressed in the same style (per unit of frequency or per unit of wavelength),
    # the \em wavelength argument may be omitted. Otherwise, the wavelength is used to convert between styles.
    def luminosityforflux(self, fluxdensity, distance, luminositydensity_unit,
                                fluxdensity_unit=None, distance_unit=None, wavelength=None):
        if 'wavelength' in _quantity[luminositydensity_unit]:
            flux_si = 'W/m3'
            lumi_si = 'W/m'
        else:
            flux_si = 'W/m2/Hz'
            lumi_si = 'W/Hz'
        fluxdensity = self.convert(fluxdensity, to_unit=flux_si, from_unit=fluxdensity_unit,
                                   quantity='fluxdensity', wavelength=wavelength)
        distance = self.convert(distance, to_unit='m', from_unit=distance_unit, quantity='distance')
        luminosity = 4.*np.pi * distance*distance * fluxdensity
        return self.convert(luminosity, to_unit=luminositydensity_unit, from_unit=lumi_si)

# -----------------------------------------------------------------
#  Private conversion facilities
# -----------------------------------------------------------------

# --- fundamental physical and astronomical constants ---

_c = 2.99792458e8     # light speed in m/s
_AU = 1.49597871e11   # astronomical unit in m
_pc = 3.08567758e16   # parsec in m
_Lsun = 3.839e26      # solar bolometric luminosity in W (without solar neutrino radiation)
_Msun = 1.9891e30     # solar mass in kg
_arcsec2 = 2.350443053909789e-11  # solid angle of 1 square arc second in steradian

# --- library providing the physical quantity corresponding to each supported unit ---

# key: unit; value: physical quantity for this unit
_quantity = { 'A': 'length', 'nm': 'length', 'micron': 'length', 'mm': 'length', 'cm': 'length',
              'm': 'length', 'km': 'length', 'AU': 'length', 'kpc': 'length', 'pc': 'length', 'Mpc': 'length',
              'm3': 'volume', 'AU3': 'volume', 'pc3': 'volume',
              'g': 'mass', 'kg': 'mass', 'Msun': 'mass',
              'W': 'luminosity', 'Lsun': 'luminosity',
              'W/m': 'wavelengthluminositydensity',
              'W/micron': 'wavelengthluminositydensity',
              'Lsun/micron': 'wavelengthluminositydensity',
              'W/Hz': 'frequencyluminositydensity',
              'erg/s/Hz': 'frequencyluminositydensity',
              'Lsun/Hz': 'frequencyluminositydensity',
              'W/m2': 'neutralfluxdensity',
              'W/m2/sr': 'neutralsurfacebrightness',
              'W/m2/arcsec2': 'neutralsurfacebrightness',
              'W/m3': 'wavelengthfluxdensity',
              'W/m2/micron': 'wavelengthfluxdensity',
              'W/m3/sr': 'wavelengthsurfacebrightness',
              'W/m2/micron/sr': 'wavelengthsurfacebrightness',
              'W/m2/micron/arcsec2': 'wavelengthsurfacebrightness',
              'W/m2/Hz': 'frequencyfluxdensity',
              'Jy': 'frequencyfluxdensity',
              'mJy': 'frequencyfluxdensity',
              'MJy': 'frequencyfluxdensity',
              'erg/s/cm2/Hz': 'frequencyfluxdensity',
              'W/m2/Hz/sr': 'frequencysurfacebrightness',
              'W/m2/Hz/arcsec2': 'frequencysurfacebrightness',
              'Jy/sr': 'frequencysurfacebrightness',
              'Jy/arcsec2': 'frequencysurfacebrightness',
              'MJy/sr': 'frequencysurfacebrightness',
              'MJy/arcsec2': 'frequencysurfacebrightness'
            }

# --- library providing the conversion factor to SI units for each supported unit ---

# key: unit; value: conversion factor to corresponding SI unit
_conversion = { 'A': 1e-10, 'nm': 1e-9, 'micron': 1e-6, 'mm': 1e-3, 'cm': 1e-2,
                'm': 1., 'km': 1e3, 'AU': _AU, 'pc': _pc, 'kpc': 1e3*_pc, 'Mpc': 1e6*_pc,
                'm3': 1., 'AU3': _AU**3, 'pc3': _pc**3,
                'g': 1e-3, 'kg': 1., 'Msun': _Msun,
                'W': 1., 'Lsun': _Lsun,
                'W/m': 1., 'W/micron': 1e6, 'Lsun/micron': _Lsun*1e6,
                'W/Hz': 1., 'Lsun/Hz': _Lsun, "erg/s/Hz": 1e-7,
                'W/m2': 1.,
                'W/m2/sr': 1.,
                'W/m2/arcsec2': 1./_arcsec2,
                'W/m3': 1.,
                'W/m2/micron': 1e6,
                'W/m3/sr': 1.,
                'W/m2/micron/sr': 1e6,
                'W/m2/micron/arcsec2': 1e6/_arcsec2,
                'W/m2/Hz': 1.,
                'Jy': 1e-26,
                'mJy': 1e-29,
                'MJy': 1e-20,
                'erg/s/cm2/Hz': 1e-3,
                'W/m2/Hz/sr': 1.,
                'W/m2/Hz/arcsec2': 1./_arcsec2,
                'Jy/sr': 1e-26,
                'Jy/arcsec2': 1e-26/_arcsec2,
                'MJy/sr': 1e-20,
                'MJy/arcsec2': 1e-20/_arcsec2
              }

# -----------------------------------------------------------------
