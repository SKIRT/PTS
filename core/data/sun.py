#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.sun Contains the Sun class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.data.sed import SED
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

# Determine the path to the Sun SED file
sun_sed_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "Sun", "SunSED.dat")

# -----------------------------------------------------------------

Lsun = 3.839e26 * u("W")                 # solar luminosity without solar neutrino radiation

# -----------------------------------------------------------------

# Effective wavelengths (in m)
eff_wavelengths = {"FUV": 152e-9 * u("m"),
                   "NUV": 231e-9 * u("m"),
                   "U": 365e-9 * u("m"),
                   "B": 445e-9 * u("m"),
                   "V": 551e-9 * u("m"),
                   "R": 658e-9 * u("m"),
                   "I": 806e-9 * u("m"),
                   "J": 1.22e-6 * u("m"),
                   "H": 1.63e-6 * u("m"),
                   "K": 2.19e-6 * u("m"),
                   "SDSS u": 354e-9 * u("m"),
                   "SDSS g": 477e-9 * u("m"),
                   "SDSS r": 623e-9 * u("m"),
                   "SDSS i": 763e-9 * u("m"),
                   "SDSS z": 913e-9 * u("m"),
                   "IRAC1": 3.56e-6 * u("m"),
                   "IRAC2": 4.51e-6 * u("m"),
                   "WISE1": 3.35e-9 * u("m"),
                   "WISE2": 4.60e-6 * u("m")}

# -----------------------------------------------------------------

class Sun(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Load the intrinsic SED of the sun
        self.sed = SED.from_text_file(sun_sed_path, wavelength_unit="micron", photometry_unit="W/micron", skiprows=4) # wavelength in micron, luminosity in W/micron

        # The total luminosity
        self.luminosity = Lsun

    # -----------------------------------------------------------------

    def total_luminosity(self, unit="W"):

        """
        This function ...
        :param unit:
        :return:
        """

        return Lsun.to(unit)

    # -----------------------------------------------------------------

    def total_luminosity_as_unit(self):

        """
        This function ...
        :return:
        """

        # Create and return the new unit
        return u("Lsun", self.total_luminosity())

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, fltr, unit="W/micron"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        # Convole the Sun SED over the filter transmission curve
        luminosity = fltr.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.photometry(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * u("W/micron")

        # Return the luminosity
        return luminosity.to(unit)

    # -----------------------------------------------------------------

    def luminosity_for_filter_as_unit(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Determine a name for this unit
        unit_name = "Lsun_" + fltr.band

        # Create and return the new unit
        from astropy.units import Unit
        return Unit(unit_name, represents=self.luminosity_for_filter(fltr))

# -----------------------------------------------------------------
