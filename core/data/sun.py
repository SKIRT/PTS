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

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.data.sed import IntrinsicSED
from ...core.tools import introspection
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

# Determine the path to the Sun SED file
sun_sed_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "Sun", "SunSED.dat")

# -----------------------------------------------------------------

Lsun = 3.839e26 * Unit("W")                 # solar luminosity without solar neutrino radiation

# -----------------------------------------------------------------

# Effective wavelengths (in m)
eff_wavelengths = {"FUV": 152e-9 * Unit("m"),
                   "NUV": 231e-9 * Unit("m"),
                   "U": 365e-9 * Unit("m"),
                   "B": 445e-9 * Unit("m"),
                   "V": 551e-9 * Unit("m"),
                   "R": 658e-9 * Unit("m"),
                   "I": 806e-9 * Unit("m"),
                   "J": 1.22e-6 * Unit("m"),
                   "H": 1.63e-6 * Unit("m"),
                   "K": 2.19e-6 * Unit("m"),
                   "SDSS u": 354e-9 * Unit("m"),
                   "SDSS g": 477e-9 * Unit("m"),
                   "SDSS r": 623e-9 * Unit("m"),
                   "SDSS i": 763e-9 * Unit("m"),
                   "SDSS z": 913e-9 * Unit("m"),
                   "IRAC1": 3.56e-6 * Unit("m"),
                   "IRAC2": 4.51e-6 * Unit("m"),
                   "WISE1": 3.35e-9 * Unit("m"),
                   "WISE2": 4.60e-6 * Unit("m")}

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
        self.sed = IntrinsicSED.from_file(sun_sed_path, skiprows=4) # wavelength in micron, luminosity in W/micron

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
        return Unit("Lsun", self.total_luminosity())

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, filter, unit="W/micron"):

        """
        This function ...
        :param filter:
        :param unit:
        :return:
        """

        # Convole the Sun SED over the filter transmission curve
        luminosity = filter.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.luminosities(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * Unit("W/micron")

        # Return the luminosity
        return luminosity.to(unit)

    # -----------------------------------------------------------------

    def luminosity_for_filter_as_unit(self, filter):

        """
        This function ...
        :param filter:
        :return:
        """

        # Determine a name for this unit
        unit_name = "Lsun_" + filter.band

        # Create and return the new unit
        return Unit(unit_name, self.luminosity_for_filter(filter))

# -----------------------------------------------------------------
