#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.wavelengths Provides ...

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import astropy.units as u

# Import the relevant PTS classes and modules
from ...core.basics.map import Map

# -----------------------------------------------------------------

spectrum_wavelengths = {("UV", "EUV"): (0.01, 0.121),
                        ("UV", "Lyman-alpha"): (0.121, 0.122),
                        ("UV", "FUV"): (0.122, 0.2),
                        ("UV", "MUV"): (0.2, 0.3),
                        ("UV", "NUV"): (0.3, 0.39),
                        ("Optical", "Violet"): (0.39, 0.45),
                        ("Optical", "Blue"): (0.45, 0.495),
                        ("Optical", "Green"): (0.495, 0.570),
                        ("Optical", "Yellow"): (0.57, 0.59),
                        ("Optical", "Orange"): (0.59, 0.62),
                        ("Optical", "Red"): (0.62, 0.75),
                        ("Optical/IR", "Red/NIR"): (0.75, 1.0),
                        #("IR", "NIR"): (0.75, 5.0),
                        ("IR", "NIR"): (1.0, 5.0),
                        ("IR", "MIR"): (5.0, 25.0),
                        ("IR", "FIR"): (25.0, 350.0),
                        ("Submm", "Submm"): (350.0, 1000.),
                        ("Radio", "Radio"): (1000., 1e9)}

# -----------------------------------------------------------------

ranges = Map()
for key in spectrum_wavelengths:

    division = key[0].lower()
    subdivision = key[1].lower()

    if division not in ranges: ranges[division] = Map()
    ranges[division][subdivision] = Map()
    ranges[division][subdivision].min = spectrum_wavelengths[key][0] * u.Unit("micron")
    ranges[division][subdivision].max = spectrum_wavelengths[key][1] * u.Unit("micron")

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

# -----------------------------------------------------------------
