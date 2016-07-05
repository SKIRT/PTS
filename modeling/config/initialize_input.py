#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration(log_path="log")

# Add optional arguments
config.add_section("wavelengths")
config.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
config.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.09)
config.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 2000)
config.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 100)
config.sections["wavelengths"].add_optional("min_zoom", float, "the minimum wavelength of the zoomed-in grid", 1)
config.sections["wavelengths"].add_optional("max_zoom", float, "the maximum wavelength of the zoomed-in grid", 30)
config.sections["wavelengths"].add_optional("npoints_zoom", int, "the number of wavelength points in the zoomed-in grid", 100)

config.add_optional("packages", float, "the number of photon packages per wavelength", 2e5)
config.add_flag("selfabsorption", "enable dust self-absorption")
config.add_optional("dust_grid", str, "the type of dust grid to use (bintree, octtree or cartesian)", "bintree")

# -----------------------------------------------------------------
