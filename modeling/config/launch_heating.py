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

config.add_section("wavelengths")
config.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
config.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.1)
config.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 10)
config.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 25)

config.add_optional("packages", float, "the number of photon packages per wavelength", 1e7)
config.add_optional("remotes", "string_list", "the list of remote hosts on which to launch the simulations", ["nancy"])

# -----------------------------------------------------------------
