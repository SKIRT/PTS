#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

definition.add_section("wavelengths")
definition.sections["wavelengths"].add_optional("unit", "string", "the unit of the wavelengths", "micron")
definition.sections["wavelengths"].add_optional("min", "real", "the minimum wavelength", 0.1)
definition.sections["wavelengths"].add_optional("max", "real", "the maximum wavelength", 10)
definition.sections["wavelengths"].add_optional("npoints", "integer", "the number of wavelength points", 25)

definition.add_optional("packages", "real", "the number of photon packages per wavelength", 1e7)
definition.add_optional("remotes", "string_list", "the list of remote hosts on which to launch the simulations", ["nancy"])

# -----------------------------------------------------------------
