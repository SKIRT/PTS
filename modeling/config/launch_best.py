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
definition.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
definition.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.05)
definition.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 1000)
definition.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 300)
definition.sections["wavelengths"].add_optional("min_zoom", float, "the minium wavelength of the zoomed-in grid", 1)
definition.sections["wavelengths"].add_optional("max_zoom", float, "the maximum wavelength of the zoomed-in grid", 30)
definition.sections["wavelengths"].add_optional("npoints_zoom", int, "the number of wavelength points in the zoomed-in grid", 300)

definition.add_optional("packages", float, "the number of photon packages per wavelength", 1e6)
definition.add_optional("selfabsorption", bool, "whether self-absorption should be enabled", True)
definition.add_optional("remote", str, "the remote host on which to launch the simulations", "nancy")

# -----------------------------------------------------------------
