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

# Add optional arguments
definition.add_section("wavelengths")
definition.sections["wavelengths"].add_optional("unit", str, "the unit of the wavelengths", "micron")
definition.sections["wavelengths"].add_optional("min", float, "the minimum wavelength", 0.09)
definition.sections["wavelengths"].add_optional("max", float, "the maximum wavelength", 2000)
definition.sections["wavelengths"].add_optional("npoints", int, "the number of wavelength points", 100)
definition.sections["wavelengths"].add_optional("min_zoom", float, "the minimum wavelength of the zoomed-in grid", 1)
definition.sections["wavelengths"].add_optional("max_zoom", float, "the maximum wavelength of the zoomed-in grid", 30)
definition.sections["wavelengths"].add_optional("npoints_zoom", int, "the number of wavelength points in the zoomed-in grid", 100)

definition.add_optional("packages", float, "the number of photon packages per wavelength", 2e5)
definition.add_flag("selfabsorption", "enable dust self-absorption")
definition.add_optional("dust_grid", str, "the type of dust grid to use (bintree, octtree or cartesian)", "bintree")

# -----------------------------------------------------------------
