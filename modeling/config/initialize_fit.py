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
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add optional arguments
definition.add_section("wavelengths", "settings for the wavelength grid")
definition.sections["wavelengths"].add_optional("unit", "string", "the unit of the wavelengths", "micron")
definition.sections["wavelengths"].add_optional("min", "real", "the minimum wavelength", 0.09)
definition.sections["wavelengths"].add_optional("max", "real", "the maximum wavelength", 2000)
definition.sections["wavelengths"].add_optional("npoints", "integer", "the number of wavelength points", 100)
definition.sections["wavelengths"].add_optional("min_zoom", "real", "the minimum wavelength of the zoomed-in grid", 1)
definition.sections["wavelengths"].add_optional("max_zoom", "real", "the maximum wavelength of the zoomed-in grid", 30)
definition.sections["wavelengths"].add_optional("npoints_zoom", "integer", "the number of wavelength points in the zoomed-in grid", 100)

definition.add_optional("packages", "real", "the number of photon packages per wavelength", 2e5)
definition.add_flag("selfabsorption", "enable dust self-absorption")
definition.add_optional("dust_grid", "string", "the type of dust grid to use (bintree, octtree or cartesian)", "bintree")

# -----------------------------------------------------------------
