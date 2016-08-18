#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Wavelength grid settings
definition.add_section("wavelengths", "settings for the wavelength grid")
definition.sections["wavelengths"].add_optional("unit", "string", "the unit of the wavelengths", "micron")
definition.sections["wavelengths"].add_optional("min", "real", "the minimum wavelength", 0.05)
definition.sections["wavelengths"].add_optional("max", "real", "the maximum wavelength", 1000)
definition.sections["wavelengths"].add_optional("npoints", "integer", "the number of wavelength points", 300)
definition.sections["wavelengths"].add_optional("min_zoom", "real", "the minium wavelength of the zoomed-in grid", 1)
definition.sections["wavelengths"].add_optional("max_zoom", "real", "the maximum wavelength of the zoomed-in grid", 30)
definition.sections["wavelengths"].add_optional("npoints_zoom", "integer", "the number of wavelength points in the zoomed-in grid", 300)

# Add optional arguments
definition.add_optional("packages", "real", "the number of photon packages per wavelength", 1e6)
definition.add_optional("selfabsorption", "boolean", "whether self-absorption should be enabled", True)
definition.add_optional("remote", "string", "the remote host on which to launch the simulations", "nancy", choices=find_host_ids())
definition.add_optional("images_remote", "string", "the remote host on which to make the observed images", "nancy", choices=find_host_ids())

# -----------------------------------------------------------------
