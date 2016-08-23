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

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints_range", "integer_range", "the range of the wavelength grid size", "150>500", convert_default=True)
definition.sections["wg"].add_optional("ngrids", "integer", "the number of wavelength grids to generate", 10)

# Settings for the dust grid generation
definition.add_section("dg", "settings for the dust grids")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("scale_range", "real_range", "the range of the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", "0.5>10.", convert_default=True)
definition.sections["dg"].add_optional("level_range", "integer_range", "the range of the maximum depth level of the tree", "6>9", convert_default=True)
definition.sections["dg"].add_optional("mass_fraction_range", "real_range", "the range of the maximum mass fraction in each cell", "0.5e-6>1e-5", convert_default=True)

# Add optional arguments
definition.add_optional("npackages", "real", "the number of photon packages per wavelength", 2e5)
definition.add_flag("selfabsorption", "enable dust self-absorption")

# -----------------------------------------------------------------
