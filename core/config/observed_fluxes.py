#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_optional("output", "string", "output directory")

# Add flags
definition.add_flag("spectral_convolution", "convolve over the wavelengths to get the most accurate fluxes", True)

# -----------------------------------------------------------------

definition.add_flag("from_images", "calculate observed fluxes from images created from the output datacubes", False)
definition.add_flag("write_images", "write the images created from the output datacubes", False)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plot the mock observed SEDs", False)

# -----------------------------------------------------------------

definition.add_optional("images_nprocesses", "positive_integer", "number of parallel processes for creating images", 1)

# -----------------------------------------------------------------
