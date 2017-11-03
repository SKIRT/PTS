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

# Processing of H-alpha image
definition.add_flag("interpolate_halpha", "interpolate the negative and NaN values in the H-alpha image", True)
definition.add_flag("smooth_halpha", "smooth the H-alpha image", True)
definition.add_optional("halpha_smoothing_factor", "positive_real", "smoothing factor", 2.)

# -----------------------------------------------------------------
