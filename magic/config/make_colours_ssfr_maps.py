#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.maps.ssfr.colours import ssfr_colours

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("colours", "string_list", "colours to use", default=ssfr_colours, choices=ssfr_colours)

# -----------------------------------------------------------------
