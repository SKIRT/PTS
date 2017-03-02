#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.make_stellar_component Make a SKIRT stellar component.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema, show_parameters, merge_parameters

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

# Get the configuration parameters interactively
parameters, children = smile.prompt_parameters_for_type("PanStellarComp")

# Show the parameters
#show_parameters(parameters, children)

merged = merge_parameters(parameters, children)
print(merged)

# -----------------------------------------------------------------
