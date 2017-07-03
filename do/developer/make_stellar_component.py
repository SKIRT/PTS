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
from pts.core.prep.smile import SKIRTSmileSchema

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

# Get the configuration parameters interactively
parameters, children = smile.prompt_parameters_for_type("PanStellarComp", merge=True)

# Show
print(parameters)

# -----------------------------------------------------------------
