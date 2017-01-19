#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.test_unit_conversion Test the unit conversion things.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.unit import PhotometricUnit

# -----------------------------------------------------------------

# Create nanomaggy unit
nanomaggy = PhotometricUnit("nMgy")

# Conversion factor to Jansky
factor = nanomaggy.conversion_factor("Jy")

print(factor, 3.613e-6)

# -----------------------------------------------------------------
