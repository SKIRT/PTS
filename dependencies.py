#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dependencies Print each PTS dependency

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from core.tools import introspection

# -----------------------------------------------------------------

# Print each dependency on a seperate line
for package in introspection.get_all_dependencies(): print(package)

# -----------------------------------------------------------------
