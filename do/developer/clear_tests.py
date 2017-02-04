#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.clear_tests Clear all test output.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Determine path
path = fs.join(introspection.pts_temp_dir, "tests")

# Remove
fs.remove_directory(path)

# Create clean
fs.create_directory(path)

# -----------------------------------------------------------------
