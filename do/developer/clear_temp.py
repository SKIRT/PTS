#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.clear_temp Clear the complete PTS temporary directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Remove
fs.remove_directory(introspection.pts_temp_dir)

# Create clean
fs.create_directory(introspection.pts_temp_dir)

# -----------------------------------------------------------------
