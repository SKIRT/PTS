#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.launch.pts import execute_pts_local

# -----------------------------------------------------------------

# Interpolate
execute_pts_local("interpolate", "result.fits", "post.reg", replace=True, backup_suffix="before_post", method="pts", debug=True, show_output=True)

# Add all original planes
execute_pts_local("add_planes", "result.fits", "result_before_post.fits", replace_masks=True, debug=True, show_output=True)

# -----------------------------------------------------------------
