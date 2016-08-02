#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.verify_fits Verify a batch of FITS files in the current working directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.io import fits

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Loop over all files in the current path
for path in fs.files_in_path(fs.cwd(), extension="fits"):

    # Open the file
    hdulist = fits.open(path)

# -----------------------------------------------------------------
