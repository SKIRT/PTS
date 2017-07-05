#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.remove_keywords Remove keywords from the header of a FITS file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronmical modules
from astropy.io import fits

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("path", "file_path", "path of the image")
definition.add_required("keys", "string_list", "keys to be removed")

config = parse_arguments("remove_keywords")

# -----------------------------------------------------------------

# Open
hdulist = fits.open(config.path)

# Get header
hdu = hdulist[0]
header = hdu.header

# Remove keys
for name in config.keys:
    del header[name]

# Overwrite
fs.remove_file(config.path)
hdu.writeto(config.path)

# -----------------------------------------------------------------
