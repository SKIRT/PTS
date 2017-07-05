#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.fits_to_png Convert a FITS file to a PNG image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "galaxy name")

config = parse_arguments("list_ukidss", definition)

# -----------------------------------------------------------------



# -----------------------------------------------------------------
