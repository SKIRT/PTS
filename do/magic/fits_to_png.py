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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "galaxy name")

setter = ArgumentConfigurationSetter("list_ukidss")
config = setter.run(definition)

# -----------------------------------------------------------------



# -----------------------------------------------------------------
