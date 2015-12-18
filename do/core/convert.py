#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.convert Perform data conversion for some specific SPH data file(s).
#
# This script performs data conversion (from foreign to SKIRT format) for some specific SPH data file(s).
# The in/out filenames and other parameters are hardcoded in the script.
# Thus the script mainly serves as an example of how to use the pts.sphconvert module.
#

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.prep.sphconvert import *

# -----------------------------------------------------------------

#convert_stars_EAGLE("galaxy_stars_035_007_147.data", "eagle_stars.dat")
#convert_gas_EAGLE("galaxy_gas_035_007_147.data", "eagle_gas.dat")

#convert_stars_AWAT("s000284.data", "awat_stars.dat")
#convert_gas_AWAT("g000284.data", "awat_gas.dat")

#convert_stars_DOLAG("stars_3r.52.215.dat", "dolag_stars.dat")

#convert_gas_DOLAG("gas_3r.52.215.dat", "dolag_gas.dat")

convert_gas_ULB("yepgas.ascii", "yepgas.dat")

# -----------------------------------------------------------------
