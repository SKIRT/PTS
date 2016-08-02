#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from .fetch import definition

# -----------------------------------------------------------------

definition.add_optional("catalogs", "string_list", "names of the catalogs to query", ["GALEX", "2MASS", "SINGS", "LVL", "Spitzer", "Spitzer/IRS", "IRAS", "IRAS-FSC", "S4G", "Emission lines", "Brown", "Planck"])
#["GALEX", "2MASS", "SINGS", "LVL", "Spitzer", "Spitzer/IRS", "IRAS", "IRAS-FSC", "S4G", "Brown", "Planck"]

# -----------------------------------------------------------------
