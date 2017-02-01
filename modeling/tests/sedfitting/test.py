#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.pts Contains the PTSTestSuite class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

commands = []
settings = []

# -----------------------------------------------------------------

settings_setup = dict()
settings_setup["type"] = "other"
settings_setup["name"] = "SN1987"
settings_setup["fitting_host_ids"] = None

# -----------------------------------------------------------------

commands.append("setup")
commands.append("model_sed")

# -----------------------------------------------------------------



# -----------------------------------------------------------------

def test():


    return

# -----------------------------------------------------------------