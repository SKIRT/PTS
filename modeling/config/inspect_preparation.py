#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

definition = definition.copy()

# Levels
definition.add_optional("levels", "real_list", "levels", [1.0, 3.0, 5.0, 10., 20.])

# -----------------------------------------------------------------
