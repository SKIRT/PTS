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

# The factor that defines the truncation boundary
definition.add_required("factor", "real", "best estimate for the value of the truncation boundary factor", suggestions=[0.82])

# -----------------------------------------------------------------
