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

# Add required setting
definition.add_required("step", "string", "the modeling step for which to create the report")

# -----------------------------------------------------------------
