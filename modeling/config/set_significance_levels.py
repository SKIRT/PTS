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

# The levels
definition.add_positional_optional("levels", "string_real_dictionary", "significance level for each image")

# -----------------------------------------------------------------
