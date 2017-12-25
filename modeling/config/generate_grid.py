#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

default_scale = "logarithmic"
scales = ["linear", "logarithmic"]

# -----------------------------------------------------------------

definition = definition.copy()

# The number of models
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 80)

# Scales
definition.add_optional("scales", "string_string_dictionary", "scales (linear/logarithmic) for the different free parameters")

# Make animations?
definition.add_flag("animate", "make animations", False)

# -----------------------------------------------------------------
