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

# -----------------------------------------------------------------

# Add flags
definition.add_flag("replot", "replot", True)
definition.add_flag("show", "show the page", False)

# -----------------------------------------------------------------
