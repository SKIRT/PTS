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
definition.add_flag("regenerate", "regenerate pages", True)
definition.add_flag("show", "show", False)
definition.add_flag("replot", "make plots again", True)

definition.add_flag("regenerate_index", "regenerate_index", False)

# Add detail pages
definition.add_flag("details", "add detailed pages", True)

# -----------------------------------------------------------------
