#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.plotting.photometry import PhotometryPlotter
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

definition = definition.copy()

# Add settings
definition.add_positional_optional("features", "string_list", "features to be plotted", choices=PhotometryPlotter.features())

# -----------------------------------------------------------------
