#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.plotting.analysis import get_features
from pts.modeling.config.plot import definition
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Add settings
features = get_features(fs.cwd())
definition.pos_optional["features"].choices = features
definition.pos_optional["features"].default = features.keys()

# -----------------------------------------------------------------
