#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.step_status View the status of the modeling steps.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.modeling.component.component import load_modeling_status

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

status = load_modeling_status(modeling_path)

# -----------------------------------------------------------------

status.show()

# -----------------------------------------------------------------
