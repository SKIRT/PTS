#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.magic.config.view import definition

# -----------------------------------------------------------------

# Cretae configuration definition
definition.add_required("image", "file_path", "image path")
definition.add_optional("regions", "file_path", "regions file path")
definition.add_optional("mask", "file_path", "mask file path")

# -----------------------------------------------------------------

definition.add_optional("with", "positive_integer", "image width", 500)
definition.add_optional("height", "positive_integer", "image height")

# -----------------------------------------------------------------
