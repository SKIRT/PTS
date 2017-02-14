#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.test_coverage Determine the test coverage.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import imp

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.basics.configuration import find_command
from pts.core.tools import formatting as fmt
from pts.core.test.pts import tests_for_subproject, path_for_test
from pts.core.tools import filesystem as fs
from pts.core.tools import sequences

# -----------------------------------------------------------------

# Get arguments tables
tables = introspection.get_arguments_tables()

# -----------------------------------------------------------------

# Get all configurable classes
classes = list(introspection.all_concrete_configurable_classes())

# Initialize container
in_tests = [False] * len(classes)

# -----------------------------------------------------------------

