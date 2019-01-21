#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.testcases Perform all or part of the standard SKIRT test suite.
#
# This script performs all or part of the standard SKIRT test suite in directory \c ~/SKIRTtests using the
# SKIRT executable in directory \c ~/SKIRT/release.
#
# When invoked without command line arguments, the script performs all test cases in the suite.
# Provide the name of a subsuite as the first argument to restrict execution to that subsuite.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os.path
import argparse

# Import the relevant PTS classes and modules
from pts.core.test.skirt import SkirtTestSuite

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('subsuite', type=str, help='a name identifying the subsuite', nargs='?', default=None)
parser.add_argument('-p', '--parallel', action='store_true', help='execute the test cases in parallel mode')

# Parse the command line arguments
arguments = parser.parse_args()
subsuite = arguments.subsuite
parallel = arguments.parallel

# -----------------------------------------------------------------

# Create the full path to the SKIRTtests directory
suitename = "SKIRT/Functional8"
suitepath = os.path.join(os.getenv("HOME"), suitename)

# Check whether a development SKIRT repository is present, otherwise use the standard SKIRT path
devskirtpath = os.path.join(os.getenv("HOME"), "SKIRT", "SKIRT8", "release", "SKIRT", "main", "skirt")
skirtpath = devskirtpath if os.path.isfile(devskirtpath) else ""

# Create the test suite instance
suite = SkirtTestSuite(suitepath=suitepath, subsuite=subsuite, parallel=parallel, skirtpath=skirtpath)

# Perform the test suite
suite.perform(sleepsecs=10)

# -----------------------------------------------------------------
