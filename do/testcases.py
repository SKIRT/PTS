#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.testcases Perform all or part of the standard SKIRT test suite.
#
# This script performs all or part of the standard SKIRT test suite in directory \c ~/SKIRTtests using the
# SKIRT executable in directory \c ~/SKIRT/release.
#
# When invoked without command line arguments, the script performs all test cases in the suite.
# Provide the name of a subsuite as the first argument to restrict execution to that subsuite.
#

# -----------------------------------------------------------------

# Import standard modules
import os.path
import sys

# Import the relevant PTS class
from pts.skirttestsuite import SkirtTestSuite
from pts.log import Log

# Create the full path to the SKIRTtests directory
suitename = "SKIRTtests"
suitepath = os.path.join(os.getenv("HOME"), suitename)

# Get the command-line argument specifying the test suite subset, if any
subsuitename = sys.argv[1] if len(sys.argv) > 1 else ""

# Create the test suite
suite = SkirtTestSuite(suitepath=suitepath, subsuitename=subsuitename)

# Perform the test suite
suite.perform(sleepsecs=10)

# -----------------------------------------------------------------
