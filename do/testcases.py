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

# import standard modules
import os.path
import sys

# import the relevant PTS class
from pts.skirttestsuite import SkirtTestSuite

# get the command-line argument specifying the test suite subset, if any
subsuite = sys.argv[1] if len(sys.argv) > 1 else ""

# perform the test (sub)suite
suite = SkirtTestSuite(os.path.join("~/SKIRTtests", subsuite))
suite.performtests(skirtpath="/Volumes/DataDisk/SKIRT/release/SKIRTmain", reportpath="~/SKIRTtests", sleepsecs=10)

# -----------------------------------------------------------------
