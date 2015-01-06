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
testpath = os.path.join(os.getenv("HOME"), "SKIRTtests")

# Create a list of the subsuites to be tested
if len(sys.argv) > 1: 
    # Get the command-line arguments specifying the test subsuites
    suitenames = sys.argv[1:]
else:
    # Automatically obtain a list of all the subdirectories of the SKIRTtests directory
    suitenames = [name for name in os.listdir(testpath) if os.path.isdir(os.path.join(testpath, name))]

# Create the logger
log = Log(testpath)

# Create the test suite
suite = SkirtTestSuite(suitedirpath=testpath, log=log)

# Show which subsuites will be tested, while adding them to the suite
log.info("Subsuites that will be tested:"),
for suitename in suitenames:
    log.info("  -  " + suitename)
    suite.add(suitename, recursive=True)

# Perform the test suite
suite.perform(sleepsecs=10)

# -----------------------------------------------------------------
