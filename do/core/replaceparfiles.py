#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.replaceparfiles Copy \c parameters.xml and \c .tex files from \c out to \c ref test suite directories.
#
# This script copies the \c _parameters.xml and \c _parameters.tex files from the \c out to the \c ref test case
# directories for all or part of the standard SKIRT test suite in directory \c ~/SKIRTtests.
#
# When run without command line arguments, the script handles all test cases in the suite.
# Provide the name of a subsuite as the first argument to restrict handling to that subsuite.
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import sys

# -----------------------------------------------------------------

# get the command-line argument specifying the test suite subset, if any
subsuite = sys.argv[1] if len(sys.argv) > 1 else ""
subsuite = os.path.realpath(os.path.expanduser(os.path.join("~/SKIRTtests",subsuite)))

# iterate over all "out" directories that reside next to a ski file
for dirpath, dirs, files in os.walk(subsuite):
    if dirpath.endswith("/out"):
        for skifile in filter(lambda fn: fn.endswith(".ski"), os.listdir(dirpath[:-4])):
            prefix = skifile[:-4] + "_"
            # iterate over all output files with parameters.xml or parameters.tex extension
            print prefix
            for name in filter(lambda fn: fn.startswith(prefix) and   \
                              (fn.endswith("_parameters.xml") or fn.endswith("_parameters.tex")), files):
                fromfile = os.path.join(dirpath, name)
                tofile = os.path.join(os.path.join(dirpath[:-4],"ref"), name)
                print fromfile
                print tofile
                os.rename(fromfile, tofile)

# -----------------------------------------------------------------
