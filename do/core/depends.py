#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.depends List the dependencies for a certain PTS do script
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import inspection
from pts.core.basics import Log

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("script", type=str, help="the name of the PTS do script for which to determine the dependencies")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Find matching scripts under the 'do' directory
match = inspection.find_matching_script(arguments.script)
if match is None: exit()

# Determine the full path to the matching script
script_path = os.path.join(inspection.pts_do_dir, match[0], match[1])

# List the dependencies of the matching script
dependencies = set()
inspection.add_dependencies(dependencies, script_path)

# Loop over the packages and report their presence
for dependency in sorted(dependencies, key=str.lower):

    if inspection.is_present(dependency): log.success(dependency + ": present")
    else: log.failure(dependency + ": not found")

# -----------------------------------------------------------------
