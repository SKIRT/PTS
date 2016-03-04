#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.depends List the dependencies for a certain PTS do script.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.tools import inspection
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("script", type=str, nargs="?", help="the name of the PTS do script for which to determine the dependencies")
parser.add_argument("-m", "--modules", action="store_true", help="show the PTS modules which import a given package")
parser.add_argument("-s", "--standard", action="store_true", help="show import packages from the python standard library")
parser.add_argument("-v", "--version", action="store_true", help="show the version numbers of the required packages")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If no script name is given, execute the "list_dependencies.py" script to list all dependencies of PTS and the
# PTS modules that use them
if arguments.script is None: dependencies = inspection.get_all_dependencies()

else:

    # Find matching scripts under the 'do' directory
    match = inspection.find_matching_script(arguments.script)
    if match is None: exit()

    # Determine the full path to the matching script
    script_path = os.path.join(inspection.pts_do_dir, match[0], match[1])

    # List the dependencies of the matching script
    dependencies = defaultdict(set)
    inspection.add_dependencies(dependencies, script_path, set())

# Get the versions of all installed python packages
if arguments.version: versions = inspection.get_pip_versions()
else: versions = None

# Loop over the packages and report their presence
for dependency in sorted(dependencies, key=str.lower):

    # Get the list of PTS scripts for this dependency
    script_list = dependencies[dependency]

    # Skip packages from the standard library, unless the appropriate flag is enabled
    if inspection.is_std_lib(dependency) and not arguments.standard: continue

    # Check whether the current package is present
    if inspection.is_present(dependency):

        # Check version number
        if versions is not None and (dependency.lower() in versions): version = versions[dependency.lower()]
        else: version = None

        # Show package name, whether it's present and version number (if requested)
        if version is not None: log.success(dependency + ": present (version " + version + ")")
        else: log.success(dependency + ": present")

    # The package is not present
    else: log.error(dependency + ": not found")

    # List the PTS modules that have this dependency
    if arguments.modules:
        for script in script_list: log.info("    " + script.split("PTS/pts/")[1])

# -----------------------------------------------------------------
