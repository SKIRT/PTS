#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.count Count the number of lines of code in the PTS project, the number of functions
#  and the number of modules.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import importlib
from collections import defaultdict
from inspect import getmembers, isfunction, getdoc

# Import the relevant PTS classes and modules
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import setup_log

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_flag("verbose", "verbose output", letter="v")
config = parse_arguments("count", definition)

# -----------------------------------------------------------------

# Setup log
setup_log("ERROR")

# -----------------------------------------------------------------

# Lines
nlines = 0
nlines_per_file = dict()
nlines_per_subproject = defaultdict(int)

# Modules
nmodules = 0
nmodules_per_subproject = defaultdict(int)

# Functions
nfunctions = 0
nfunctions_per_subproject = defaultdict(int)

# Loop over all module files
for path, subproject, module in introspection.all_modules():

    nmodules += 1
    nmodules_per_subproject[subproject] += 1

    # Loop over the members
    for member in getmembers(module):

        if isfunction(member[1]):

            nfunctions += 1
            nfunctions_per_subproject[subproject] += 1

    # Open the module file
    with open(path, 'r') as pyfile:

        nlines_per_file[path] = 0

        for line in pyfile:

            line = line.rstrip("\n")

            if line == "": continue
            if line.startswith("#"): continue

            nlines += 1

            nlines_per_file[path] += 1
            nlines_per_subproject[subproject] += 1

# -----------------------------------------------------------------

print("")

# State the number of lines
print(fmt.green + fmt.bold + str(nlines) + " lines of python code" + fmt.reset)

# Verbose output
if config.verbose:

    # State the number of lines for each module
    #for path in nlines_per_file:
    #    print(" - " + path + ": " + str(nlines_per_file[path]) + " lines")
    print("")
    for subproject in nlines_per_subproject:
        print(" - " + subproject + ": " + str(nlines_per_subproject[subproject]) + " lines")

if config.verbose: print("")

# Modules
print(fmt.green + fmt.bold + str(nmodules) + " module files" + fmt.reset)

if config.verbose:

    # State number of modules per subproject
    print("")
    for subproject in nmodules_per_subproject:
        print(" - " + subproject + ": " + str(nmodules_per_subproject[subproject]) + " modules")

if config.verbose: print("")

# Functions
print(fmt.green + fmt.bold + str(nfunctions) + " functions" + fmt.reset)

if config.verbose:

    # State number of functions per subproject
    print("")
    for subproject in nfunctions_per_subproject:
        print(" - " + subproject + ": " + str(nfunctions_per_subproject[subproject]) + " functions")

print("")

# -----------------------------------------------------------------
