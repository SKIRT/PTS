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
from pts.core.tools import parsing
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools.logging import setup_log

# -----------------------------------------------------------------

# Setup log
setup_log("ERROR")

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_flag("verbose", "verbose output", letter="v")
setter = ArgumentConfigurationSetter("count")
config = setter.run(definition)

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
for path, name in fs.files_in_path(introspection.pts_package_dir, extension="py", recursive=True, returns=["path", "name"]):

    # Skip module initialization files, they are empty anyways
    if name == "__init__": continue
    if name == "__main__": continue
    if name == "run_queue": continue
    if name == "enable_qch_mathjax": continue
    if path.endswith("eagle/config.py"): continue
    if path.endswith("eagle/collections.py"): continue
    if path.endswith("eagle/database.py"): continue
    if path.endswith("eagle/galaxy.py"): continue
    if path.endswith("eagle/plotresults.py"): continue
    if path.endswith("eagle/runner.py"): continue
    if path.endswith("eagle/scheduler.py"): continue
    if path.endswith("eagle/skirtrun.py"): continue
    if name == "fit2BB_Md": continue

    # Get subproject
    subproject = path.split("pts/")[1].split("/")[0]
    nmodules_per_subproject[subproject] += 1

    nmodules += 1

    is_config = path.split(subproject + "/")[1].split("/")[0] == "config"

    #print(is_config)

    if subproject != "do" and not is_config:

        # Load the module, get number of functions
        relpath = "pts." + path.split("pts/")[1].replace("/", ".")[:-3]
        #import_name = relpath.split(".")[-1]
        #relpath = relpath.split("." + import_name)[0]
        #print(import_name)
        try:
            for member in getmembers(importlib.import_module(relpath)):

                if isfunction(member[1]):

                    nfunctions += 1
                    nfunctions_per_subproject[subproject] += 1

        except ImportError: pass

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
