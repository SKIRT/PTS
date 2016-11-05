#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.configuration_types List the different labels that can be used to specify a type in the context of PTS configuration definitions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from inspect import getmembers, isfunction, getdoc

# Import the relevant PTS classes and modules
from pts.core.tools import parsing
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Get a list of all the functions in the parsing module
function_list = [o for o in getmembers(parsing) if isfunction(o[1])]

names = []
descriptions = []
longest_name = 0

# Print all the function names
for function in function_list:

    name = function[0]
    documentation = getdoc(function[1])
    description = documentation.split(":param")[0].split(">>>")[0].rstrip("\n")

    description = description.replace("\n", " ")

    if len(name) > longest_name: longest_name = len(name)

    names.append(name)
    descriptions.append(description)

for i in range(len(names)):

    name = names[i]
    description = descriptions[i]
    nspaces = longest_name - len(name) + 2
    spaces = " " * nspaces

    print(fmt.yellow + name + fmt.reset + spaces + description)

# -----------------------------------------------------------------
