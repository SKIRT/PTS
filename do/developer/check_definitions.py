#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.check_definitions For each configurable class,
#  check whether the corresponding configuration definition is present.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.basics.configuration import find_command
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Loop over all configurable classes
for cls in introspection.all_concrete_configurable_classes():

    # Find the corresponding command
    command_name, class_name, configuration_module_path, description = find_command(cls)

    # Show
    if command_name is not None: print(fmt.green + class_name + fmt.reset)
    else: print(fmt.red + class_name + fmt.reset)

# -----------------------------------------------------------------
