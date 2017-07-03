#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.modeling_commands Show the modelling commands.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.core.steps import single_commands, repeated_commands
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

print("")
print(fmt.green + fmt.underlined + "For single invokation:" + fmt.reset)
print("")

for command in single_commands: print(" - " + command)

print("")
print(fmt.yellow + fmt.underlined + "For repeated invokation:" + fmt.reset)
print("")

for command in repeated_commands: print(" - " + command)

print("")

# -----------------------------------------------------------------
