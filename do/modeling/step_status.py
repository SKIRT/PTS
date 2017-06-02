#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.step_status View the status of the modeling steps.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.modeling.component.component import load_modeling_history
from pts.modeling.core.steps import single_commands, repeated_commands
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Get the history
history = load_modeling_history(fs.cwd())

print("")
name = fs.name(modeling_path)
print(fmt.underlined + name + fmt.reset)
print("")

nfinished = 0
ntotal = 0

for command in single_commands:

    ntotal += 1

    if history.finished(command):
        print(fmt.green + " - " + command + ": finished" + fmt.reset)
        nfinished += 1
    elif history.finished(command): print(fmt.yellow + " - " + command + ": started" + fmt.reset)
    else: print(fmt.red + " - " + command + ": not started" + fmt.reset)

print("")

fraction = float(nfinished) / float(ntotal)

print(tostr(fraction*100, decimal_places=0) + "% completed")
print("")

# -----------------------------------------------------------------
