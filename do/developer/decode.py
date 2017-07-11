#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.decode Decode a control input.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys

# -----------------------------------------------------------------

# FROM https://stackoverflow.com/questions/36205970/python-send-control-q-then-control-a-special-keys

# Is this really working?

while True:
    inp = sys.stdin.read(1)
    if len(inp) == 0:
        break
    print(ord(inp[0]))

# -----------------------------------------------------------------
