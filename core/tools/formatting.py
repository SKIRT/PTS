#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.formatting Formatting text in the terminal.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

# SOURCE: http://misc.flogisoft.com/bash/tip_colors_and_formatting

# -----------------------------------------------------------------

# Set
bold = "\033[1m"
dim = "\033[2m"
underlined = "\033[4m"
blink = "\033[5m"
inverted = "\033[7m"
hidden = "\033[8m"

# -----------------------------------------------------------------

# Reset
reset = "\033[0m"
reset_bold = "\033[21m"
reset_dim = "\033[22m"
reset_underlined = "\033[24m"
reset_blink = "\033[25m"
reset_inverted = "\033[27m"
reset_hidden = "\033[28m"

# -----------------------------------------------------------------

# Text colours
default_text = "\033[39m"
black = "\033[30m"
red = "\033[31m"
green = "\033[32m"
yellow = "\033[33m"
blue = "\033[34m"
magenta = "\033[35m"
cyan = "\033[36m"
lightgray = "\033[37m"
darkgray = "\033[90m"
lightred = "\033[91m"
lightgreen = "\033[92m"
lightyellow = "\033[93m"
lightblue = "\033[94m"
lightmagenta = "\033[95m"
lightcyan = "\033[96m"
white = "\033[97m"

# -----------------------------------------------------------------

# Background colours
default_background = "\033[49m"
black_background = "\033[40m"
red_background = "\033[41m"
green_background = "\033[42m"
yellow_background = "\033[43m"
blue_background = "\033[44m"
magenta_background = "\033[45m"
cyan_background = "\033[46m"
lightgray_background = "\033[47m"
darkgray_background = "\033[100m"
lightred_background = "\033[101m"
lightgreen_background = "\033[102m"
lightyellow_background = "\033[103m"
lightblue_background = "\033[104m"
lightmagenta_background = "\033[105m"
lightcyan_background = "\033[106m"
white_background = "\033[107m"

# -----------------------------------------------------------------
