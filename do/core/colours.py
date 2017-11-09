#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from __future__ import print_function

# -----------------------------------------------------------------

colored = [0] + [0x5f + 40 * n for n in range(0, 5)]
colored_palette = [
	"%02x/%02x/%02x" % (r, g, b) 
	for r in colored
	for g in colored
	for b in colored
]

grayscale = [0x08 + 10 * n for n in range(0, 24)]
grayscale_palette = [
	"%02x/%02x/%02x" % (a, a, a)
	for a in grayscale 
]

normal = "\033[38;5;%sm" 
bold = "\033[1;38;5;%sm"
reset = "\033[0m"

#print(colored_palette)

for (i, color) in enumerate(colored_palette + grayscale_palette, 16):

	hex_name = "#" + color.replace("/","")

	index = (bold + "%4s" + reset) % (i, str(i) + ':')
	hex   = (normal + "%s" + reset) % (i, hex_name)

	print(index, hex)

# -----------------------------------------------------------------
