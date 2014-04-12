#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.image Experiment with images.
#

# -----------------------------------------------------------------

# import standard modules
import sys
import os.path

# import relevant pts modules
from pts.rgbimage import RGBImage

# -----------------------------------------------------------------

print "Starting..."

im = RGBImage("refblue.png")
im.addbelow(RGBImage("refred.png"))
im.enlargecanvas((130,70))
im.saveto("resultredbelow.png")

im = RGBImage("refblue.png")
im.addright(RGBImage("refred.png"))
im.enlargecanvas((120,90))
im.saveto("resultredright.png")

print "Finished."

# -----------------------------------------------------------------
