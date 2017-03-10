#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.invert Invert the colors of an image.

# -----------------------------------------------------------------

# Import standard modules
import imageio
import argparse

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import invert_colors
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic options
parser.add_argument("filename", type=str, help="the name of the input image file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Open the original image
image = imageio.imread(arguments.filename)

# Invert the colours
invert_colors(image)

# Determine output name
name = fs.strip_extension(fs.name(arguments.filename))
extension = fs.get_extension(arguments.filename)
newname = name + "_inverted." + extension

# Write the inverted image
imageio.imwrite(newname, image)

# -----------------------------------------------------------------
