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

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import invert_colors
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "name of the input image file")

# Parse the command line arguments
setter = ArgumentConfigurationSetter("invert")
config = setter.run(definition)

# -----------------------------------------------------------------

# Open the original image
image = imageio.imread(config.filename)

# Invert the colours
invert_colors(image)

# Determine output name
name = fs.strip_extension(fs.name(config.filename))
extension = fs.get_extension(config.filename)
newname = name + "_inverted." + extension

# Write the inverted image
imageio.imwrite(newname, image)

# -----------------------------------------------------------------
