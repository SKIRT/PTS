#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.invert Invert the colors of an image.

# -----------------------------------------------------------------

# Import standard modules
import imageio

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import invert_colors, invert_black_and_white, invert_greys
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_positional_optional("filenames", "filepath_list", "names of the input image file")
definition.add_flag("bw", "only interchange black and white")
definition.add_flag("greys", "only invert shades of grey (including black and white)")

# Parse the command line arguments
config = parse_arguments("invert", definition)

# -----------------------------------------------------------------

if config.filenames is None: filenames = fs.files_in_cwd(extension="png")
else: filenames = config.filenames

# -----------------------------------------------------------------

# Loop over the image files
for filename in filenames:

    # Open the original image
    image = imageio.imread(filename)

    # Invert the colours
    if config.bw: invert_black_and_white(image)
    elif config.greys: invert_greys(image)
    else: invert_colors(image)

    # Determine output name
    name = fs.strip_extension(fs.name(filename))
    extension = fs.get_extension(filename)
    newname = name + "_inverted." + extension

    # Write the inverted image
    imageio.imwrite(newname, image)

# -----------------------------------------------------------------
