#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.transparent Make certain colors of an image transparent.

# -----------------------------------------------------------------

# Import standard modules
import imageio
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.rgbimage import set_transparent_color, set_transparency_gradient
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

black = (0,0,0,)
white = (255,255,255,)

# -----------------------------------------------------------------

colors = dict()
colors["black"] = black
colors["white"] = white

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_positional_optional("filenames", "filepath_list", "names of the input image file")
definition.add_optional("color", "string", "color to set to transparent", "white")
definition.add_optional("tolerance", "percentage", "percentage of tolerance of deviation from exact color")
definition.add_flag("gradient", "use gradient")
definition.add_optional("scale", "percentage", "scale for gradient", "25", convert_default=True)

# Parse the command line arguments
config = parse_arguments("transparent", definition)

# -----------------------------------------------------------------

if config.filenames is None: filenames = fs.files_in_cwd(extension="png")
else: filenames = config.filenames

# -----------------------------------------------------------------

# Loop over the image files
for filename in filenames:

    # Open the original image
    image = imageio.imread(filename)

    # Standard color
    if config.color in colors: red, blue, green = colors[config.color]

    # User-specified color
    elif "," in config.color: red, blue, green = map(int, config.color.split(","))

    # Invalid
    else: raise ValueError("Invalid input")

    # Check dimensions
    if image.shape[2] == 3:
        shape = list(image.shape)
        old_image = image
        shape[2] = 4
        shape = tuple(shape)
        image = np.ones(shape) * 255
        image[:, :, :3] = old_image

    # Set transparent
    if config.gradient: set_transparency_gradient(image, (red, blue, green,), config.scale)
    else: set_transparent_color(image, (red, blue, green,), tolerance=config.tolerance)

    # Determine output name
    name = fs.strip_extension(fs.name(filename))
    #extension = fs.get_extension(filename)
    extension = "png" # because transparancy has to be supported
    newname = name + "_transparent." + extension

    # Write the inverted image
    imageio.imwrite(newname, image)

# -----------------------------------------------------------------
