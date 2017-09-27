#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.animate Make an animation from a series of images in the current working directory.

# -----------------------------------------------------------------

import imageio

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.animation import Animation
from pts.core.basics.apng import APNG
from pts.core.basics.numpngw import write_apng

# -----------------------------------------------------------------

formats = ["avi", "gif", "apng"]
default_format = "gif"

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()

# File format
definition.add_optional("format", "string", "output format", default_format, choices=formats)
definition.add_optional("fps", "positive_integer", "frames per second", 10)

# Parse the command line arguments
config = parse_arguments("animate", definition)

# -----------------------------------------------------------------

# Add the frames
paths = fs.files_in_cwd(extension="png", sort=int)

# Determine path
path = "animation." + config.format

# APNG
if config.format == "apng":

    delay = 1000/config.fps  # in miliseconds

    # CRASHES
    #animation = APNG.from_files(paths, delay=delay)
    # Save the animation
    #animation.save(path)

    # SLOW :(
    seq = [imageio.imread(frame_path) for frame_path in paths]
    #print(seq)
    write_apng(path, seq, delay=delay) #default_image=im_all, use_palette=True)

# OTHER
else:

    # Initialize the animation
    animation = Animation()
    animation.fps = config.fps

    for frame_path in paths:

        # Add frame
        animation.add_frame_from_file(frame_path)

    # Write the animation
    animation.saveto(path)

# -----------------------------------------------------------------
