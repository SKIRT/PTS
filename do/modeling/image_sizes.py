#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.image_sizes

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.modeling.component.galaxy import get_data_images_path
from pts.magic.core.frame import Frame
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

data_images_path = get_data_images_path(modeling_path)

# -----------------------------------------------------------------

shapes = dict()

# Loop over the images
for image_path, image_name in fs.files_in_path(data_images_path, extension="fits", not_contains="poisson", returns=["path", "name"], recursive=True, recursion_level=1):

    # Load the image
    frame = Frame.from_file(image_path, silent=True)

    # Determine filter name
    filter_name = frame.filter_name

    # Set pixel shape
    shapes[filter_name] = frame.shape

# -----------------------------------------------------------------

# Sort on shape
sorted_filter_names = sorted(shapes.keys(), key=lambda name: shapes[name])

# -----------------------------------------------------------------

print("")
fmt.print_columns(sorted_filter_names, [shapes[name] for name in sorted_filter_names])
print("")

# -----------------------------------------------------------------
