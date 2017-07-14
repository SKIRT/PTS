#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.image_pixelscales

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.modeling.component.galaxy import get_data_images_path
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.tools import headers
from pts.magic.core.fits import get_header
from pts.core.tools import formatting as fmt
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

data_images_path = get_data_images_path(modeling_path)

# -----------------------------------------------------------------

pixelscales = dict()

# Loop over the images
for image_path, image_name in fs.files_in_path(data_images_path, extension="fits", not_contains="poisson", returns=["path", "name"], recursive=True, recursion_level=1):

    # Load the images
    wcs = CoordinateSystem.from_file(image_path)

    # Get pixelscale
    pixelscale = wcs.average_pixelscale

    # Get filter name
    header = get_header(image_path)
    fltr = headers.get_filter(image_name, header)
    filter_name = str(fltr)

    # Set pixelscale
    pixelscales[filter_name] = pixelscale

# -----------------------------------------------------------------

# Sort on pixelscale
sorted_filter_names = sorted(pixelscales.keys(), key=lambda name: pixelscales[name])

# -----------------------------------------------------------------

print("")
fmt.print_columns(sorted_filter_names, [pixelscales[name] for name in sorted_filter_names])
print("")

# -----------------------------------------------------------------
