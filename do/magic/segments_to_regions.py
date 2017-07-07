#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.segments_to_regions Convert a segmentation map to a region list.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.tools.segments import segments_to_regions
from pts.magic.core.segmentationmap import SegmentationMap

# -----------------------------------------------------------------

# Initialize definition
definition = ConfigurationDefinition()

# Add options
definition.add_required("segments", "file_path", "segmentation map")
definition.add_required("regions", "string", "path of the output regions file")

# Add optional
definition.add_optional("mask", "string", "path to an output mask file")
definition.add_optional("offset", "real", "offset to make the regions larger or smaller, in pixels")

#definition.add_optional("", )
# parser.add_argument("--xc", nargs='?', const=1, help="Optional: If you don't want to mask the galaxy with the given x-coordinate of the center",type=str, default=None)
# parser.add_argument("--yc", nargs='?', const=1, help="Optional: If you don't want to mask the galaxy with the given y-coordinate of the center",type=str, default=None)

# parser.add_argument("--fits_slice", nargs='?', const=1, help="Optional: Fits slice in the segmentation map to be used. Default 0.",type=str, default='0')

# Get the configuration
config = parse_arguments("pix_to_sky", definition)

# -----------------------------------------------------------------

# Open the segmentation map
segments = SegmentationMap.from_file(config.segments)
output_region_file = config.regions
output_mask_image = config.mask

#fits_slice = config.fits_slice
fits_slice = None
offset = config.offset
#xc = config.xc
#yc = config.yc
xc = yc = None

# Call the segments to region function
segments_to_regions(segments, output_region_file, output_mask_image=output_mask_image, fits_slice = fits_slice, offset=offset, xc=xc, yc=yc)

# -----------------------------------------------------------------
