#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.add_planes Add planes from a FITS file into another FITS file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.magic.core.image import Image

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_required("add_to", "file_path", "add to this FITS file")
definition.add_required("add_from", "file_path", "add from this FITS file")
definition.add_flag("frames", "add frames", True)
definition.add_flag("masks", "add masks", True)
definition.add_flag("segments", "add segmentation maps", True)
definition.add_flag("replace", "replace planes", False)
definition.add_flag("replace_frames", "replace frames", False)
definition.add_flag("replace_masks", "replace masks", False)
definition.add_flag("replace_segments", "replace segmentation maps", False)
definition.add_flag("backup", "make backup", False)
config = parse_arguments("interpolate", definition)

# -----------------------------------------------------------------

if config.replace: config.replace_frames = config.replace_masks = config.replace_segments = True

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the images ...")

# Load
to_image = Image.from_file(config.add_to)
from_image = Image.from_file(config.add_from)

# -----------------------------------------------------------------

# Add the frames
if config.frames:

    # Inform the user
    log.info("Adding the frames ...")

    # Loop over the frames
    for frame_name in from_image.frames:

        # Check
        if frame_name in to_image.frames:
            if config.replace_frames: to_image.remove_frame(frame_name)
            else: continue

        # Add the frame
        to_image.add_frame(from_image.frames[frame_name], frame_name)

# -----------------------------------------------------------------

# Add the masks
if config.masks:

    # Inform the user
    log.info("Adding the masks ...")

    # Loop over the masks
    for mask_name in from_image.masks:

        # Check
        if mask_name in to_image.masks:
            if config.replace_masks: to_image.remove_mask(mask_name)
            else: continue

        # Add the mask
        to_image.add_mask(from_image.masks[mask_name], mask_name)

# -----------------------------------------------------------------

# Add the segmentation maps
if config.segments:

    # Inform the user
    log.info("Adding the segmentation maps ...")

    # Add the segmentation maps
    for segments_name in from_image.segments:

        # Check
        if segments_name in to_image.segments:
            if config.replace_segments: to_image.remove_segments(segments_name)
            else: continue

        # Add the map
        to_image.add_segments(from_image.segments[segments_name], segments_name)

# -----------------------------------------------------------------

# Save
to_image.save()

# -----------------------------------------------------------------
