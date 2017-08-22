#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.interpolate Interpolate an image within regions defined by the user.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.tools import interpolation
from pts.magic.region.list import load_as_pixel_region_list
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.magic.core.cutout import interpolation_methods
from pts.magic.core.detection import Detection
from pts.magic.tools import plotting

# -----------------------------------------------------------------

default_method = "biharmonic"
default_interpolation_method = "pts"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_required("image", "string", "name of the input image")
definition.add_required("regions", "string", "name of the regions file over which to interpolate")
definition.add_flag("write_mask", "write out the mask")
definition.add_optional("color", "string", "only interpolate over the shapes with this color")
definition.add_optional("ignore_color", "string", "ignore shapes with this particular color")
definition.add_optional("shapes", "string_list", "only interpolate over these kinds of shapes")
definition.add_optional("input", "string", "name of the input directory", letter="i")
definition.add_optional("output", "string", "name of the output directory", letter="o")
definition.add_optional("method", "string", "interpolation method to use", default=default_method)
definition.add_optional("interpolation_method", "string", "interpolation method", default_interpolation_method, choices=interpolation_methods)
definition.add_flag("sigma_clip", "apply sigma clipping before interpolation", True)
definition.add_optional("source_outer_factor", "real", "outer factor", 1.4)
definition.add_flag("plot", "plot after interpolation", False)
definition.add_flag("replace", "allow the original image to be replaced", False)
definition.add_flag("backup", "backup if replaced", True)
definition.add_optional("backup_suffix", "string", "backup suffix", "backup")
config = parse_arguments("interpolate", definition)

# -----------------------------------------------------------------

# Set input path
if config.input is not None: input_path = fs.absolute_or_in_cwd(config.input)
else: input_path = fs.cwd()

# Set output path
if config.output is not None: output_path = fs.absolute_or_in_cwd(config.output)
else: output_path = fs.cwd()

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the image ...")

# Determine the full path to the image
image_path = fs.absolute_path(config.image)

# Import the image
importer = ImageImporter()
importer.run(image_path)

# Get the primary image frame
frame = importer.image.primary

# Get the original header
header = importer.image.original_header

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the region ...")

# Load the region
region_path = fs.join(input_path, config.regions)
regions = load_as_pixel_region_list(region_path, frame.wcs, only=config.shapes, color=config.color, ignore_color=config.ignore_color)

# -----------------------------------------------------------------

# Inform the user
log.info("Creating a mask from the region ...")

# Create a mask from the region list
mask = regions.to_mask(frame.xsize, frame.ysize)

# -----------------------------------------------------------------

# Inform the user
log.info("Interpolating the frame within the masked pixels ...")

# Create a mask of the pixels that are NaNs
nans = frame.nans()

# Set the NaN pixels to zero in the frame
frame[nans] = 0.0

# Interpolate the frame in the masked pixels
if config.method == "biharmonic": frame._data = interpolation.inpaint_biharmonic(frame, mask)
elif config.method == "local_mean": frame._data = interpolation.in_paint(frame, mask, method="localmean")
elif config.method == "pts":

    sources = []

    # Create sources
    for shape in regions:

        # Create a source
        source = Detection.from_shape(frame, shape, config.source_outer_factor)
        sources.append(source)

    # Loop over the sources
    for source in sources:

        # Estimate the background
        source.estimate_background(config.interpolation_method, sigma_clip=config.sigma_clip)

        # Replace the pixels by the background
        source.background.replace(frame, where=source.mask)

# Invalid
else: raise ValueError("Invalid interpolation method")

# -----------------------------------------------------------------

# Set the original NaN pixels back to NaN
frame[nans] = float("nan")

# -----------------------------------------------------------------

# Plot
if config.plot: plotting.plot_box(frame)

# -----------------------------------------------------------------

# Inform the user
log.info("Saving the result ...")

# Determine the path
path = fs.join(output_path, config.image)
if fs.is_file(path):
    if config.replace:
        if config.backup:
            fs.backup_file(path, suffix=config.backup_suffix)
            fs.remove_file(path)
        else: fs.remove_file(path)
    else: raise ValueError("The image already exists")

# Save
frame.saveto(path, header=header)

# Write the mask
if config.write_mask:

    path = fs.join(output_path, "mask.fits")
    mask.saveto(path, header=header)

# -----------------------------------------------------------------
