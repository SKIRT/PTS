#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_residuals Show residuals of photometry images and mock observed images in current directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import find_modeling_environment_up_cwd
from pts.magic.plot.imagegrid import ResidualImageGridPlotter
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_optional("not_filters", "lazy_broad_band_filter_list", "don't show these filters")

definition.add_flag("downsample", "perform downsampling")
definition.add_optional("max_npixels", "positive_integer", "maximum number of pixels to enabled downsampling", 400)
definition.add_flag("truncate", "truncate the images", True)

definition.add_flag("distributions", "show residual distributions", False)

# Extra
definition.add_flag("normalize", "normalize the images")
definition.add_optional("share_scale_with", "string", "share the scale of all other images with this image")
definition.add_optional("colormap", "string", "colormap", "viridis")

# Extra
definition.add_flag("write_data", "write data")

# Output
definition.add_optional("output", "directory_path", "output directory")

# Create the configuration
config = parse_arguments("show_residuals", definition, "Show residuals of photometry images and mock observed images in current director")

# -----------------------------------------------------------------

# Get modeling path
environment = find_modeling_environment_up_cwd()

# -----------------------------------------------------------------

# Create plotter
plotter = ResidualImageGridPlotter()

# Plot residual distributions
plotter.config.distributions = config.distributions

# Downsampling
#plotter.config.downsample = config.downsample
#plotter.config.max_npixels = config.max_npixels
plotter.config.output = config.output

# Extra
#plotter.config.normalize = config.normalize
#if config.share_scale_with is not None:
#    plotter.config.share_scale = True
#    plotter.config.scale_reference = config.share_scale_with
#plotter.config.colormap = config.colormap

#plotter.config.coordinates = True

plotter.config.max_nrows = 3
plotter.config.ngrids = 4

# Write data
plotter.config.write = config.write_data

# Crop to
plotter.crop_to = environment.truncation_box

# -----------------------------------------------------------------

# Loop over the frames in the current working directory
for path, name in fs.files_in_cwd(extension="fits", returns=["path", "name"]):

    # Load the frame
    frame = Frame.from_file(path)

    # Get the filter
    fltr = frame.filter
    filter_name = str(fltr)

    # Do we have a photometry image for this filter
    if not environment.has_photometry_for_filter(fltr): continue

    # Get the photometry image
    reference_path = environment.photometry_image_paths_for_filters[fltr]
    reference = Frame.from_file(reference_path)

    # Add row
    plotter.add_row(reference, frame, filter_name)

# -----------------------------------------------------------------

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
