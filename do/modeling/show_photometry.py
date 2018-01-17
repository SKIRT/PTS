#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_photometry Show photometry images/masks.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.magic.core.image import Image
from pts.core.tools import sequences
from pts.magic.plot.imagegrid import StandardImageGridPlotter

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("filters", "lazy_broad_band_filter_list", "filters for which to show the photometry")
definition.add_optional("fitting", "string", "show only for the fitting filters of this fitting run")
definition.add_optional("not_filters", "lazy_broad_band_filter_list", "don't show these filters")
definition.add_flag("downsample", "perform downsampling")
definition.add_optional("max_npixels", "positive_integer", "maximum number of pixels to enabled downsampling", 400)
definition.add_flag("truncate", "truncate the images", True)

# Extra
definition.add_flag("normalize", "normalize the images")
definition.add_optional("share_scale_with", "string", "share the scale of all other images with this image")
definition.add_optional("colormap", "string", "colormap", "viridis")

# Extra
definition.add_flag("write_data", "write data")

# Output
definition.add_optional("output", "directory_path", "output directory")

# Create the configuration
config = parse_arguments("show_photometry", definition, "Show photometry images/masks")

# -----------------------------------------------------------------

# Get modeling path
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Set the filters
if config.fitting is not None:

    # Load the fitting run
    fitting_run = runs.load(config.fitting)

    # Get the fitting filters
    filters = fitting_run.fitting_filters

# Filters are specified
elif config.filters is not None: filters = config.filters

# All photometry filters
else: filters = environment.photometry_filters

# Remove filters?
if config.not_filters is not None: filters = sequences.removed(filters, config.not_filters)

# -----------------------------------------------------------------

# Create plotter
plotter = StandardImageGridPlotter()

# Downsampling
plotter.config.downsample = config.downsample
plotter.config.max_npixels = config.max_npixels
plotter.config.output = config.output

# Extra
plotter.config.normalize = config.normalize
if config.share_scale_with is not None:
    plotter.config.share_scale = True
    plotter.config.scale_reference = config.share_scale_with
plotter.config.colormap = config.colormap

# Write data
plotter.config.write = config.write_data

# Crop to
plotter.crop_to=environment.truncation_box

# -----------------------------------------------------------------

# Loop over the filters
for fltr in filters:

    # Get filter name
    filter_name = str(fltr)

    # Get path
    path = environment.photometry_image_paths_for_filters[fltr]

    # Add to plot
    plotter.add_image_from_file(path, masks=False, regions=False)

# -----------------------------------------------------------------

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
