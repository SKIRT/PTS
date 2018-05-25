#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.plot_frame Plot a frame.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.frame import Frame, regularize_frame
from pts.magic.tools import plotting
from pts.core.basics.distribution import Distribution
from pts.core.plot.distribution import plot_distribution
from pts.core.tools import filesystem as fs
from pts.magic.region.list import load_region_list
from pts.core.basics.plot import all_colormaps

# -----------------------------------------------------------------

default_interval = "pts"
default_scale = "log"
default_symmetric_method = "mean"
default_cmap = "inferno"

# -----------------------------------------------------------------

definition = ConfigurationDefinition()

# The frame
definition.add_required("frame", "file_path", "path of the frame")

# Save to path
definition.add_optional("path", "string", "save plot file")

# The colormap and scale
definition.add_positional_optional("cmap", "string", "colormap", default_cmap, choices=all_colormaps)
definition.add_optional("scale", "string", "scale", default_scale)
definition.add_optional("min", "real", "minimum value")
definition.add_optional("max", "real", "maximum value")

# Take absolute value
definition.add_flag("absolute", "take absolute values")

# Flags
definition.add_flag("no_nans", "remove NaNs (replace by zero)")
definition.add_flag("no_infs", "remove infinities (replace by zero)")
definition.add_flag("no_negatives", "remove negative values (replace by zero)")
definition.add_flag("no_positives", "remove positive values (replace by zero)")
definition.add_optional("cutoff_above", "real", "cutoff above this value")
definition.add_optional("cutoff_below", "real", "cutoff below this value")
definition.add_flag("dilate_nans", "dilate the NaNs (erode the image)")
definition.add_optional("dilation_radius", "real", "radius of dilation disk (in number of pixels)", 5)
definition.add_optional("dilation_niterations", "positive_integer", "number of iterations of dilation", 1)

# Show distribution
definition.add_flag("distribution", "show the distribution of pixel values")
definition.add_flag("sigma_clip", "apply sigma-clipping on the distribution values", False)

# Interval
definition.add_optional("interval", "string", "interval", default_interval)
#definition.add_optional("minmax", "real_pair", "minimum and maximum value")
definition.add_flag("soft_min", "interpret minimum of minmax as soft")
definition.add_flag("soft_max", "interpret maximum of minmax as soft")

# Contours
definition.add_flag("contours", "show contours", False)
definition.add_optional("ncontours", "positive_integer", "number of contour levels", 5)
definition.add_optional("contours_color", "string", "color for the contour lines", "white")

# Other options
definition.add_flag("around_zero", "around zero", False)
definition.add_flag("check_around_zero", "check around zero", False)
definition.add_flag("symmetric", "symmetric", False)
definition.add_optional("symmetric_method", "string", "symmetric method", default_symmetric_method)
definition.add_optional("normalize_in", "file_path", "region file for normalization to calculate within")

# Add colorbar
definition.add_flag("colorbar", "add colorbar", True)

# Get the configuration
config = parse_arguments("plot_frame", definition)

# -----------------------------------------------------------------

# Get the image name
name = fs.strip_extension(fs.name(config.frame))

# -----------------------------------------------------------------

# Open the frame
frame = Frame.from_file(config.frame)

# -----------------------------------------------------------------

if config.min is not None and config.max is not None: minmax = (config.min, config.max)
else: minmax = None

# -----------------------------------------------------------------

# Soft min
if config.soft_min:
    #if config.minmax is None: raise ValueError("Cannot enable soft_min flag if minmax is not specified")
    if config.min is None: raise ValueError("Cannot enable soft_min flag if min is not specified")
    soft_zmin = True
    #zmin = config.minmax[0]
    zmin = config.min
else: zmin, soft_zmin = None, False

# Soft max
if config.soft_max:
    #if config.minmax is None: raise ValueError("Cannot enable soft_max flag if minmax is not specified")
    if config.max is None: raise ValueError("Cannot enable soft_max flag if max is not specified")
    soft_zmax = True
    #zmax = config.minmax[1]
    zmax = config.max
else: zmax, soft_zmax = None, False

# -----------------------------------------------------------------

# Set interval
if soft_zmin or soft_zmax:
    interval = config.interval
    #if zmin is None and config.minmax is not None: zmin = config.minmax[0]
    if zmin is None and config.min is not None: zmin = config.min
    #if zmax is None and config.minmax is not None: zmax = config.minmax[1]
    if zmax is None and config.max is not None: zmax = config.max

# Both min and max are specified
elif minmax is not None: interval = minmax

# Maybe min or max is specified, or none
else:

    # min specified?
    if config.min is not None: zmin = config.min
    if config.max is not None: zmax = config.max

    # Set interval
    interval = config.interval

# -----------------------------------------------------------------

# Load regions
if config.normalize_in is not None: normalize_in = load_region_list(config.normalize_in)
else: normalize_in = None

# -----------------------------------------------------------------

#print(interval)
if config.symmetric: config.around_zero = True

# Determine min and max
vmin, vmax = plotting.get_vmin_vmax(frame.data, interval=interval, around_zero=config.around_zero,
                                    symmetric=config.symmetric, normalize_in=normalize_in, symmetric_method=config.symmetric_method,
                                    check_around_zero=config.check_around_zero, wcs=frame.wcs,
                                    zmin=zmin, zmax=zmax, soft_zmin=soft_zmin, soft_zmax=soft_zmax)
#print(vmin, vmax)

# -----------------------------------------------------------------

# Set the interval
interval = (vmin, vmax)

# -----------------------------------------------------------------

# Show the distribution of pixel values
if config.distribution:

    # Create distribution
    distribution = Distribution.from_data(name, frame, sigma_clip=config.sigma_clip, min_value=vmin, max_value=vmax)

    # Show the distribution
    plot_distribution(distribution, x_limits=interval)

# -----------------------------------------------------------------

# Adjust the frame for plotting
if config.symmetric or config.around_zero: raise ValueError("Cannot specifiy 'symmetric' or 'around_zero' if absolute values are taken")
if config.no_positives: raise ValueError("Cannot enable 'no_positives' if absolute values are taken")
regularize_frame(frame, absolute=config.absolute, dilate_nans=config.dilate_nans, dilation_radius=config.dilation_radius,
                 dilation_niterations=config.dilation_niterations, no_nans=config.no_nans, no_infs=config.no_infs,
                 no_negatives=config.no_negatives, no_positives=config.no_positives, cutoff_above=config.cutoff_above,
                 cutoff_below=config.cutoff_below)

# -----------------------------------------------------------------

# Plot
plotting.plot_map(frame, contours=config.contours, interval=interval, scale=config.scale, colorbar=config.colorbar,
                  cmap=config.cmap, ncontours=config.ncontours, contours_color=config.contours_color, path=config.path)

# -----------------------------------------------------------------
