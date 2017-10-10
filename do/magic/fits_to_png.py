#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.fits_to_png Convert a FITS file to a PNG image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.basics.log import log
from pts.magic.core.rgba import alpha_methods

# -----------------------------------------------------------------

default_alpha_method = "absolute"
scales = ["log", "sqrt"]

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "path of the FITS file")
definition.add_optional("scale", "string", "scaling", "log", scales)
definition.add_optional("interval", "string", "interval", "pts")
definition.add_optional("colours", "string", "colour or colour scale", "red")
definition.add_optional("alpha", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_flag("no_alpha", "no alpha", False)
definition.add_optional("output", "string", "output filepath", letter="o")
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.)
definition.add_optional("max_npixels", "positive_integer", "maximum number of pixels")
definition.add_optional("downsample", "positive_real", "downsample with this factor")
definition.add_flag("show", "show after creating", False)
config = parse_arguments("fits_to_png", definition)

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the FITS file ...")

# Load the FITS file
frame = Frame.from_file(config.filename)

# -----------------------------------------------------------------

if config.output is not None: filepath = fs.absolute_or_in(config.output, fs.cwd())
else:

    # Determine the path
    name = fs.strip_extension(fs.name(config.filename))
    filepath = fs.absolute_path(name + ".png")

# -----------------------------------------------------------------

# Max npixels
if config.max_npixels is not None:

    # Determine downsample factor
    if frame.xsize > config.max_npixels or frame.ysize > config.max_npixels:

        factor = max(frame.xsize, frame.ysize) / float(config.max_npixels)
        # Make integer
        #numbers.round_up_to_int(factor)
        #if factor < 2.: factor = 2
        config.downsample = factor

# Downsample
if config.downsample is not None:

    # Inform the user
    log.info("Downsampling the image ...")

    # Downsample
    frame.downsample(config.downsample)

# -----------------------------------------------------------------

# Inform the user
log.info("Converting the image to a png file ...")

if config.no_alpha: alpha = None
else: alpha = config.alpha

# Save as PNG image
frame.saveto_png(filepath, interval=config.interval, scale=config.scale, alpha=alpha, peak_alpha=config.peak_alpha, colours=config.colours)

# -----------------------------------------------------------------

# Open the file
if config.show: fs.open_file(filepath)

# -----------------------------------------------------------------
