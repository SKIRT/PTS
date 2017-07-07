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

# Import standard modules
import imageio
import numpy as np

# Import astronomical modules
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import MinMaxInterval, ZScaleInterval

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.magic.core.frame import Frame
from pts.core.tools import parsing

# -----------------------------------------------------------------

scales = ["log", "sqrt"]

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "path of the FITS file")
definition.add_optional("scale", "string", "scaling", "log", scales)
definition.add_optional("interval", "string", "interval", "pts")
definition.add_optional("colours", "string", "colour or colour scale", "red")
definition.add_flag("alpha", "use alpha", True)
definition.add_optional("output", "string", "output filepath", letter="o")
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.)
config = parse_arguments("fits_to_png", definition)

# -----------------------------------------------------------------

# Load the FITS file
frame = Frame.from_file(config.filename)

# -----------------------------------------------------------------

# INTERVAL
if config.interval == "zscale": vmin, vmax = ZScaleInterval().get_limits(frame.data)
elif config.interval == "pts":
    # Determine the maximum value in the box and the mimimum value for plotting
    #print("here")
    vmin = max(np.nanmin(frame.data), 0.)
    vmax = 0.5 * (np.nanmax(frame.data) + vmin)
elif config.interval == "minmax": vmin, vmax = MinMaxInterval().get_limits(frame.data)
else:
    try: vmin, xmax = parsing.real_tuple(config.interval)
    except ValueError: raise ValueError("Cannot interpret the interval")

# Normalization
if config.scale == "log": norm = ImageNormalize(stretch=LogStretch(), vmin=vmin, vmax=vmax)
elif config.scale == "sqrt": norm = ImageNormalize(stretch=SqrtStretch(), vmin=vmin, vmax=vmax)
else: raise ValueError("Invalid option for 'scale'")

# -----------------------------------------------------------------

# Set min and max for normalizer
#norm.vmin = vmin
#norm.vmax = vmax

# -----------------------------------------------------------------

frame.replace_nans(0.0)
frame.replace_infs(0.0)

# -----------------------------------------------------------------

# Normalize
normalized = norm(frame.data)

# REmove nans and infinites
#frame.data[np.isnan(normalized)] = 0.0
#normalized[np.isinf(normalized)] = 0.0

# Determine transparency
if config.alpha: transparency = config.peak_alpha * normalized / np.nanmax(normalized)
else: transparency = np.ones_like(normalized)

# -----------------------------------------------------------------

# CREATE THE CHANNEL ARRAYS

# Red image
if config.colours == "red":

    # NxMx4
    red = normalized * 255
    blue = np.zeros_like(red)
    green = np.zeros_like(red)
    alpha = transparency * 255

# Blue image
elif config.colours == "blue":

    red = np.zeros_like(normalized)
    blue = normalized * 255
    green = np.zeros_like(red)
    alpha = transparency * 255

# Green image
elif config.colours == "green":

    red = np.zeros_like(normalized)
    blue = np.zeros_like(normalized)
    green = normalized * 255
    alpha = transparency * 255

else: raise NotImplementedError("Not impl")

# -----------------------------------------------------------------

# MAKE THE IMAGE ARRAY
# Stack, create the image array
arrays = [red, green, blue, alpha]
image = np.stack(arrays, axis=-1)

# APPLY COLOUR MAP:
# A = self.cmap(A, alpha=self.get_alpha(), bytes=True)

# -----------------------------------------------------------------

if config.output is not None: filepath = fs.absolute_or_in(config.output, fs.cwd())
else:

    # Determine the path
    name = fs.strip_extension(fs.name(config.filename))
    filepath = fs.absolute_path(name + ".png")

# -----------------------------------------------------------------

# Write
imageio.imwrite(filepath, image)

# -----------------------------------------------------------------

# Open the file
fs.open_file(filepath)

# -----------------------------------------------------------------
