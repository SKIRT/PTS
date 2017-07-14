#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_truncated Plot the truncated images for a certain factor.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.magic.plot.imagegrid import StandardImageGridPlotter
from pts.core.filter.filter import parse_filter
from pts.core.tools.parsing import real
from pts.magic.core.frame import Frame
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
definition.add_required("factor", "real", "truncation ellipse factor for which to plot")

# Get configuration
config = parse_arguments("plot_truncated", definition)

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Create plotter
plotter = StandardImageGridPlotter()

# Loop over the directories in the truncation path
for path, name in fs.directories_in_path(environment.truncation_path, returns=["path", "name"]):

    # Determine the path to the lowres directory
    lowres_path = fs.join(path, "lowres")

    # Determine the filter
    fltr = parse_filter(name)
    filter_name = str(fltr)

    # Initializ variable
    the_image_path = None

    # Find the image corresponding to the specified factor
    for image_path, image_name in fs.files_in_path(lowres_path, extension="fits", returns=["path", "name"]):

        # Determine the factor
        factor = real(image_name)

        # If the factor corresponds to the specified factor, take this image
        if np.isclose(factor, config.factor, rtol=0.01):
            the_image_path = image_path
            break

    # Check
    if the_image_path is None: raise ValueError("No truncated " + filter_name + " image found for a factor of " + str(config.factor))

    # Add the image
    frame = Frame.from_file(the_image_path)
    plotter.add_image(frame, filter_name)

# -----------------------------------------------------------------

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
