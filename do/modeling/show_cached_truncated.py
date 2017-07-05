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
from pts.core.tools import time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.core.remote.remote import Remote
from pts.core.tools.parsing import real
from pts.core.tools.introspection import pts_temp_dir
from pts.magic.core.frame import Frame
from pts.magic.tools import plotting
from pts.modeling.core.steps import cached_directory_name_for_single_command

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
definition.add_required("filter", "filter", "filter for which to show the truncated image")
definition.add_required("factor", "real", "truncation ellipse factor for which to plot")

# Get configuration
config = parse_arguments("plot_truncated", definition)

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Setup the remote
remote = Remote(host_id=environment.cache_host_id)

# -----------------------------------------------------------------

directory_name = cached_directory_name_for_single_command(environment, "truncate")
remote_truncation_path = fs.join(remote.home_directory, directory_name)

# -----------------------------------------------------------------

filter_name = str(config.filter)
remote_truncation_path_filter = fs.join(remote_truncation_path, filter_name)
if not remote.is_directory(remote_truncation_path_filter): raise ValueError("Could not find cached data for the " + filter_name + " image")

# -----------------------------------------------------------------

# Initialize variable
the_image_path = None

# Find the truncated image with the factor corresponding to the specified factor
for image_path, image_name in remote.files_in_path(remote_truncation_path_filter, extension="fits", returns=["path", "name"]):

    #print(image_path, image_name)

    # Determine the truncation factor
    factor = real(image_name)

    # Check the factor
    if np.isclose(factor, config.factor, rtol=0.01):
        the_image_path = image_path
        break

# Check
if the_image_path is None: raise ValueError("Could not find the truncated " + filter_name + " image with a truncation factor of " + str(config.factor))

# -----------------------------------------------------------------

# Download the image to a temporary path
filename = time.unique_name("truncated_" + filter_name + "_" + str(config.factor))
filepath = remote.download_file_to(the_image_path, pts_temp_dir, new_name=filename)

# -----------------------------------------------------------------

# Open the image
frame = Frame.from_file(filepath)

# -----------------------------------------------------------------

# Plot the truncated image
plotting.plot_box(frame)

# -----------------------------------------------------------------
