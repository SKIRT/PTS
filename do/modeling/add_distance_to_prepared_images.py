#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.add_distance_to_prepared_images Set the galaxy distance in the header of each prepared image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties, get_prepared_dataset
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Initialize log
log.start("Starting add_distance_to_prepared_images ...")

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Get galaxy distance
properties = get_galaxy_properties(modeling_path)
distance = properties.distance

# -----------------------------------------------------------------

# Loop over all prepared images
dataset = get_prepared_dataset(modeling_path)

# -----------------------------------------------------------------

# Loop over the images
for name in dataset.names:

    # Inform
    log.info("Adding the distance to the " + name + " image ...")
    
    # Load the image
    image = dataset.get_image(name)
    
    # Set the distance
    image.distance = distance
    
    # Save
    image.save()

# -----------------------------------------------------------------