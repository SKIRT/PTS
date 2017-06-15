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
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties

# -----------------------------------------------------------------

# Initialize log
log.start("Starting add_distance_to_prepared_images ...")

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Get galaxy distance
properties = get_galaxy_properties(modeling_path)
distance = properties.distance

# -----------------------------------------------------------------

# Loop over all prepared images


# -----------------------------------------------------------------