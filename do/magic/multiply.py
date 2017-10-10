#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.multiply Multiply a FITS image with a certain factor

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.image import Image
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_required("file_path", "file_path", "path of the image")
definition.add_required("factor", "real", "multiplication factor")
definition.add_flag("backup", "make a backup of each image that is multiplied", False)
config = parse_arguments("multiply", definition)

# -----------------------------------------------------------------

# BACKUP FIRST
if config.backup:

    # Inform the user
    log.info("Making a backup of the original image ...")

    # Determine new filepath and copy
    new_filepath = fs.appended_filepath(config.filepath, "_backup")
    fs.copy_file(config.filepath, new_filepath)

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the image '" + config.filepath + "' ...")

# Load the image
image = Image.from_file(config.filepath)

# -----------------------------------------------------------------

# Debugging
sum_before = np.nansum(image.primary.data)
log.debug("Sum of the primary frame before multiplication is " + str(sum_before))

# -----------------------------------------------------------------

# Inform th euser
log.info("Multiplying the image with a factor of " + str(config.factor) + " ...")

# Multiply
image *= config.factor

# -----------------------------------------------------------------

# Debugging
sum_after = np.nansum(image.primary.data)
log.debug("Sum of the primary frame after multiplication is " + str(sum_after))

# -----------------------------------------------------------------

# CHECK
if not np.isclose(sum_before * config.factor, sum_after): raise RuntimeError("Something went wrong: " + str(sum_before) + ", " + str(sum_after) + ", " + str(config.factor))

# -----------------------------------------------------------------

# Inform the user
log.info("Saving the image ...")

# Save
image.save()

# -----------------------------------------------------------------
