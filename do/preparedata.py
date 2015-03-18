#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.preparedata Prepare data for use with SKIRT or FitSKIRT
#

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import numpy as np

# Import relevant PTS modules
from pts.image import Image
from pts.log import Log

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Get the full path to the data directory
datadir = os.path.join(os.getcwd(), "Data")


## 0. READING IMAGES

log.info("Reading images...")

# Make a list of the different images
images = []

for filename in os.listdir(datadir):

    # Check whether this file is a fits file and not hidden
    if filename.endswith(".fits") and not filename.startswith("."):

        filepath = os.path.join(datadir, filename)
        images.append(Image(filepath))

for image in images:

    image.maskregions()

## 1. SKY-SUBTRACTION

log.info("Sky-subtraction...")

# Subtract sky for each image
for image in images:

    # Check whether the image has already been sky-subtracted first
    if not image.subtracted: image.subtractsky()


## 2. UNIT CONVERSIONS

log.info("Unit conversions...")


# -----------------------------------------------------------------

## OLD STARS: IRAC, YOUNG NI STARS: FUV, YOUNG I STARS: Ha + 24micron, DUST: H + PACS70 + PACS160

# -----------------------------------------------------------------






