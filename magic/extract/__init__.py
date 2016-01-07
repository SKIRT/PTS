#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.magic.extract TO DO
#
# This package ...
#

# -----------------------------------------------------------------

# Import classes to make them available at the level of this subpackage
from .extraction import Extractor
from .galaxyextraction import GalaxyExtractor
from .starextraction import StarExtractor
from .trainedextractor import TrainedExtractor
