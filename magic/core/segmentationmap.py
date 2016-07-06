#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.segmentationmap Contains the SegmentationMap class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from scipy import ndimage

# Import the relevant PTS classes and modules
from ..basics.vector import Position, Extent
from ..basics.geometry import Rectangle
from ..tools import cropping, fitting, interpolation, plotting, statistics
from ...core.tools.logging import log
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

class SegmentationMap(Frame):
    
    

# -----------------------------------------------------------------
