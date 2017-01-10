#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.pointsource Contains the PointSource class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy
from skimage import morphology

# Import astronomical modules
from astropy.table import Table
from photutils import find_peaks
from photutils import detect_sources
from photutils import detect_threshold
from astropy.coordinates import Angle
from astropy import units as u
from astropy.convolution import convolve, convolve_fft

# Import the relevant PTS classes and modules
from .image import Image
from .frame import Frame
from .box import Box
from ..basics.mask import Mask
from ..region.ellipse import PixelEllipseRegion
from ..region.rectangle import PixelRectangleRegion
from ..region.circle import PixelCircleRegion
from ..tools import plotting, statistics
from ..basics.stretch import PixelStretch
from ..basics.coordinate import PixelCoordinate

# -----------------------------------------------------------------

class PointSource(object):

    """
    This class...
    """

    def __init__(self, position):

        """
        The constructor ...
        :param position:
        """

        # Set attributes
        self.position = None

# -----------------------------------------------------------------
