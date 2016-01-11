#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.apertures Contains functions for dealing with apertures.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from ..basics import Position

# -----------------------------------------------------------------

def position(aperture):
    
    """
    This function ...
    """

    x, y = aperture.positions[0]
    return Position(x, y)

# -----------------------------------------------------------------
    
def ellipticity(aperture):
    
    """
    This function ...
    """

    # Calculate the ellipticity
    return (aperture.a - aperture.b) / aperture.a

# -----------------------------------------------------------------
