#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.geometry Contains geometry related classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class Ellipse(object):
    
    """
    This class ...
    """
    
    def __init__(self, center, radius, angle=0.0):

        """
        The constructor ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        self.center = center
        self.radius = radius
        self.angle = angle

# -----------------------------------------------------------------

class Rectangle(object):
    
    """
    This class ...
    """
    
    def __init__(self, lowerleft, upperright, angle=0.0):
        
        """
        The constructor ...
        :param lowerleft:
        :param upperright:
        :param angle:
        :return:
        """

        pass


        
# -----------------------------------------------------------------
