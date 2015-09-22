#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np

# Import astronomical units
from astropy import units as u

# *****************************************************************

class Box(object):
    
    def __init__(self):
        
        self.data = None
        self.background = None
        self.x_shift = None
        self.y_shift = None
        
    # *****************************************************************
    
    def plot(self):
        
        pass
        
# *****************************************************************
