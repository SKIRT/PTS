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

class Star(object):
    
    def __init__(self):
        
        self.radius = None
        self.amplitude = None
        self.flux = None
        self.center = None
        
    # *****************************************************************
    
    def plot(self):
        
        pass
        
# *****************************************************************