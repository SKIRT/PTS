#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.polygonmovie Create some movies for a particular set of polygon data.
#
# This script creates some movies for a particular set of polygon data.
# The in/out filenames and other parameters are hardcoded in the script.
# Thus the script mainly serves as an example of how to use the pts.polygonmovie module.
#

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.prep.polygonmovie import *

# -----------------------------------------------------------------

def optim(x,y):
    return (16*x*(1-x)*y*(1-y)*np.sin(np.pi*x)*np.sin(np.pi*y)) ** 2

createpolygonmovie("polygon.dat", "polygon1.mov", rate=5, shape=(800,800),
                    background=optim)

createpolygonmovie("polygon.dat", "polygon2.mov", rate=5, shape=(800,800),
                    background=optim, xrange=(0.49,0.51), yrange=(0.49,0.51))

createpolygonmovie("polygon.dat", "polygon3.mov", rate=5, shape=(800,800),
                    background=optim, xrange=(0.499,0.501), yrange=(0.499,0.501))


createpolygonmovie("polygon.dat", "polygon1.mpg", shape=(800,600),
                    background=optim)

createpolygonmovie("polygon.dat", "polygon2.mpg", shape=(800,600),
                    background=optim, xrange=(0.49,0.51), yrange=(0.49,0.51))

createpolygonmovie("polygon.dat", "polygon3.mpg", shape=(800,600),
                    background=optim, xrange=(0.499,0.501), yrange=(0.499,0.501))

# -----------------------------------------------------------------
