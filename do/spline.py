#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.spline Experiment with splines.
#

# -----------------------------------------------------------------

# import standard modules
import sys
import os.path
import numpy as np

# import relevant pts modules
from pts.rgbimage import CubicSpline
import matplotlib.pyplot as plt

# -----------------------------------------------------------------

print "Starting..."

def plos(p1,p2):
    curve = CubicSpline(p1,p2)
    x = np.arange(0,1.005,0.01)
    y = np.arange(0,1.005,0.01)
    for index in range(len(x)): y[index] = curve.y(x[index])
    plt.figure()
    plt.plot((0,1),(0,1))
    plt.plot(x,y)
    plt.show()

plos( (0.25,0.16),(0.80,0.86) )

print "Finished."



# -----------------------------------------------------------------
