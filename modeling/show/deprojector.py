#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.show.deprojector Contains the Deprojector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import ShowComponent
from ...core.tools.logging import log

# -----------------------------------------------------------------

class Deprojector(ShowComponent):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        super(Deprojector, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        pass

### double x,y,z;
    #bfr.cartesian(x,y,z);

    #// Project and rotate the x and y coordinates
    #project(x);
    #rotate(x,y);

    #// Find the corresponding pixel in the image
    #int i = static_cast<int>(floor(x-_xmin)/_deltay);
    #int j = static_cast<int>(floor(y-_ymin)/_deltay);
    #if (i<0 || i>=_Nx || j<0 || j>=_Ny) return 0.0;

    #// Return the density
    #return _image(i,j) * exp(-fabs(z)/_hz) /(2*_hz)/(_deltax*_deltay);

### void
#ReadFitsGeometry::rotate(double &x, double &y)
#const
#{
#    // Cache the original values of x and y
#    double xorig = x;
#    double yorig = y;
#
#    // Calculate the coordinates in the plane of the image
#    x = (_sinpa * xorig)  + (_cospa * yorig);
#    y = (-_cospa * xorig) + (_sinpa * yorig);
#}

#////////////////////////////////////////////////////////////////////
#
#void
#ReadFitsGeometry::derotate(double &x, double &y)
#const
#{
#    // Cache the original values of x and y
#    double xorig = x;
#    double yorig = y;
#
#    // Calculate the coordinates in the rotated plane
#    x = (_sinpa * xorig) - (_cospa * yorig);
#    y = (_cospa * xorig) + (_sinpa * yorig);
#}

#////////////////////////////////////////////////////////////////////

#void
#ReadFitsGeometry::project(double &x)
#const
#{
#    x = x*_cosi;
#}

#////////////////////////////////////////////////////////////////////

#void
#ReadFitsGeometry::deproject(double &x)
#const
#{
#    x = x/_cosi;
#}

#////////////////////////////////////////////////////////////////////

### // Calculate the boundaries of the image in physical coordinates
#    _xmax = ((_Nx-_xc)*_deltay);
#    _xmin = -_xc*_deltay;
#    _ymax = ((_Ny-_yc)*_deltay);
#    _ymin = -_yc*_deltay;

#    // Calculate the sines and cosines of the position angle and inclination
#    _cospa = cos(_positionangle - (M_PI/2.0));
#    _sinpa = sin(_positionangle - (M_PI/2.0));
#    _cosi = cos(_inclination);
#    _sini = sin(_inclination);

#    // Calculate the physical pixel size in the x direction of the galactic plane
#    _deltax = _deltay / _cosi;

# -----------------------------------------------------------------
