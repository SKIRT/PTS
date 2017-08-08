#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.viewangles Interactive preview of instrument angles
#
# The class in this module provides an interactive preview of a simple 3D scene
# as if recorded through a SKIRT instrument with certain specified attributes.

# -----------------------------------------------------------------

# Import standard modules
import warnings
import numpy as np
import math

# Import the relevant PTS classes and modules
from ..basics.log import log

# Use an interactive back-end that supports animation
import matplotlib
with warnings.catch_warnings():
    warnings.filterwarnings('error')
    try:
        if matplotlib.get_backend().lower() != "tkagg": matplotlib.use("tkagg")
    except Warning as w: log.warning("An failed attempt of setting the Matplotlib backend has been made because it has already been set")
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ..tools.geometry import Transform

# -----------------------------------------------------------------

## This class provides an interactive preview of a simple 3D scene as if recorded through
# a SKIRT instrument with certain specified attributes. It supports both distant instruments
# (with given inclination, azimuth and position angles) and perspective instruments
# (with given viewport origin, crosshair and upwards positions, focal length, and viewport width).
#
# The constructor creates a plot window displaying the scene with a default instrument.
# Various member functions update the scene to a new instrument, optionally using animation
# to simulate the camera motion for a fly-by movie (see the flybymovie module).
#
class ViewAngles:

    ## The constructor constructs a scene representing the three axes and a flat disk in the xy plane,
    # and then creates a plot window displaying the scene using a distant instrument in the first octant.
    def __init__(self, figsize=(8,8)):

        # construct our scene
        self.scene = Scene()

        # setup the figure, using a distant instrument with default angles
        plt.ion()
        self.figure = plt.figure(figsize=figsize)
        plt.xlim(-1.25,1.25)
        plt.ylim(-1.25,1.25)
        self.scene.plot(DistantInstrument(55,45,0))
        plt.show()

    ## This function updates the scene in the plot window to reflect the distant instrument
    # with the specified inclination, azimuth and position angles (in degrees).
    def distant(self, inclination, azimuth, positionAngle):
        # remove any previous "lines" from the plot
        while len(self.figure.axes[0].lines)>0: self.figure.axes[0].lines.pop()
        # plot the scene and update the window
        self.scene.plot(DistantInstrument(inclination, azimuth, positionAngle))
        plt.draw()

    ## This function updates the scene in the plot window to reflect the perspective instrument with the
    # specified viewport origin, crosshair position, upwards position, focal length, and viewport width.
    # The first three attributes are given as (x,y,z) tuples in world coordinates.
    def perspective(self, viewport, crosshair, upwards, focal, width=1.):
        # remove any previous "lines" from the plot
        while len(self.figure.axes[0].lines)>0: self.figure.axes[0].lines.pop()
        # plot the scene and update the window
        self.scene.plot(PerspectiveInstrument(viewport, crosshair, upwards, focal, width))
        plt.draw()

    ## This function updates the scene in steps, moving through a list of frames obtained from an instance
    # of the Timeline class. For good results, the timeline's viewport should be a square sized around 1.
    def flyby(self, timeline):
        delay = 1. / timeline.getrate()
        for shape,size,viewport,crosshair,upwards,focal in timeline.getframes():
            self.perspective(viewport,crosshair,upwards,focal,size[0])
            plt.pause(delay)

    ## This function closes the figure.
    def close(self):
        plt.close()

# -----------------------------------------------------------------

## This class represents a scene for use in the ViewAngles class. The constructor creates
# an object containing the scene and offering functions to plot it using a given instrument.
#
class Scene:

    ## This constructor constructs a scene representing the three axes and an asymmetric figure.
    # All points fit inside the unit cube.
    def __init__(self):

        # initially empty lists for points to be painted black, red, green or blue.
        self.k = []
        self.r = []
        self.g = []
        self.b = []

        # origin
        self.k += [ (0.0, 0.0, 0.0) ]
        # x-axis
        self.g += [ (0.1, 0.0, 0.0) ]
        self.g += [ (0.2, 0.0, 0.0) ]
        self.g += [ (0.3, 0.0, 0.0) ]
        self.g += [ (0.4, 0.0, 0.0) ]
        self.g += [ (0.5, 0.0, 0.0) ]
        self.g += [ (0.6, 0.0, 0.0) ]
        self.g += [ (0.7, 0.0, 0.0) ]
        self.g += [ (0.8, 0.0, 0.0) ]
        self.g += [ (0.9, 0.0, 0.0) ]
        self.r += [ (-0.1, 0.0, 0.0) ]
        self.r += [ (-0.2, 0.0, 0.0) ]
        self.r += [ (-0.3, 0.0, 0.0) ]
        self.r += [ (-0.4, 0.0, 0.0) ]
        self.r += [ (-0.5, 0.0, 0.0) ]
        self.r += [ (-0.6, 0.0, 0.0) ]
        self.r += [ (-0.7, 0.0, 0.0) ]
        self.r += [ (-0.8, 0.0, 0.0) ]
        self.r += [ (-0.9, 0.0, 0.0) ]
        # y-axis
        self.g += [ (0.0, 0.1, 0.0) ]
        self.g += [ (0.0, 0.2, 0.0) ]
        self.g += [ (0.0, 0.3, 0.0) ]
        self.g += [ (0.0, 0.4, 0.0) ]
        self.g += [ (0.0, 0.5, 0.0) ]
        self.g += [ (0.0, 0.6, 0.0) ]
        self.r += [ (0.0, -0.1, 0.0) ]
        self.r += [ (0.0, -0.2, 0.0) ]
        self.r += [ (0.0, -0.3, 0.0) ]
        self.r += [ (0.0, -0.4, 0.0) ]
        self.r += [ (0.0, -0.5, 0.0) ]
        self.r += [ (0.0, -0.6, 0.0) ]
        # z-axis
        self.g += [ (0.0, 0.0, 0.1) ]
        self.g += [ (0.0, 0.0, 0.2) ]
        self.g += [ (0.0, 0.0, 0.3) ]
        self.r += [ (0.0, 0.0, -0.1) ]
        self.r += [ (0.0, 0.0, -0.2) ]
        self.r += [ (0.0, 0.0, -0.3) ]

        # a disk in the xy plane
        npoints = 50
        disk = np.zeros((npoints,3))
        for i in range(npoints):
            alpha = i*2*math.pi/npoints
            self.k += [ ( 0.9*math.cos(alpha), 0.6*math.sin(alpha), 0.0 ) ]

        # an asymmetric feature (floor)
        self.b += [ (-0.3, -0.2, 0.1) ]
        self.b += [ (-0.2, -0.2, 0.1) ]
        self.b += [ (-0.1, -0.2, 0.1) ]
        self.b += [ (-0.0, -0.2, 0.1) ]
        self.b += [ ( 0.1, -0.2, 0.1) ]
        self.b += [ ( 0.2, -0.2, 0.1) ]
        self.b += [ ( 0.3, -0.2, 0.1) ]
        self.b += [ ( 0.4, -0.2, 0.1) ]
        self.b += [ ( 0.5, -0.2, 0.1) ]
        self.b += [ ( 0.6, -0.2, 0.1) ]
        self.b += [ ( 0.7, -0.2, 0.1) ]
        self.b += [ ( 0.7, -0.1, 0.1) ]
        self.b += [ ( 0.7,  0.0, 0.1) ]
        self.b += [ ( 0.7,  0.1, 0.1) ]
        self.b += [ ( 0.7,  0.2, 0.1) ]
        self.b += [ ( 0.7,  0.3, 0.1) ]
        self.b += [ ( 0.7,  0.4, 0.1) ]
        self.b += [ ( 0.6,  0.4, 0.1) ]
        self.b += [ ( 0.5,  0.4, 0.1) ]
        self.b += [ ( 0.4,  0.4, 0.1) ]
        self.b += [ ( 0.3,  0.4, 0.1) ]
        self.b += [ ( 0.2,  0.4, 0.1) ]
        self.b += [ ( 0.1,  0.4, 0.1) ]
        self.b += [ ( 0.0,  0.4, 0.1) ]
        self.b += [ (-0.1,  0.4, 0.1) ]
        self.b += [ (-0.2,  0.4, 0.1) ]
        self.b += [ (-0.3,  0.4, 0.1) ]
        self.b += [ (-0.3,  0.3, 0.1) ]
        self.b += [ (-0.3,  0.2, 0.1) ]
        self.b += [ (-0.3,  0.1, 0.1) ]
        self.b += [ (-0.3,  0.0, 0.1) ]
        self.b += [ (-0.3, -0.1, 0.1) ]
        self.b += [ (-0.3, -0.2, 0.1) ]

        # an asymmetric feature (roof)
        self.b += [ (-0.2, -0.14, 0.2) ]
        self.b += [ (-0.1, -0.08, 0.3) ]
        self.b += [ (-0.0, -0.02, 0.4) ]
        self.b += [ ( 0.1, +0.04, 0.5) ]
        self.b += [ ( 0.2, +0.10, 0.6) ]
        self.b += [ ( 0.3, +0.16, 0.5) ]
        self.b += [ ( 0.4, +0.22, 0.4) ]
        self.b += [ ( 0.5, +0.28, 0.3) ]
        self.b += [ ( 0.6, +0.34, 0.2) ]
        self.b += [ ( 0.6, -0.14, 0.2) ]
        self.b += [ ( 0.5, -0.08, 0.3) ]
        self.b += [ ( 0.4, -0.02, 0.4) ]
        self.b += [ ( 0.3, +0.04, 0.5) ]
        self.b += [ ( 0.1, +0.16, 0.5) ]
        self.b += [ ( 0.0, +0.22, 0.4) ]
        self.b += [ (-0.1, +0.28, 0.3) ]
        self.b += [ (-0.2, +0.34, 0.2) ]

        # replace lists by numpy arrays
        self.k = np.array(self.k).T
        self.r = np.array(self.r).T
        self.g = np.array(self.g).T
        self.b = np.array(self.b).T

    ## This function plots the scene to the current figure using the specified instrument
    def plot(self, instrument):
        # transform and plot the points for each color in turn
        x,y = instrument.transform(self.k)
        plt.plot(x, y, 'k.')
        x,y = instrument.transform(self.r)
        plt.plot(x, y, 'r.')
        x,y = instrument.transform(self.g)
        plt.plot(x, y, 'go')
        x,y = instrument.transform(self.b)
        plt.plot(x, y, 'b*')

# -----------------------------------------------------------------

## This class represents a distant instrument and the corresponding coordinate transformation.
#
class DistantInstrument:

    ## This constructor creates a distant instrument with the specified inclination,
    # azimuth and position angles (in degrees).
    def __init__(self, inclination, azimuth, positionAngle):
        # calculate sin/cos of the angles
        self.costheta = math.cos(inclination*math.pi/180.)
        self.sintheta = math.sin(inclination*math.pi/180.)
        self.cosphi = math.cos(azimuth*math.pi/180.)
        self.sinphi = math.sin(azimuth*math.pi/180.)
        self.cospa = math.cos(positionAngle*math.pi/180.)
        self.sinpa = math.sin(positionAngle*math.pi/180.)

    ## This function transforms the specified 3D point to viewport coordinates for this
    # instrument, just like in SKIRT.
    def transform(self, points):
        x,y,z = points
        xpp = - self.sinphi*x + self.cosphi*y
        ypp = - self.cosphi*self.costheta*x - self.sinphi*self.costheta*y + self.sintheta*z
        xp = self.cospa * xpp - self.sinpa * ypp
        yp = self.sinpa * xpp + self.cospa * ypp
        return (xp, yp)

# -----------------------------------------------------------------

## This class represents a perspective instrument and the corresponding coordinate transformation.
#
class PerspectiveInstrument:

    ## This constructor creates a perspective instrument with the specified
    # viewport origin, crosshair position, upwards position, focal length, and viewport width.
    # The first three attributes are given as (x,y,z) tuples in world coordinates.
    def __init__(self, viewport, crosshair, upwards, focal, width=1):
        # parameters
        V = np.array(viewport, np.float)
        C = np.array(crosshair, np.float)
        U = np.array(upwards, np.float)
        F = float(focal)
        S = float(width)
        # derived positions and directions
        VC = V - C
        G = VC / np.sqrt(np.sum(VC*VC))
        E = V + F * G
        # the perspective transformation
        self.transf = Transform()
        # from world to eye coordinates
        self.transf.translate(-E[0], -E[1], -E[2])
        a,b,c = G
        v = np.sqrt(b*b+c*c)
        if v > 0.3:
            self.transf.rotateX(c/v, -b/v)
            self.transf.rotateY(v, -a)
            k = (b*b+c*c)*U[0] - a*b*U[1] - a*c*U[2]
            l = c*U[1] - b*U[2]
            u = np.sqrt(k*k+l*l)
            self.transf.rotateZ(l/u, -k/u)
        else:
            v = np.sqrt(a*a+c*c)
            self.transf.rotateY(c/v, -a/v)
            self.transf.rotateX(v, -b)
            k = c*U[0] - a*U[2]
            l = (a*a+c*c)*U[1] - a*b*U[0] - b*c*U[2]
            u = np.sqrt(k*k+l*l)
            self.transf.rotateZ(l/u, -k/u)
        self.transf.scale(1, 1, -1)
        # from eye to viewport coordinates
        self.transf.perspectiveZ(F)
        self.transf.scale(1./S, 1./S, 1)

    ## This function transforms the specified 3D points to viewport coordinates for this
    # instrument, just like in SKIRT.
    def transform(self, points):
        x,y,z = points
        xout = []
        yout = []
        for i in range(len(x)):
            xv, yv, zv, wv = self.transf.transform(x[i], y[i], z[i], 1)
            if zv > 0:
                xout += [xv / wv]
                yout += [yv / wv]
        return (np.array(xout), np.array(yout))

# -----------------------------------------------------------------
