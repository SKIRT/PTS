#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.geometry Basic 3D geometry calculations
#
# The functions and classes in this module offer some basic 3D geometry calculations.

# -----------------------------------------------------------------

# Import standard modules
import numpy as np
import math

# -----------------------------------------------------------------

## This function returns the norm of a vector, specified as a sequence of 3 coordinates.
def norm(vector):
    x,y,z = vector
    return math.sqrt(x*x+y*y+z*z)

## This function returns the unit vector in the same direction as an arbitrary vector,
# specified as a sequence of 3 coordinates.
def normalize(vector):
    x,y,z = vector
    norm = math.sqrt(x*x+y*y+z*z)
    return np.array((x/norm, y/norm, z/norm), np.float)

## This function returns the cross product (vector product) of two arbitrary vectors,
# each specified as a sequence of 3 coordinates.
def cross(vector1, vector2):
    u,v,w = vector1
    x,y,z = vector2
    return np.array((v*z-w*y, w*x-u*z, u*y-v*x), np.float)

## This function rotates a 3D position about an arbitrary rotation axis over some angle.
# The function arguments are:
#  - position: the coordinates of the position to be transformed as a 3-sequence.
#  - axispoint: the coordinates of a point on the rotation axis, as a 3-sequence.
#  - axisdirection: the coordinates of a direction vector for the positive rotation axis, as a 3-sequence.
#  - angle: the rotation angle in degrees; positive values indicate clockwise rotation
#    when looking down from the positive rotation axis.
#
def rotate(position, axispoint=(0.,0.,0.), axisdirection=(0.,0.,1.), angle=0.):
    # See "Principles of Interactive Computer Graphics", second edition, page 346
    tf = Transform();
    tf.translate(-axispoint[0], -axispoint[1], -axispoint[2])
    a,b,c = normalize(axisdirection)
    v = np.sqrt(b*b+c*c)
    if v > 0.3:
        tf.rotateX(c/v, -b/v)
        tf.rotateY(v, -a)
        tf.rotateZ(math.cos(angle*math.pi/180.), math.sin(angle*math.pi/180.))
        tf.rotateY(v, a)
        tf.rotateX(c/v, b/v)
    else:
        v = np.sqrt(a*a+c*c)
        tf.rotateY(c/v, -a/v)
        tf.rotateX(v, -b)
        tf.rotateZ(math.cos(angle*math.pi/180.), math.sin(angle*math.pi/180.))
        tf.rotateX(v, b)
        tf.rotateY(c/v, a/v)
    tf.translate(axispoint[0], axispoint[1], axispoint[2])
    return np.array(tf.transform(position[0], position[1], position[2], 1.)[0:3], np.float)

# -----------------------------------------------------------------

## This class represents a homogeneous coordinate transform, specified by a 4x4 matrix.
#
class Transform:

    ## This constructor creates an identity transform.
    def __init__(self):
        self.T = np.matrix(np.identity(4))

    ## This function adds a translation to the transformation.
    def translate(self, x, y, z):
        M = np.matrix(np.identity(4))
        M[3,0] = x
        M[3,1] = y
        M[3,2] = z
        self.T = self.T * M

    ## This function adds a scaling to the transformation.
    def scale(self, x, y, z):
        M = np.matrix(np.identity(4))
        M[0,0] = x
        M[1,1] = y
        M[2,2] = z
        self.T = self.T * M

    ## This function adds a rotation about the x axis to the transformation.
    def rotateX(self, cos, sin):
        M = np.matrix(np.identity(4))
        M[1,1] = cos
        M[2,2] = cos
        M[1,2] = -sin
        M[2,1] = sin
        self.T = self.T * M

    ## This function adds a rotation about the y axis to the transformation.
    def rotateY(self, cos, sin):
        M = np.matrix(np.identity(4))
        M[0,0] = cos
        M[2,2] = cos
        M[0,2] = -sin
        M[2,0] = sin
        self.T = self.T * M

    ## This function adds a rotation about the z axis to the transformation.
    def rotateZ(self, cos, sin):
        M = np.matrix(np.identity(4))
        M[0,0] = cos
        M[1,1] = cos
        M[0,1] = -sin
        M[1,0] = sin
        self.T = self.T * M

    ## This function adds perspective along the z axis to the transformation.
    def perspectiveZ(self, f):
        M = np.matrix(np.identity(4))
        M[2,3] = 1/f
        M[3,2] = -f
        M[3,3] = 0
        self.T = self.T * M

    ## This function returns the transformed coordinates for the specified input coordinates.
    def transform(self, x, y, z, w):
        P = np.matrix(( x, y, z, w )) * self.T
        return (P[0,0], P[0,1], P[0,2], P[0,3])

    ## This function returns the transformed coordinates for the specified sets input coordinates,
    # specified as numpy arrays.
    # TODO: improve efficiencey via either np.vectorize or explicit elementwise matrix multiplication
    def transform_vec(self, x, y, z, w):
        n = x.size
        crds = np.zeros((n, 4))
        for i in range(n):
            P = np.matrix(( x[i], y[i], z[i], w[i] )) * self.T
            crds[i] = np.asarray(P)[0,:]
        return crds[:,:3], crds[:,3]

# -----------------------------------------------------------------
