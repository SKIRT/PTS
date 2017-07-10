#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.interpolation This module is used for interpolating masked areas in an image.
#   Parts of it are adapted from the "lib.pyx" module, which can be found at
#   https://github.com/gasagna/openpiv-python/blob/master/openpiv/src/lib.pyx
#   written by Davide Lasagna, with slight changes applied described at
#   http://astrolitterbox.blogspot.be/2012/03/healing-holes-in-arrays-in-python.html

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

def inpaint_biharmonic(frame, mask):

    """
    This function ...
    :param frame:
    :param mask:
    :return:
    """

    from skimage.restoration import inpaint

    maximum = np.nanmax(frame)
    normalized = frame / maximum
    data = inpaint.inpaint_biharmonic(normalized, mask, multichannel=False)

    return data * maximum

# -----------------------------------------------------------------

# TODO: create a better inpainting function. OpenCV has one, but this is a terrible dependency because it's hard to
# install. Options:
#  - The below replace_nans function can be replaced by the more up to date version at:
#    https://github.com/OpenPIV/openpiv-python/blob/master/openpiv/src/lib.pyx
#    We may want to keep it in cython so that it runs faster. However, this original does not have the inverse distance
#    weighing as in the code below, but we can maybe add this ourselves in the cython code
#  - Write our own code.
# SOLUTION: SEE FUNCTION ABOVE, GENERALLY, IT IS MUCH BETTER

def in_paint(data, mask, method="localmean"):

    """
    This function ...
    :param data:
    :param mask:
    :param method:
    :return:
    """

    # Fill the data with nans according to the mask
    #data_ma = np.ma.array(data.astype(float), mask=mask)
    #data_nans = data_ma.filled(np.NaN)

    data_with_nans = np.copy(data)
    data_with_nans[mask] = np.NaN

    interpolated = replace_nans(data_with_nans, 5, 0.5, 2, method)

    # If the interpolated box contains nans, do not fill in the corresponding pixels of the data with these nans,
    # therefore set the pixels that are nan to False in the box_mask (take the difference between the box_mask
    # and the np.isnan(interpolated_box) mask). Then, set the nans to zero in the interpolated_box because
    # False * nan would otherwise still equal to nan.
    #box_mask = masks.subtract(box_mask, np.isnan(interpolated_box))
    #interpolated_box[np.isnan(interpolated_box)] = 0.0

    interpolated[np.isnan(interpolated)] = data[np.isnan(interpolated)]

    # Return the interpolated data
    return interpolated

# -----------------------------------------------------------------

def replace_nans(array, max_iter, tol, kernel_size=1, method='localmean'):

    """
    This function is used to replace NaN elements in an array using an iterative image inpainting algorithm.
    The algorithm is the following:
    1) For each element in the input array, replace it by a weighted average
       of the neighbouring elements which are not NaN themselves. The weights depends
       of the method type. If ``method=localmean`` weight are equal to 1/( (2*kernel_size+1)**2 -1 )

    2) Several iterations are needed if there are adjacent NaN elements.
      If this is the case, information is "spread" from the edges of the missing
      regions iteratively, until the variation is below a certain threshold.

    This function takes the following arguments:
    :param array: 2d np.ndarray; an array containing NaN elements that have to be replaced
    :param max_iter: int; the number of iterations
    :param tol:
    :param kernel_size: int; the size of the kernel, default is 1
    :param method: the method used to replace invalid values. Valid options are "localmean" or "idw" (inverse distance weighing).
    :return:
    """

    # Initialize arrays
    #filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)

    filled = np.zeros_like(array)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64)

    # Indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )

    # Number of NaN elements
    n_nans = len(inans)

    # Arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)

    # Depending on kernel type, fill kernel array
    if method == 'localmean':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1.

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5],
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        #print kernel, 'kernel'

    else:
        raise ValueError("Method not valid. Should be one of 'localmean' and 'idw'")

    # Fill new array with input elements
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            filled[i,j] = array[i,j]

    # Make several passes
    # until we reach convergence
    for it in range(max_iter):

        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]

            # Initialize to zero
            filled[i,j] = 0.0
            n = 0

            # Loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):

                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:

                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :

                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:

                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1

            # Divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = np.nan

        # Check if mean square difference between values of replaced
        # elements is below a certain tolerance
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]

    return filled

# -----------------------------------------------------------------

def sincinterp(image, x, y, kernel_size=3 ):

    """
    This function re-samples an image at intermediate positions between pixels. It uses a cardinal interpolation
    formula which limits the loss of information in the resampling process. It uses a limited number of neighbouring
    pixels. The new image \f$im^+\f$ at fractional locations \f$x\f$ and \f$y\f$ is computed as:
    \f[
    im^+(x,y) = \sum_{i=-\mathtt{kernel\_size}}^{i=\mathtt{kernel\_size}} \sum_{j=-\mathtt{kernel\_size}}^{j=\mathtt{kernel\_size}} \mathtt{image}(i,j)  sin[\pi(i-\mathtt{x})]  sin[\pi(j-\mathtt{y})]  / \pi(i-\mathtt{x}) / \pi(j-\mathtt{y})
    \f]
    This function takes the following arguments:
    :param image: np.darray, dtype np.int32, the image array
    :param x: two dimensions np.ndarray of floats; an array containing fractional pixel row positions at which to interpolate the image
    :param y: two dimensions np.darray of floats; an array containing fractional pixel column positions at which to interpolate the image
    :param kernel_size: int, interpolation is performed over a (2*kernel_size+1)*(2*kernel_size+1) submatrix in the neighbourhood of each interpolation point.
    :return: im, an np.darray, dtype np.float64: the interpolated value of image at the points specified by x and y.
    """

    # The output array
    r = np.zeros( [x.shape[0], x.shape[1]], dtype=np.float64)

    # Fast pi
    pi = 3.1419

    # For each point of the output array
    for I in range(x.shape[0]):
        for J in range(x.shape[1]):

            #loop over all neighbouring grid points
            for i in range( int(x[I,J])-kernel_size, int(x[I,J])+kernel_size+1 ):
                for j in range( int(y[I,J])-kernel_size, int(y[I,J])+kernel_size+1 ):
                    # check that we are in the boundaries
                    if i >= 0 and i <= image.shape[0] and j >= 0 and j <= image.shape[1]:
                        if (i-x[I,J]) == 0.0 and (j-y[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j]
                        elif (i-x[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(j-y[I,J]) )/( pi*(j-y[I,J]) )
                        elif (j-y[I,J]) == 0.0:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(i-x[I,J]) )/( pi*(i-x[I,J]) )
                        else:
                            r[I,J] = r[I,J] + image[i,j] * np.sin( pi*(i-x[I,J]) )*np.sin( pi*(j-y[I,J]) )/( pi*pi*(i-x[I,J])*(j-y[I,J]))

    # Return the interpolated image
    return r

# -----------------------------------------------------------------
