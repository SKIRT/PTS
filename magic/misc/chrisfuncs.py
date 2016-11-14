#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.chrisfuncs Contains parts of the ChrisFuncs package 
#  (https://github.com/Stargrazer82301/ChrisFuncs)

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import pdb
import math
import numpy as np
import scipy.spatial
import scipy.ndimage

# -----------------------------------------------------------------

# IDIOT'S GUIDE TO ELLIPTICAL APERTURES
# I assume you also know the location of the ellipse's center. Call that (x0,y0).
# Let t be the counterclockwise angle the major axis makes with respect to the
# x-axis. Let a and b be the semi-major and semi-minor axes, respectively. If

# P = (x,y) is an arbitrary point then do this:
# X = (x-x0)*cos(t)+(y-y0)*sin(t); % Translate and rotate coords.
# Y = -(x-x0)*sin(t)+(y-y0)*cos(t); % to align with ellipse
# If
# X^2/a^2+Y^2/b^2
# is less than 1, the point P lies inside the ellipse. If it equals 1, it is right on
# the ellipse. If it is greater than 1, P is outside.

# -----------------------------------------------------------------

def EllipseCentre(a):

    """
    Function to calculate the coordinates of the centre of an ellipse produced by EllipseFit
    Input: Ellipse produced by EllipseFit
    Output: Array of x & y coordinates of ellipse centre
    """

    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    b,c,d,f,g,a = b,c,d,f,g,a
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

# -----------------------------------------------------------------

def EllipseAxes(a):

    """
    Function to calculate the lengths of the axes of an ellipse produced by EllipseFit
    Input: Ellipse produced by EllipseFit
    Output: Array of ellipse's major & minor axes
    """

    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

# -----------------------------------------------------------------

def EllipseAngle(a):

    """
    # Function to calculate the position angle of the centre of an ellipse produced by EllipseFit
    # Input: Ellipse produced by EllipseFit
    # Output: Ellipse's position angle
    """

    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    b,c,d,f,g,a = b,c,d,f,g,a
    return 0.5*np.arctan(2*b/(a-c))

# -----------------------------------------------------------------

def SigmaClip(values, tolerance=0.001, median=False, sigma_thresh=3.0, no_zeros=False):

    """
    Function to perform a sigma clip upon a set of values
    Input: Array of values, convergence tolerance, state if median instead of mean should be used for clip centrepoint, clipping threshold, boolean for whether sigma of zero can be accepted
    Returns: List containing the clipped standard deviation, the average, and the values themselves
    """

    # Remove NaNs from input values
    values = np.array(values)
    values = values[ np.where(np.isnan(values)==False) ]
    values_original = np.copy(values)

    # Continue loop until result converges
    diff = 10E10
    while diff > tolerance:

        # Assess current input iteration
        if median == False: average = np.mean(values)
        elif median == True: average = np.median(values)
        sigma_old = np.std(values)

        # Mask those pixels that lie more than 3 stdev away from mean
        check = np.zeros([len(values)])
        check[ np.where( values>(average+(sigma_thresh*sigma_old)) ) ] = 1
        check[ np.where( values<(average-(sigma_thresh*sigma_old)) ) ] = 1
        values = values[ np.where(check<1) ]

        # Re-measure sigma and test for convergence
        sigma_new = np.std(values)
        diff = abs(sigma_old-sigma_new) / sigma_old

    # Perform final mask
    check = np.zeros([len(values)])
    check[ np.where( values>(average+(sigma_thresh*sigma_old)) ) ] = 1
    check[ np.where( values<(average-(sigma_thresh*sigma_old)) ) ] = 1
    values = values[ np.where(check<1) ]

    # If required, check if calculated sigma is zero
    if no_zeros:

        if sigma_new==0.0:

            sigma_new = np.std(values_original)
            if median==False: average = np.mean(values)
            elif median==True: average = np.median(values)

    # Return results
    return [sigma_new, average, values]

# -----------------------------------------------------------------

def EllipseMask(array, rad, axial_ratio, angle, i_centre, j_centre):

    """
    Function to return a mask identifying all pixels within an ellipse of given parameters
    Input: Array, semi-major axis (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Output: Mask array of same dimensions as input array where pixels that lie within ellipse have value 1
    """

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj = float(rad)
    semi_min = float(rad) / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, array.shape[0]-1, array.shape[0])
    j_linespace = np.linspace(0, array.shape[1]-1, array.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check = (j_trans**2 / semi_maj**2) + (i_trans**2 / semi_min**2 )

    # Create ellipse mask
    ellipse_mask = np.zeros([array.shape[0], array.shape[1]])
    ellipse_mask[ np.where( ellipse_check<=1 ) ] = 1.0

    # Return array
    return ellipse_mask

# -----------------------------------------------------------------

def EllipseSum(array, rad, axial_ratio, angle, i_centre, j_centre):

    """
    Function to sum all elements in an ellipse centred on the middle of a given array
    Input: Array, semi-major axis (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Returns: Numpy array containing the sum of the pixel values in the ellipse, total number of pixels counted, and an array containing the pixel values
    """

    # Create slice of input array, containing only the region of interest
    i_cutout_min = int(np.floor(max([0, i_centre-rad])))
    i_cutout_max = int(np.ceil(min([(array.shape)[0], i_centre+rad])))
    j_cutout_min = int(np.floor(max([0, j_centre-rad])))
    j_cutout_max = int(np.ceil(min([(array.shape)[1], j_centre+rad])))
    array_slice = array[ int(round(i_cutout_min)):int(round(i_cutout_max))+1, int(round(j_cutout_min)):int(round(j_cutout_max))+1 ]
    i_centre_slice = i_centre - i_cutout_min
    j_centre_slice = j_centre - j_cutout_min
    if array[int(i_centre),int(j_centre)]!=array_slice[int(i_centre_slice),int(j_centre_slice)]:
        if np.isnan(array[int(i_centre),int(j_centre)]==False) and np.isnan(array_slice[int(i_centre_slice),int(j_centre_slice)]==False):
            print('SEVERE ERROR: EllipseSum check failed.')
            pdb.set_trace()
    else:
        array = array_slice
        i_centre = i_centre_slice
        j_centre = j_centre_slice

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj = float(rad)
    semi_min = float(rad) / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, array.shape[0]-1, array.shape[0])
    j_linespace = np.linspace(0, array.shape[1]-1, array.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check = (j_trans**2 / semi_maj**2) + (i_trans**2 / semi_min**2 )

    # Calculate flux & pixels in aperture, and store pixel values
    ellipse_where = np.where( (ellipse_check<=1) & (np.isnan(array)==False) )
    ellipse_tot = sum( array[ ellipse_where ] )
    ellipse_count = ellipse_where[0].shape[0]
    ellipse_pix = array[ ellipse_where ]
    ellipse_nan = np.where( (ellipse_check<=1) & (np.isnan(array)==True) )

    # Return results
    return [ellipse_tot, ellipse_count, ellipse_pix, ellipse_nan]

# -----------------------------------------------------------------

def EllipseFit(x,y):

    """
    Function that combines all of the ellipse-fitting steps (finds convex hull, fits ellipse to this, then finds properties of ellipse)
    Input: x & y coordinates to which the ellipse is to be fitted
    Output: Array of x & y coordinates of ellipse centre, array of ellipse's major & minor axes, ellipse's position angle
    """

    # Find convex hull of points
    p = np.zeros([x.shape[0],2])
    p[:,0], p[:,1] = x, y
    h = []
    for s in scipy.spatial.ConvexHull(p).simplices:
        h.append(p[s[0]])
        h.append(p[s[1]])
    h = np.array(h)
    x, y = h[:,0], h[:,1]

    # Carry out ellipse-fitting witchcraft
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]

    # Calculate and return properties of ellipse
    centre = np.real(EllipseCentre(a))
    axes = np.real(EllipseAxes(a))
    angle = (180/3.14159) * np.real(EllipseAngle(a))
    if axes[0]<axes[1]:
        angle += 90.0
    return np.array([centre, axes, angle, [x,y]])

# -----------------------------------------------------------------

def ToPrecision(x, p):

    """
    Function from https://github.com/randlet/to-precision
    returns a string representation of x formatted with a precision of p
    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0" * (p - 1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x / tens)

    if n < math.pow(10, p - 1):
        e = e - 1
        tens = math.pow(10, e - p + 1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n = n + 1

    if n >= math.pow(10, p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e + 1])
        if e + 1 < len(m):
            out.append(".")
            out.extend(m[e + 1:])
    else:
        out.append("0.")
        out.extend(["0"] * -(e + 1))
        out.append(m)

    return "".join(out)

# -----------------------------------------------------------------

def AnnulusSum(array, rad_inner, width, axial_ratio, angle, i_centre, j_centre):

    """
    Function to sum all elements in an annulus centred upon the middle of the given array
    Input: Array, semi-major axis of inside edge of annulus (pix), width of annulus (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Returns: Numpy array containing the sum of the pixel values in the annulus, the total number of pixels counted, and an array containing the pixel values
    """

    # Create slice of input array, containing only the region of interest
    i_cutout_min = int(np.floor(max([0, i_centre-(rad_inner+width)])))
    i_cutout_max = int(np.ceil(min([(array.shape)[0], i_centre+(rad_inner+width)])))
    j_cutout_min = int(np.floor(max([0, j_centre-(rad_inner+width)])))
    j_cutout_max = int(np.ceil(min([(array.shape)[1], j_centre+(rad_inner+width)])))
    array_slice = array[ int(round(i_cutout_min)):int(round(i_cutout_max))+1, int(round(j_cutout_min)):int(round(j_cutout_max))+1 ]
    i_centre_slice = i_centre - i_cutout_min
    j_centre_slice = j_centre - j_cutout_min
    if array[int(i_centre),int(j_centre)]!=array_slice[int(i_centre_slice),int(j_centre_slice)]:
        if np.isnan(array[int(i_centre),int(j_centre)]==False) and np.isnan(array_slice[int(i_centre_slice),int(j_centre_slice)]==False):
            print('SEVERE ERROR: AnnulusSum check failed.')
            pdb.set_trace()
    else:
        array = array_slice
        i_centre = i_centre_slice
        j_centre = j_centre_slice

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj_inner = float(rad_inner)
    semi_min_inner = float(semi_maj_inner) / float(axial_ratio)
    semi_maj_outer = float(rad_inner) + float(width)
    semi_min_outer  = float(semi_maj_outer) / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, array.shape[0]-1, array.shape[0])
    j_linespace = np.linspace(0, array.shape[1]-1, array.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within inner ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_inner = (j_trans**2 / semi_maj_inner**2) + (i_trans**2 / semi_min_inner**2 )

    # Use meshgrids to create array identifying which coordinates lie within outer ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_outer = (j_trans**2 / semi_maj_outer**2) + (i_trans**2 / semi_min_outer**2 )

    # Calculate flux & pixels in aperture, and store pixel values
    annulus_where = np.where( (ellipse_check_outer<=1) & (ellipse_check_inner>1) & (np.isnan(array)==False) )
    annulus_tot = sum( array[ annulus_where ] )
    annulus_count = annulus_where[0].shape[0]
    annulus_pix = array[ annulus_where ]
    annulus_nan = np.where( (ellipse_check_outer<=1) & (ellipse_check_inner>1) & (np.isnan(array)==True) )

    # Return results
    return [annulus_tot, annulus_count, annulus_pix, annulus_nan]

# -----------------------------------------------------------------

def LogError(value, error):

    """
    New function to convert an uncertainty to log space
    Input: Value, uncertainty
    Output: Logarithmic uncertainty
    """

    value, error = np.array(value), np.array(error)
    frac = 1.0 + (error/value)
    error_up = value * frac
    error_down = value / frac
    log_error_up = np.abs( np.log10(error_up) - np.log10(value) )
    log_error_down = np.abs( np.log10(value) - np.log10(error_down) )
    return 0.5*(log_error_up+log_error_down)

# -----------------------------------------------------------------

def EllipseSumUpscale(cutout, rad, axial_ratio, angle, i_centre, j_centre, upscale=1):

    """
    Function to sum all elements in an ellipse centred on the middle of an array that has been resized to allow better pixel sampling
    Input: Array, semi-major axis (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse, upscaling factor
    Returns: Numpy array containing the sum of the pixel values in the ellipse, the total number of pixels counted, and an array containing the pixel values
    """

    # Create slice of input array, containing only the region of interest
    i_cutout_min = np.floor(max([0, i_centre-rad]))
    i_cutout_max = np.ceil(min([(cutout.shape)[0], i_centre+rad]))
    j_cutout_min = np.floor(max([0, j_centre-rad]))
    j_cutout_max = np.ceil(min([(cutout.shape)[1], j_centre+rad]))
    cutout_slice = cutout[ int(round(i_cutout_min)):int(round(i_cutout_max))+1, int(round(j_cutout_min)):int(round(j_cutout_max))+1 ]
    i_centre_slice = i_centre - i_cutout_min
    j_centre_slice = j_centre - j_cutout_min
    if cutout[i_centre,j_centre]!=cutout[i_centre_slice,j_centre_slice]:
        if np.isnan(cutout[i_centre,j_centre]==False) and np.isnan(cutout_slice[i_centre_slice,j_centre_slice]==False):
            print('SEVERE ERROR: EllipseSumUpscale check failed.')
            pdb.set_trace()
    else:
        cutout = cutout_slice
        i_centre = i_centre_slice
        j_centre = j_centre_slice

    # Resize array to increase pixel sampling, cupdate centre coords, and downscale pixel values accordinly to preserve flux
    cutout_inviolate = np.copy(cutout)
    cutout = np.zeros([cutout_inviolate.shape[0]*upscale, cutout_inviolate.shape[1]*upscale])
    scipy.ndimage.zoom(cutout_inviolate, upscale, output=cutout, order=0)
    cutout *= float(upscale)**-2.0
    i_centre = float(i_centre) * float(upscale)
    j_centre = float(j_centre) * float(upscale)

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj = float(rad) * float(upscale)
    semi_min = semi_maj / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, cutout.shape[0]-1, cutout.shape[0])
    j_linespace = np.linspace(0, cutout.shape[1]-1, cutout.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check = (j_trans**2 / semi_maj**2) + (i_trans**2 / semi_min**2 )

    # Calculate flux & pixels in aperture, and store pixel values
    ellipse_where = np.where( (ellipse_check<=1) & (np.isnan(cutout)==False) )
    ellipse_tot = sum( cutout[ ellipse_where ] )
    ellipse_count = ellipse_where[0].shape[0]
    ellipse_pix = cutout[ ellipse_where ]

    # Scale output values down to what they would've been for original array
    ellipse_count *= float(upscale)**-2.0

    # Return results
    return [ellipse_tot, ellipse_count, ellipse_pix]

# -----------------------------------------------------------------

def AnnulusSumUpscale(cutout, rad_inner, width, axial_ratio, angle, i_centre, j_centre, upscale=1):

    """
    Function to sum all elements in an annulus centred upon the middle of an array that has been resized to allow better pixel sampling
    Input: Array, semi-major axis of inside edge of annulus (pix), width of annulus (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse, upscaling factor
    Returns: Numpy array containing the sum of the pixel values in the annulus, the total number of pixels counted, and an array containing the pixel values
    """

    # Create slice of input array, containing only the region of interest
    i_cutout_min = np.floor(max([0, i_centre-(rad_inner+width)]))
    i_cutout_max = np.ceil(min([(cutout.shape)[0], i_centre+(rad_inner+width)]))
    j_cutout_min = np.floor(max([0, j_centre-(rad_inner+width)]))
    j_cutout_max = np.ceil(min([(cutout.shape)[1], j_centre+(rad_inner+width)]))
    cutout_slice = cutout[ int(round(i_cutout_min)):int(round(i_cutout_max))+1, int(round(j_cutout_min)):int(round(j_cutout_max))+1 ]
    i_centre_slice = i_centre - i_cutout_min
    j_centre_slice = j_centre - j_cutout_min
    if cutout[i_centre,j_centre]!=cutout_slice[i_centre_slice,j_centre_slice]:
        if np.isnan(cutout[i_centre,j_centre]==False) and np.isnan(cutout_slice[i_centre_slice,j_centre_slice]==False):
            print('SEVERE ERROR: AnnulusQuickSum check failed.')
            pdb.set_trace()
    else:
        cutout = cutout_slice
        i_centre = i_centre_slice
        j_centre = j_centre_slice

    # Resize array to increase pixel sampling, update centre coords, and downscale pixel values accordinly to preserve flux
    cutout_inviolate = np.copy(cutout)
    cutout = np.zeros([cutout_inviolate.shape[0]*upscale, cutout_inviolate.shape[1]*upscale])
    scipy.ndimage.zoom(cutout_inviolate, upscale, output=cutout, order=0)
    cutout *= float(upscale)**-2.0
    i_centre = float(i_centre) * float(upscale)
    j_centre = float(j_centre) * float(upscale)

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj_inner = float(rad_inner) * float(upscale)
    semi_min_inner = semi_maj_inner / float(axial_ratio)
    semi_maj_outer = ( float(rad_inner) * float(upscale) ) + ( float(width) * float(upscale) )
    semi_min_outer  = semi_maj_outer / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, cutout.shape[0]-1, cutout.shape[0])
    j_linespace = np.linspace(0, cutout.shape[1]-1, cutout.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create cutout identifying which coordinates lie within inner ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_inner = (j_trans**2 / semi_maj_inner**2) + (i_trans**2 / semi_min_inner**2 )

    # Use meshgrids to create cutout identifying which coordinates lie within outer ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_outer = (j_trans**2 / semi_maj_outer**2) + (i_trans**2 / semi_min_outer**2 )

    # Calculate flux & pixels in aperture, and store pixel values
    annulus_where = np.where( (ellipse_check_outer<=1) & (ellipse_check_inner>1) & (np.isnan(cutout)==False) )
    annulus_tot = sum( cutout[ annulus_where ] )
    annulus_count = annulus_where[0].shape[0]
    annulus_pix = cutout[ annulus_where ]

    # Scale output values down to what they would've been for original array
    annulus_count *= float(upscale)**-2.0

    # Return results
    return [annulus_tot, annulus_count, annulus_pix]

# -----------------------------------------------------------------
