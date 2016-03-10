#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.box Contains the Box class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np
from scipy import ndimage
from skimage.restoration.inpaint import inpaint_biharmonic

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant AstroMagic classes and modules
from ..basics import Position, Region, Rectangle, Extent
from ..tools import cropping, fitting, interpolation, plotting

# -----------------------------------------------------------------

class Box(np.ndarray):

    """
    This class ...
    """

    def __new__(cls, data, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param cls:
        :param input_array:
        :param info:
        :return:
        """

        # Create an object of the Box class
        obj = np.asarray(data).view(cls)

        # Set attributes of the object
        obj.x_min = x_min
        obj.x_max = x_max
        obj.y_min = y_min
        obj.y_max = y_max

        # Return the object
        return obj

    # -----------------------------------------------------------------

    @classmethod
    def cutout(cls, frame, center, x_radius, y_radius):

        """
        This class method ...
        :param frame:
        :param center:
        :param x_radius:
        :param y_radius:
        :return:
        """

        # Crop the frame
        cropped, x_min, x_max, y_min, y_max = cropping.crop(frame, center.x, center.y, x_radius, y_radius)

        # Check that the center position lies within the box
        assert (x_min <= center.x < x_max and y_min <= center.y < y_max)

        # Return a new box
        return cls(cropped, x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    @classmethod
    def cutout_limits(cls, frame, x_min, x_max, y_min, y_max, absolute=False):

        """
        This function ...
        :return:
        """

        # Crop the frame
        if absolute: cropped = cropping.crop_absolute(frame, x_min, x_max, y_min, y_max)
        else: cropped, x_min, x_max, y_min, y_max = cropping.crop_direct(frame, x_min, x_max, y_min, y_max)

        # Return a new box
        return cls(cropped, x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    @classmethod
    def from_rectangle(cls, frame, rectangle, absolute=False):

        """
        This function ...
        :param frame:
        :param rectangle:
        :return:
        """

        # TODO: fix this for rotated rectangles

        # Convert into integers
        x_min = int(round(rectangle.x_min))
        x_max = int(round(rectangle.x_max))
        y_min = int(round(rectangle.y_min))
        y_max = int(round(rectangle.y_max))

        if x_min == x_max:
            x_min = int(math.floor(rectangle.x_min))
            x_max = int(math.ceil(rectangle.x_max))
        if y_min == y_max:
            y_min = int(math.floor(rectangle.y_min))
            y_max = int(math.ceil(rectangle.y_max))

        if x_max - x_min < 2: x_max = x_min + 2
        if y_max - y_min < 2: y_max = y_min + 2

        # Create cutout
        return cls.cutout_limits(frame, x_min, x_max, y_min, y_max, absolute)

    # -----------------------------------------------------------------

    @classmethod
    def from_ellipse(cls, frame, ellipse, shape=None):

        """
        This function ...
        :param frame:
        :param ellipse:
        :param shape:
        :return:
        """

        # Get bouding box
        rectangle = ellipse.bounding_box

        # ...
        if shape is not None: rectangle = Rectangle(rectangle.center, Extent(0.5 * shape.x, 0.5 * shape.y))

        return cls.from_rectangle(frame, rectangle, absolute=(shape is not None))

    # -----------------------------------------------------------------

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        rel_y_min = box.y_min - self.y_min
        rel_y_max = box.y_max - self.y_min

        rel_x_min = box.x_min - self.x_min
        rel_x_max = box.x_max - self.x_min

        data = self[rel_y_min:rel_y_max, rel_x_min:rel_x_max]

        # Create the new box and return it
        return Box(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # -----------------------------------------------------------------

    def zoom(self, center, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        # Calculate the size of the smaller box
        new_xsize = int(round(0.5 * self.xsize / factor))
        new_ysize = int(round(0.5 * self.ysize / factor))

        # Calculate the relative coordinate of the center
        rel_center = self.rel_position(center)

        # Create a smaller box
        data, rel_x_min, rel_x_max, rel_y_min, rel_y_max = cropping.crop(self, rel_center.x, rel_center.y, new_xsize, new_ysize)

        # Create the new box
        return Box(data, rel_x_min+self.x_min, rel_x_max+self.x_min, rel_y_min+self.y_min, rel_y_max+self.y_min)

    # -----------------------------------------------------------------

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        # see InfoArray.__array_finalize__ for comments
        if obj is None: return

        self.x_min = getattr(obj, 'x_min', None)
        self.x_max = getattr(obj, 'x_max', None)
        self.y_min = getattr(obj, 'y_min', None)
        self.y_max = getattr(obj, 'y_max', None)

    # -----------------------------------------------------------------

    def plot(self, frame=None):

        """
        This function ...
        :param frame:
        :return:
        """

        plotting.plot_box(self)

        # If frame is not None, plot 'zoom' plot

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        return self.shape[1]

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        return self.shape[0]

    # -----------------------------------------------------------------

    @property
    def x_slice(self):

        return slice(self.x_min, self.x_max)

    # -----------------------------------------------------------------

    @property
    def y_slice(self):

        return slice(self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def rel_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the relative position
        return Position(position.x - self.x_min, position.y - self.y_min)

    # -----------------------------------------------------------------

    def abs_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the absolute position
        return Position(position.x + self.x_min, position.y + self.y_min)

    # -----------------------------------------------------------------

    def fit_polynomial(self, order, mask=None):

        """
        This function ...
        :return:
        """

        # Do the fitting
        polynomial = fitting.fit_polynomial(self, order, mask=mask)

        # Evaluate the polynomial
        data = fitting.evaluate_model(polynomial, 0, self.xsize, 0, self.ysize)

        # Return a new box
        return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def gaussian_filter(self, sigma):

        """
        This function ...
        :return:
        """

        # Apply the filter
        data = ndimage.filters.gaussian_filter(self, sigma=sigma)

        # Return a new box
        return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def fit_model(self, center, model_name, sigma=None, max_center_offset=None, amplitude=None):

        """
        This function ...
        :return:
        """

        # Calculate the relative coordinate of the center
        rel_center = self.rel_position(center)

        # Fit a 2D Gaussian to the data
        if model_name == "Gaussian":

            # Do the fitting
            model = fitting.fit_2D_Gaussian(self, rel_center, sigma=sigma, max_center_offset=max_center_offset, amplitude=amplitude)

        # Fit an Airy Disk model to the data
        elif model_name == "Airy":

            # See https://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile and
            # http://astropy.readthedocs.org/en/latest/api/astropy.modeling.functional_models.AiryDisk2D.html#astropy.modeling.functional_models.AiryDisk2D
            #sigma = 0.42 * airy.radius * 0.81989397882
            radius = sigma / (0.42 * 0.81989397882) if sigma is not None else None

            # Do the fitting
            model = fitting.fit_2D_Airy(self, rel_center, radius=radius, max_center_offset=max_center_offset, amplitude=amplitude)

        # Fit a 2D (vertically) shifted Gaussian to the data
        elif model_name == "ShiftedGaussian":

            # Do the fitting
            model = fitting.fit_2D_ShiftedGaussian(self, rel_center, sigma=sigma, max_center_offset=max_center_offset, amplitude=amplitude)

        # Unknown model name
        else: raise ValueError("Model name should be 'Gaussian' or 'Airy'")

        # Shift the position of the model so that it represents its absolute position in the frame
        fitting.shift_model(model, self.x_min, self.y_min)

        # Return the model
        return model

    # -----------------------------------------------------------------

    def evaluate_model(self, model):

        """
        This function ...
        :param model:
        :return:
        """

        # Make a local copy of the model so that we can adapt its position to be relative to this box
        rel_model = fitting.shifted_model(model, -self.x_min, -self.y_min)

        # Create x and y meshgrid for evaluating
        y_values, x_values = np.mgrid[:self.ysize, :self.xsize]

        # Evaluate the model
        data = rel_model(x_values, y_values)

        # Return a new box
        return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

    # -----------------------------------------------------------------

    def interpolated(self, mask, method):

        """
        This function ...
        :param mask:
        :param method:
        :return:
        """

        # Fit a polynomial to the data
        if method == "polynomial":
            try:
                return self.fit_polynomial(3, mask=mask)
            except TypeError:
                from ..tools import plotting
                print("cutout=", type(self), mask.shape)
                print("mask=", type(mask), mask.shape)
                #plotting.plot_box(self)
                #plotting.plot_box(mask)
                #plotting.plot_box(np.ma.masked_array(self, mask=mask))
                mask = mask.eroded(2, 1)
                plotting.plot_box(np.ma.masked_array(self, mask=mask))
                return self.fit_polynomial(3, mask=mask)

        # Interpolate using the local mean method
        elif method == "local_mean":

            try:
                # Calculate the interpolated data
                data = interpolation.in_paint(self, mask)
            except IndexError:
                print(self)
                print(mask)

            # Create and return a new box
            return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Interpolate using inverse distance weighing
        elif method == "idw":

            # Calculate the interpolated data
            data = interpolation.in_paint(self, mask, method="idw")

            # Create and return a new box
            return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Calculate the mean value of the data
        elif method == "mean":

            mean = np.ma.mean(np.ma.masked_array(self, mask=mask))
            return self.full(mean)

        # Calculate the median value of the data
        elif method == "median":

            median = np.ma.median(np.ma.masked_array(self, mask=mask))
            return self.full(median)

        elif method == "biharmonic":

            data = inpaint_biharmonic(self, mask, multichannel=False)
            return Box(data, self.x_min, self.x_max, self.y_min, self.y_max)

        # Invalid option
        else: raise ValueError("Unknown interpolation method")

    # -----------------------------------------------------------------

    def value(self, position):

        """
        This function ...
        :param position:
        :param interpolate:
        :return:
        """

        # Get the x and y coordinate of the corresponding pixel
        x = int(round(position.x - self.x_min))
        y = int(round(position.y - self.y_min))

        # Return the pixel value
        return self[y, x]

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Convert to relative position
        rel_position = self.rel_position(position)

        # Calculate the pixel coordinates of the position
        x_pixel = int(round(rel_position.x))
        y_pixel = int(round(rel_position.y))

        # Check whether this box contains the position
        if x_pixel < 0 or y_pixel < 0 or x_pixel >= self.xsize or y_pixel >= self.ysize: return False
        else: return True

    # -----------------------------------------------------------------

    def replace(self, frame, where=None):

        """
        This function ...
        :param frame:
        :return:
        """

        # Replace the pixel values in the frame
        if where is None: frame[self.y_min:self.y_max, self.x_min:self.x_max] = self
        else: frame[self.y_min:self.y_max, self.x_min:self.x_max][where] = self[where]

        #from ..tools import plotting
        #plotting.plot_box(frame[self.y_min:self.y_max, self.x_min:self.x_max])
        #plotting.plot_box(np.ma.masked_array(frame[self.y_min:self.y_max, self.x_min:self.x_max], mask=where))

# -----------------------------------------------------------------
