#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import Astromagic modules
from ..tools import cropping
from ..tools import fitting
from .vector import Position, Extent

# *****************************************************************

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

    # *****************************************************************

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
        try:
            assert (x_min <= center.x < x_max and y_min <= center.y < y_max)
        except AssertionError:
            print(x_min, x_max, center.x)
            print(y_min, y_max, center.y)

        # Return a new box
        return cls(cropped, x_min, x_max, y_min, y_max)

    # *****************************************************************

    @classmethod
    def from_ellipse(cls, frame, center, radius, angle):

        """
        This function ...
        :param frame:
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        # Set the x and y radius
        if isinstance(radius, Extent):

            x_radius = np.abs(np.sin(angle)*radius.x + np.cos(angle)*radius.y)
            y_radius = np.abs(np.cos(angle)*radius.x - np.sin(angle)*radius.y)

        else:

            x_radius = radius
            y_radius = radius

        # Return a new box
        return cls.cutout(frame, center, x_radius, y_radius)

    # *****************************************************************

    def box_like(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        data = self[box.y_min - self.y_min:box.y_max - self.y_min,
                    box.x_min - self.x_min:box.x_max - self.x_min]

        # Create the new box
        return Box(data, box.x_min, box.x_max, box.y_min, box.y_max)

    # *****************************************************************

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

    # *****************************************************************

    def plot(self, frame=None):

        """
        This function ...
        :param frame:
        :return:
        """

        pass

        # If frame is not None, plot 'zoom' plot

    # *****************************************************************

    @property
    def xsize(self):

        return self.shape[1]

    # *****************************************************************

    @property
    def ysize(self):

        return self.shape[0]

    # *****************************************************************

    def rel_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the relative position
        return Position(x=position.x-self.x_min, y=position.y-self.y_min)

    # *****************************************************************

    def abs_position(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        # Return the absolute position
        return Position(x=position.x+self.x_min, y=position.y+self.y_min)

    # *****************************************************************

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

# *****************************************************************
