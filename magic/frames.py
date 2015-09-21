#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np

# *****************************************************************

class Frame(object):

    """
    This class ...
    """

    # *****************************************************************

    def __init__(self, data, coordinates, description):

        """
        The constructor ...
        """

        # Copy the data
        self.data = data

        # Copy the coordinate system
        self.coordinates = coordinates

        # Set as unactive initially
        self.selected = False

        # Set the description
        self.description = description

    # *****************************************************************

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # *****************************************************************

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # *****************************************************************

    @property
    def xsize(self): return self.data.shape[1]

    # *****************************************************************

    @property
    def ysize(self): return self.data.shape[0]

    # *****************************************************************

    @property
    def dtype(self): return self.data.dtype.name

    # *****************************************************************

    @property
    def mean(self): return np.mean(self.data)

    # *****************************************************************

    @property
    def median(self): return np.median(self.data)

    # *****************************************************************

    @property
    def min(self): return np.min(self.data)

    # *****************************************************************

    @property
    def max(self): return np.max(self.data)

    # *****************************************************************

    @property
    def stddev(self):

        # Set the delta degrees of freedom
        ddof = 1

        # Return the standard deviation of the data
        return np.std(self.data, ddof=ddof)
        
    @property
    def sum(self):
        
        return np.sum(self.data)

# *****************************************************************