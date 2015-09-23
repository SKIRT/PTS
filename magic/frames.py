#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np

# *****************************************************************

class Frame(np.ndarray):

    """
    This class ...
    """

    # *****************************************************************

    def __new__(cls, input_array, coordinates=None, description=None):

        """
        This function ...
        :param cls:
        :param input_array:
        :param info:
        :return:
        """

        obj = np.asarray(input_array).view(cls)
        obj.coordinates = coordinates
        obj.description = description
        obj.selected = False

        return obj

    # *****************************************************************

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.coordinates = getattr(obj, 'coordinates', None)
        self.description = getattr(obj, 'description', None)

    # *****************************************************************

    def __init__(self, input_array, coordinates=None, description=None):

        self = np.asarray(input_array).view(self.__class__)
        self.coordinates = coordinates
        self.description = description
        self.selected = False

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
    def xsize(self): return self.shape[1]

    # *****************************************************************

    @property
    def ysize(self): return self.shape[0]

# *****************************************************************