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

    def __new__(cls, input_array, wcs=None, pixelscale=None, description=None, selected=False):

        """
        This function ...
        :param cls:
        :param input_array:
        :param info:
        :return:
        """

        obj = np.asarray(input_array).view(cls)
        obj.wcs = wcs
        obj.pixelscale = pixelscale
        obj.description = description
        obj.selected = selected

        return obj

    # *****************************************************************

    def __array_finalize__(self, obj):

        """
        This function ...
        :param obj:
        :return:
        """

        if obj is None: return
        self.wcs = getattr(obj, 'wcs', None)
        self.description = getattr(obj, 'description', None)
        self.selected = getattr(obj, 'selected', False)

    # *****************************************************************

    #def __init__(self, input_array, coordinates=None, description=None, selected=False):

        #self = np.asarray(input_array).view(self.__class__)
        #self.coordinates = coordinates
        #self.description = description
        #self.selected = selected

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