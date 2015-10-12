#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
from abc import ABCMeta
from abc import abstractmethod
import numpy as np

# Import Astromagic modules
from ..core.masks import Mask

# Import astronomical modules
from astropy import log
import pyregion

# *****************************************************************

class ObjectExtractor(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # *****************************************************************

    def __init__(self, config):

        """
        The constructor ...
        """

        # Set attributes
        self.config = config

        # Set the log level
        log.setLevel(self.config.log_level)

        # Initialize an empty list for the sky objects
        self.objects = []

        # Set special mask
        #self.special_mask = None

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the list of sky objects
        self.objects = []

    # *****************************************************************

    def find_sources(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Looking for sources near the object positions")

        # Loop over all sky objects in the list
        for skyobject in self.objects:

            # Find a source
            skyobject.find_source(frame, self.config.detection)

        # Inform the user
        log.debug("Success ratio: {0:.2f}%".format(self.have_source/len(self.objects)*100.0))

    # *****************************************************************

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for skyobject in self.objects: count += skyobject.has_source
        return count

    # *****************************************************************

    def set_special(self, frame):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the region and create a mask from it
        region = pyregion.open(self.config.special_region).as_imagecoord(frame.wcs.to_header())
        special_mask = region.get_mask(shape=frame.shape)

        # Loop over all objects
        for skyobject in self.objects:

            x_center, y_center = skyobject.position.to_pixel(frame.wcs, origin=0)

            # Calculate x and y of the pixel corresponding to the object's position
            x = int(round(x_center))
            y = int(round(y_center))

            # Set special if position is covered by the mask
            if special_mask[y,x]: skyobject.special = True

    # *****************************************************************

    def create_mask(self, frame):

        """
        This function ...

        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(frame))

        # Loop over all sky objects
        for skyobject in self.objects:

            # If no source was found for the object, skip it
            if not skyobject.has_source: continue

            # Add this sky object to the mask
            if self.config.mask.use_aperture and skyobject.has_aperture:

                object_mask_frame = Mask.from_aperture(frame.xsize, frame.ysize, skyobject.aperture)
                object_mask = object_mask_frame[skyobject.source.cutout.y_min:skyobject.source.cutout.y_max, skyobject.source.cutout.x_min:skyobject.source.cutout.x_max]

            else: object_mask = skyobject.source.mask

            # Add this galaxy to the total mask
            mask[skyobject.source.cutout.y_min:skyobject.source.cutout.y_max, skyobject.source.cutout.x_min:skyobject.source.cutout.x_max] += object_mask

        # Return the mask
        return mask

# *****************************************************************
