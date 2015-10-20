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
import copy

# Import Astromagic modules
from ..core.masks import Mask

# Import astronomical modules
from astropy import log
import astropy.logger
import pyregion

# *****************************************************************

class ObjectExtractor(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # *****************************************************************

    def __init__(self):

        """
        The constructor ...
        """

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ###

        # Initialize an empty list for the sky objects
        self.objects = []

        # Set the frame to None
        self.frame = None

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the list of sky objects
        self.objects = []

        # Clear the frame
        self.frame = None

    # *****************************************************************

    def find_sources(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Looking for sources near the object positions")

        # Loop over all sky objects in the list
        for skyobject in self.objects:

            # If this sky object should be ignored, skip it
            if skyobject.ignore: continue

            # Find a source
            try:
                skyobject.find_source(self.frame, self.config.detection)

            except Exception as e:

                import traceback

                log.error("Error when finding source")
                print(type(e))
                print(e)
                traceback.print_exc()

                if self.config.plot_track_record_if_exception:

                    if skyobject.has_track_record: skyobject.track_record.plot()
                    else: print("Track record is not enabled")

                log.error("Continuing with next source")

        # Inform the user
        log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.objects), self.have_source/len(self.objects)*100.0))

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

    def set_special(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the region and create a mask from it
        region = pyregion.open(self.config.special_region).as_imagecoord(self.frame.wcs.to_header())
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all objects
        for skyobject in self.objects:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(self.frame.wcs)

            # Set special if position is covered by the mask
            if special_mask.masks(position): skyobject.special = True

    # *****************************************************************

    def set_ignore(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Load the region and create a mask from it
        region = pyregion.open(self.config.ignore_region).as_imagecoord(frame.wcs.to_header())
        ignore_mask = Mask(region.get_mask(shape=frame.shape))

        # Loop over all objects
        for skyobject in self.objects:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(frame.wcs)

            # Ignore if position is covered by the mask
            if ignore_mask.masks(position): skyobject.ignore = True

    # *****************************************************************

    @property
    def aperture_mask(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(self.frame))

        # Loop over all sky objects
        for skyobject in self.objects:

            # If the object does not have an aperture, skip it
            if not skyobject.has_aperture: continue

            # Create a mask from the aperture of the object
            object_mask_frame = Mask.from_aperture(self.frame.xsize, self.frame.ysize, skyobject.aperture, expansion_factor=self.config.aperture_mask.expansion_factor)

            # Now, we don't limit setting the mask within the source's cutout, because we expanded the apertures to perhaps a size larger than this cutout,
            # so just add the object_mask_frame to the total frame
            mask += object_mask_frame

        # Return the mask
        return mask

    # *****************************************************************

    def save_region(self):

        """
        This function ...
        """

        self.write_region(self.config.saving.region_path, self.config.saving.region_annotation)

    # *****************************************************************

    def save_masked_frame(self):

        """
        This function ...
        """

        # Create a frame where the objects are masked
        frame = copy.deepcopy(self.frame)
        frame[self.mask] = 0.0

        # Save the masked frame
        frame.save(self.config.saving.masked_frame_path)

    # *****************************************************************

    def save_result(self):

        """
        This function ...
        :return:
        """

        # Save the resulting frame
        self.frame.save(self.config.saving.result_path)

    # *****************************************************************

    @property
    def positions(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the galaxies
        for skyobject in self.objects:

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(skyobject.pixel_position(self.frame.wcs))

        # Return the list
        return positions

# *****************************************************************
