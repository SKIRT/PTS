#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import Astromagic modules
from ..core.masks import Mask
from ..core.regions import ellipse_parameters
from ..core.source import Source
from ..core.vector import Position, Extent
from ..tools import logging

# Import astronomical modules
from astropy import log
import pyregion
import astropy.units as u
from astropy.coordinates import Angle

# -----------------------------------------------------------------

class ObjectExtractor(object):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self):

        """
        The constructor ...
        """

        ## Attributes

        # Initialize an empty list for the sky objects
        self.objects = []

        # Initialize an empty list to contain the manual sources
        self.manual_sources = []

        # Set the frame to None
        self.frame = None

        # Set the mask to None
        self.mask = None

    # -----------------------------------------------------------------

    def setup(self, frame):

        """
        This function ...
        """

        # Make a local reference to the passed frame
        self.frame = frame

        # Set-up the logging system
        log.setLevel(self.config.logging.level)  # the logging level
        if self.config.logging.path is not None: logging.link_file_log(log, self.config.logging.path, self.config.logging.level)

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask(np.zeros_like(self.frame))

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the extractor")

        # Clear the list of sky objects
        self.objects = []

        # Clear the list of manual sources
        self.manual_sources = []

        # Clear the frame and the mask
        self.frame = None
        self.mask = None

    # -----------------------------------------------------------------

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for skyobject in self.objects: count += skyobject.has_source
        return count

    # -----------------------------------------------------------------

    def set_special(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Setting special region from " + self.config.special_region)

        # Load the region and create a mask from it
        region = pyregion.open(self.config.special_region).as_imagecoord(self.frame.wcs.to_header())
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all objects
        for skyobject in self.objects:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(self.frame.wcs)

            # Set special if position is covered by the mask
            if special_mask.masks(position): skyobject.special = True

    # -----------------------------------------------------------------

    def set_ignore(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        log.info("Setting region to ignore for subtraction from " + self.config.ignore_region)

        # Load the region and create a mask from it
        region = pyregion.open(self.config.ignore_region).as_imagecoord(self.frame.wcs.to_header())
        ignore_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all objects
        for skyobject in self.objects:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(self.frame.wcs)

            # Ignore if position is covered by the mask
            if ignore_mask.masks(position): skyobject.ignore = True

    # -----------------------------------------------------------------

    def set_manual(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Setting region for manual star extraction from " + self.config.manual_region)

        # Load the region and create a mask from it
        region = pyregion.open(self.config.manual_region).as_imagecoord(self.frame.wcs.to_header())

        # Loop over the shapes in the region
        for shape in region:

            # Get the center and radius of the shape (can be a circle or an ellipse)
            x_center, y_center, x_radius, y_radius = ellipse_parameters(shape)

            # Create a source
            source = Source(self.frame, Position(x_center, y_center), Extent(x_radius, y_radius), Angle(0.0, u.deg), self.config.manual.background_outer_factor)

            # Add the source to the list of manual sources
            self.manual_sources.append(source)

    # -----------------------------------------------------------------

    def remove_manual(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing manually located stars from the frame")

        # Loop over each item in the list of manual sources
        for source in self.manual_sources:

            # Estimate the background for the source
            source.estimate_background(self.config.manual.interpolation_method, self.config.manual.sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing table to " + self.config.writing.table_path)

        # Write the table to file
        self.table.write(self.config.writing.table_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Writing masked frame to " + self.config.writing.masked_frame_path)

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = 0.0

        # Write out the masked frame
        frame.save(self.config.writing.masked_frame_path)

    # -----------------------------------------------------------------

    def write_result(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing resulting frame to " + self.config.writing.result_path)

        # Write out the resulting frame
        self.frame.save(self.config.writing.result_path)

    # -----------------------------------------------------------------

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

# -----------------------------------------------------------------
