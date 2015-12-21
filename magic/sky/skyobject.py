#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.skyobject Contains the abstract SkyObject class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta
from abc import abstractmethod

# Import astronomical modules
from astropy.wcs.wcs import NoConvergence
from photutils import segment_properties, properties_table
from photutils import EllipticalAperture

# Import the relevant AstroMagic classes and modules
from ..basics import Position, TrackRecord
from ..analysis import sources

# -----------------------------------------------------------------

class SkyObject(object):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, position):

        """
        The constructor ...
        """

        # Position
        self.position = position

        # Set the source attribute to None initially
        self.source = None

        # Set the aperture (for saturation) to None initially
        self.aperture = None

        # Initialize a track record of sources
        self.track_record = None

        # Set the special flag to False initially
        self.special = False

        # Set the ignore flag to False initially
        self.ignore = False

    # -----------------------------------------------------------------

    @property
    def has_source(self):

        """
        This function ...
        :return:
        """

        return self.source is not None

    # -----------------------------------------------------------------

    @property
    def has_aperture(self):

        """
        This function ...
        :return:
        """

        return self.aperture is not None

    # -----------------------------------------------------------------

    @property
    def has_track_record(self):

        """
        This function ...
        :return:
        """

        return self.track_record is not None

    # -----------------------------------------------------------------

    def enable_track_record(self):

        """
        This function ...
        :return:
        """

        # Create a new track record
        self.track_record = TrackRecord("Start")

    # -----------------------------------------------------------------

    @abstractmethod
    def ellipse_parameters(self, wcs, pixelscale, initial_radius):

        """
        This function ...
        :param wcs:
        :param pixelscale:
        :param initial_radius:
        :return:
        """

        return

    # -----------------------------------------------------------------

    def pixel_position(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        #print(self.position.ra.value, self.position.dec.value)
        #print(self.position.ra.hms)
        #print(self.position.to_string('hmsdms'))

        # Get the x and y coordinate of the object's position
        #try:
        #    x, y = self.position.to_pixel(wcs, origin=0)
        #except NoConvergence:
        #    x, y = self.position.to_pixel(wcs, origin=0, mode='wcs')  # Ignore distortions

        x, y = self.position.to_pixel(wcs, origin=0, mode='wcs')

        # Return the position in pixel coordinates
        return Position(x, y)

    # -----------------------------------------------------------------

    def find_source(self, frame, config):

        """
        This function ...
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("detection")

        # Get the parameters of the circle
        center, radius, angle = self.ellipse_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

        # Find a source
        self.source = sources.find_source(frame, center, radius, angle, config, self.track_record, special=self.special)

    # -----------------------------------------------------------------

    def find_aperture(self, frame, config):

        """
        This function ...
        :return:
        """

        #print(self.source.cutout.xsize)
        #print(self.source.cutout.ysize)

        #from ..tools import plotting
        #plotting.plot_box(self.source.cutout)

        props = segment_properties(self.source.cutout, self.source.mask)
        #tbl = properties_table(props)

        x_shift = self.source.cutout.x_min
        y_shift = self.source.cutout.y_min

        # Since there is only one segment in the self.source.mask (the center segment), the props
        # list contains only one entry (one galaxy)
        properties = props[0]

        # Obtain the position, orientation and extent
        position = (properties.xcentroid.value + x_shift, properties.ycentroid.value + y_shift)
        a = properties.semimajor_axis_sigma.value * config.sigma_level
        b = properties.semiminor_axis_sigma.value * config.sigma_level
        theta = properties.orientation.value

        # Calculate the difference (in number of pixels) between the aperture center and the position of the sky object
        difference = self.pixel_position(frame.wcs) - Position(position[0], position[1])

        # Set the aperture if the offset is smaller than or equal to the specified maximum
        if difference.norm <= config.max_offset or config.max_offset is None:

            # Create the aperture
            self.aperture = EllipticalAperture(position, a, b, theta=theta)

# -----------------------------------------------------------------
