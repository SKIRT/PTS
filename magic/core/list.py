#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.list Contains the CoordinateSystemList, FrameList, and ImageList classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.filter.filter import parse_filter, Filter
from ..region.list import SkyRegionList
from ...core.units.parsing import parse_unit as u
from .frame import Frame
from ..basics.coordinatesystem import CoordinateSystem
from ..tools import coordinates
from ..basics.coordinate import SkyCoordinate
from ..basics.stretch import SkyStretch
from ..region.rectangle import SkyRectangleRegion

# -----------------------------------------------------------------

class CoordinateSystemList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        THe constructor ...
        """

        # The coordinate systems
        self.systems = OrderedDict()

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return: 
        """

        return len(self.systems)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.systems.keys()

    # -----------------------------------------------------------------

    def __contains__(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if isinstance(fltr, basestring): fltr = parse_filter(fltr)
        return fltr in self.systems

    # -----------------------------------------------------------------

    def append(self, frame_or_wcs, fltr=None):

        """
        This function ...
        :param frame_or_wcs: 
        :param fltr: 
        :return: 
        """

        if isinstance(frame_or_wcs, Frame):
            wcs = frame_or_wcs.wcs
            fltr = frame_or_wcs.filter
        elif isinstance(frame_or_wcs, CoordinateSystem):
            wcs = frame_or_wcs
            if fltr is None: raise ValueError("Filter must be specified")
        else: raise ValueError("Invalid input")

        # Check the key
        if fltr in self.systems: raise ValueError("Already a coordinate system for the '" + str(fltr) + "' filter")

        # Add to the dictionary
        self.systems[fltr] = wcs

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        THis function ...
        :return: 
        """

        for name in self.systems: yield self.systems[name]

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Get the filter
        if isinstance(index_or_filter, basestring): fltr = parse_filter(index_or_filter)
        elif isinstance(index_or_filter, Filter): fltr = index_or_filter
        else: fltr = self.systems.keys()[index_or_filter]

        # Return the coordinate system
        return self.systems[fltr]

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def min_ra(self):

        """
        This function ...
        :return:
        """

        min_ra = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            wcs_min_ra = wcs.min_ra

            if min_ra is None or min_ra > wcs_min_ra: min_ra = wcs_min_ra

        return min_ra

    # -----------------------------------------------------------------

    @property
    def min_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_ra.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def max_ra(self):

        """
        This function ...
        :return:
        """

        max_ra = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            wcs_max_ra = wcs.max_ra

            if max_ra is None or max_ra < wcs_max_ra: max_ra = wcs_max_ra

        return max_ra

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def max_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_ra.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def ra_center(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.min_ra + self.max_ra)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def ra_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.ra_center.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def ra_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]

            if the_range is None: the_range = wcs.ra_range
            else: the_range.adjust(wcs.ra_range)

        # Return the range
        return the_range

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def min_dec(self):

        """
        This function ...
        :return:
        """

        min_dec = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            wcs_min_dec = wcs.min_dec

            if min_dec is None or min_dec > wcs_min_dec: min_dec = wcs_min_dec

        # Return
        return min_dec

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def min_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_dec.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def max_dec(self):

        """
        This function ...
        :return:
        """

        max_dec = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            wcs_max_dec = wcs.max_dec

            if max_dec is None or max_dec < wcs_max_dec: max_dec = wcs_max_dec

        # Return
        return max_dec

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def dec_center(self):

        """
        This function ...
        :return:
        """

        dec_center = 0.5 * (self.min_dec + self.max_dec)
        return dec_center

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def dec_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.dec_center.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def max_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_dec.to("deg").value

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def dec_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]

            if the_range is None: the_range = wcs.dec_range
            else: the_range.adjust(wcs.dec_range)

        # Return the dec range
        return the_range

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(self.ra_center_deg, self.dec_center_deg, unit="deg")

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(self.dec_center_deg, self.min_ra_deg, self.max_ra_deg))
        dec_distance = abs(self.max_dec_deg - self.min_dec_deg)

        # Create RA and DEC span as quantities
        ra_span = ra_distance * u("deg")
        dec_span = dec_distance * u("deg")

        # Return the center coordinate and the RA and DEC span
        return self.center, ra_span, dec_span

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range

        # Create box
        radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)
        box = SkyRectangleRegion(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for fltr in self.filters:

            wcs = self[fltr]
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

# -----------------------------------------------------------------

class FrameList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        THe constructor ...
        """

        # The frames
        self.frames = OrderedDict()

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return: 
        """

        return len(self.frames)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        THis function ...
        :return: 
        """

        return self.frames.keys()

    # -----------------------------------------------------------------

    def append(self, frame):

        """
        This function ...
        :param frame: 
        :return: 
        """

        # Check keys
        if frame.fltr in self.frames: raise ValueError("Already a frame for the '" + str(frame.filter) + "' filter")

        # Add the frame
        self.frames[frame.filter] = frame

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        THis function ...
        :return: 
        """

        for name in self.frames: yield self.frames[name]

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Get the filter
        if isinstance(index_or_filter, basestring): fltr = parse_filter(index_or_filter)
        elif isinstance(index_or_filter, Filter): fltr = index_or_filter
        else: fltr = self.frames.keys()[index_or_filter]

        # Return the frame
        return self.frames[fltr]

    # -----------------------------------------------------------------

    @property
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the frames
        for fltr in self.frames:
            if fwhm is None or self.frames[fltr].fwhm < fwhm: fwhm = self.frames[fltr].fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the frames
        for fltr in self.frames:
            if fwhm is None or self.frames[fltr].fwhm > fwhm: fwhm = self.frames[fltr].fwhm

        # Return the maximum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the frames
        for fltr in self.frames:

            wcs = self.frames[fltr].wcs
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the frames
        for fltr in self.frames:

            wcs = self.frames[fltr].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.frames: boxes_region.append(self.frames[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    @property
    def center_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.bounding_box.center

    # -----------------------------------------------------------------

    @property
    def coordinate_systems(self):

        """
        THis function ...
        :return: 
        """

        for fltr in self.frames: yield self.frames[fltr].wcs

    # -----------------------------------------------------------------

    @property
    def min_ra(self):

        """
        This function ...
        :return:
        """

        min_ra = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]
            wcs_min_ra = wcs.min_ra

            if min_ra is None or min_ra > wcs_min_ra: min_ra = wcs_min_ra

        return min_ra

    # -----------------------------------------------------------------

    @property
    def min_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_ra.to("deg").value

    # -----------------------------------------------------------------

    @property
    def max_ra(self):

        """
        This function ...
        :return:
        """

        max_ra = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]
            wcs_max_ra = wcs.max_ra

            if max_ra is None or max_ra < wcs_max_ra: max_ra = wcs_max_ra

        return max_ra

    # -----------------------------------------------------------------

    @property
    def max_ra_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_ra.to("deg").value

    # -----------------------------------------------------------------

    @property
    def ra_center(self):

        """
        This function ...
        :return:
        """

        return 0.5 * (self.min_ra + self.max_ra)

    # -----------------------------------------------------------------

    @property
    def ra_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.ra_center.to("deg").value

    # -----------------------------------------------------------------

    @property
    def ra_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]

            if the_range is None: the_range = wcs.ra_range
            else: the_range.adjust(wcs.ra_range)

        # Return the range
        return the_range

    # -----------------------------------------------------------------

    @property
    def min_dec(self):

        """
        This function ...
        :return:
        """

        min_dec = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]
            wcs_min_dec = wcs.min_dec

            if min_dec is None or min_dec > wcs_min_dec: min_dec = wcs_min_dec

        # Return
        return min_dec

    # -----------------------------------------------------------------

    @property
    def min_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.min_dec.to("deg").value

    # -----------------------------------------------------------------

    @property
    def max_dec(self):

        """
        This function ...
        :return:
        """

        max_dec = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]
            wcs_max_dec = wcs.max_dec

            if max_dec is None or max_dec < wcs_max_dec: max_dec = wcs_max_dec

        # Return
        return max_dec

    # -----------------------------------------------------------------

    @property
    def dec_center(self):

        """
        This function ...
        :return:
        """

        dec_center = 0.5 * (self.min_dec + self.max_dec)
        return dec_center

    # -----------------------------------------------------------------

    @property
    def dec_center_deg(self):

        """
        This function ...
        :return:
        """

        return self.dec_center.to("deg").value

    # -----------------------------------------------------------------

    @property
    def max_dec_deg(self):

        """
        This function ...
        :return:
        """

        return self.max_dec.to("deg").value

    # -----------------------------------------------------------------

    @property
    def dec_range(self):

        """
        This function ...
        :return:
        """

        the_range = None

        # Loop over the coordinate systems
        #for fltr in self.coordinate_systems:
        for wcs in self.coordinate_systems:

            #wcs = self.coordinate_systems[fltr]

            if the_range is None: the_range = wcs.dec_range
            else: the_range.adjust(wcs.dec_range)

        # Return the dec range
        return the_range

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(self.ra_center_deg, self.dec_center_deg, unit="deg")

    # -----------------------------------------------------------------

    @property
    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(self.dec_center_deg, self.min_ra_deg, self.max_ra_deg))
        dec_distance = abs(self.max_dec_deg - self.min_dec_deg)

        # Create RA and DEC span as quantities
        ra_span = ra_distance * u("deg")
        dec_span = dec_distance * u("deg")

        # Return the center coordinate and the RA and DEC span
        return self.center, ra_span, dec_span

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Get coordinate range
        center, ra_span, dec_span = self.coordinate_range

        # Create box
        radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)
        box = SkyRectangleRegion(center, radius)

        # Return the box
        return box

    # -----------------------------------------------------------------

    def converted_to_same_unit(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return: 
        """

        # Inform the user
        log.info("Converting frames to the same unit ...")

        # Check if the unit is defined
        if "unit" in kwargs: unit = kwargs.pop("unit")
        else: unit = self[0].unit

        # Debugging
        log.debug("Converting to unit '" + str(unit) + "' ...")

        # Initialize list for converted frames
        new_frames = FrameList()

        # Convert all
        for fltr in self.frames:
            frame = self.frames[fltr]
            new_frames.append(frame.converted_to(unit, **kwargs))

        # Return the new set of frames
        return new_frames

# -----------------------------------------------------------------

class NamedFrameList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        self.frames = OrderedDict()

    # -----------------------------------------------------------------

    @classmethod
    def from_dictionary(cls, dictionary):

        """
        This function ...
        :param dictionary: 
        :return: 
        """

        new = cls()
        for name in dictionary: new.append(dictionary[name], name=name)
        return new

    # -----------------------------------------------------------------

    def append(self, frame, name=None):

        """
        This function ...
        :param frame: 
        :param name: 
        :return: 
        """

        if name is None: name = frame.name
        if name is None: raise ValueError("Frame does not have a name")

        if name in self.frames: raise ValueError("Already a frame with the name '" + name + "'")

        self.frames[name] = frame

# -----------------------------------------------------------------

class ImageList(object):
        
    """
    This class ...
    """

    def __init__(self):

        """
        THe constructor ...
        """

        self.images = OrderedDict()

# -----------------------------------------------------------------

class NamedImageList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.images = OrderedDict()

    # -----------------------------------------------------------------

    @classmethod
    def from_dictionary(cls, dictionary):

        """
        This function ...
        :param dictionary: 
        :return: 
        """

        new = cls()
        for name in dictionary: new.append(dictionary[name], name=name)
        return new

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        THis function ...
        :return: 
        """

        return self.images.keys()

    # -----------------------------------------------------------------

    @property
    def items(self):

        """
        This function ...
        :return: 
        """

        return self.images.items()

    # -----------------------------------------------------------------

    def append(self, image, name=None):

        """
        This function ...
        :param image: 
        :param name: 
        :return: 
        """

        if name is None: name = image.name
        if name is None: raise ValueError("Image does not have a name")

        if name in self.images: raise ValueError("Already an image with the name '" + name + "'")

        self.images[name] = image

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return: 
        """

        return len(self.images)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        THis function ...
        :return: 
        """

        for name in self.images: yield self.images[name]

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_name):

        """
        This function ...
        :param index_or_name:
        :return: 
        """

        # Get the filter
        if isinstance(index_or_name, basestring): name = index_or_name
        else: name = self.images.keys()[index_or_name]

        # Return the image
        return self.images[name]

    # -----------------------------------------------------------------

    def __contains__(self, name):

        """
        This function ...
        :param name: 
        :return: 
        """

        return name in self.images

    # -----------------------------------------------------------------

    @property
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.names:
            if fwhm is None or self.images[name].fwhm < fwhm: fwhm = self.images[name].fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.names:
            if fwhm is None or self.images[name].fwhm > fwhm: fwhm = self.images[name].fwhm

        # Return the maximum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.names:

            wcs = self.images[name].wcs
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.names:

            wcs = self.images[name].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.names: boxes_region.append(self.images[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    @property
    def center_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.bounding_box.center

# -----------------------------------------------------------------