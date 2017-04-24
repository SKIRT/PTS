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

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.filter.filter import parse_filter, Filter
from ..region.list import SkyRegionList
from ...core.units.parsing import parse_unit as u
from .frame import Frame
from .image import Image
from ..basics.coordinatesystem import CoordinateSystem
from ..tools import coordinates
from ..basics.coordinate import SkyCoordinate
from ..basics.stretch import SkyStretch
from ..region.rectangle import SkyRectangleRegion
from ...core.basics.containers import KeyList, NamedList
from ...core.tools import types
from ...core.tools import filesystem as fs
from ..convolution.aniano import AnianoKernels

# -----------------------------------------------------------------

class CoordinateSystemList(KeyList):

    """
    This class ...
    """

    def __init__(self):

        """
        THe constructor ...
        """

        # Call the constructor of the base class
        super(CoordinateSystemList, self).__init__()

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        new = cls()
        for path in fs.files_in_path(path): new.append(Frame.from_file(path))
        return new

    # -----------------------------------------------------------------

    @property
    def systems(self): # an alias for the contents for this subclass

        """
        This function ...
        :return: 
        """

        return self.contents

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.keys

    # -----------------------------------------------------------------

    def __contains__(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.keys

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
        if fltr in self.keys: raise ValueError("Already a coordinate system for the '" + str(fltr) + "' filter")

        # Call the append function of the base class
        super(CoordinateSystemList, self).append(fltr, wcs)

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Get the filter
        if types.is_string_type(index_or_filter): index_or_filter = parse_filter(index_or_filter)

        # Call the function of the base class
        return super(CoordinateSystemList, self).__getitem__(index_or_filter)

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

class FrameList(KeyList):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        THe constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FrameList, self).__init__()

        # Add frames
        for frame in args: self.append(frame)
        for filter_name in kwargs: self.append(kwargs[filter_name], fltr=parse_filter(filter_name))

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        new = cls()
        for path in fs.files_in_path(path, extension="fits"): new.append(Frame.from_file(path))
        return new

    # -----------------------------------------------------------------

    @property
    def frames(self): # an alias for the contents for this subclass

        """
        This function ...
        :return: 
        """

        return self.contents

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        THis function ...
        :return: 
        """

        return self.keys

    # -----------------------------------------------------------------

    def append(self, frame, fltr=None):

        """
        This function ...
        :param frame: 
        :param fltr:
        :return: 
        """

        if fltr is None: fltr = frame.filter

        # Check keys
        if fltr in self.frames: raise ValueError("Already a frame for the '" + str(fltr) + "' filter")

        # Call the function of the base class
        super(FrameList, self).append(fltr, frame)

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Get the filter
        if types.is_string_type(index_or_filter): index_or_filter = parse_filter(index_or_filter)

        # Call the function of the base class
        return super(FrameList, self).__getitem__(index_or_filter)

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

    def convolve_to_highest_fwhm(self):

        """
        This function ...
        :return: 
        """

        new_frames = convolve_to_highest_fwhm(*self.values)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def rebin_to_highest_pixelscale(self):

        """
        This function ...
        :return: 
        """

        new_frames = rebin_to_highest_pixelscale(*self.values)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convert_to_same_unit(self, unit):

        """
        This function ...
        :param unit:
        :return: 
        """
        
        new_frames = convert_to_same_unit(*self.values, unit=unit)
        self.remove_all()
        for frame in new_frames: self.append(frame)

# -----------------------------------------------------------------

class NamedFrameList(NamedList):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(NamedFrameList, self).__init__()

    # -----------------------------------------------------------------

    @property
    def frames(self): # an alias for the contents for this subclass

        """
        This function ...
        :return: 
        """

        return self.contents

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, contains=None):

        """
        This function ...
        :param path: 
        :param contains:
        :return: 
        """

        new = cls()
        for path, name in fs.files_in_path(path, returns=["path", "name"], extension="fits", contains=contains): new.append(name, Frame.from_file(path))
        return new

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

        # Call the append function of the base class
        super(NamedFrameList, self).append(name, frame)

# -----------------------------------------------------------------

class ImageList(KeyList):
        
    """
    This class ...
    """

    def __init__(self):

        """
        THe constructor ...
        """

        # Call the constructor of the base class
        super(ImageList, self).__init__()

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        new = cls()
        for path in fs.files_in_path(path): new.append(Image.from_file(path))
        return new

    # -----------------------------------------------------------------

    def append(self, image):

        """
        This function ...
        :param image: 
        :return: 
        """

        # Check keys
        if image.fltr in self.images: raise ValueError("Already an image for the '" + str(image.filter) + "' filter")

        # Call the function of the base class
        super(ImageList, self).append(image.filter, image)

    # -----------------------------------------------------------------

    @property
    def images(self): # an alias for the contents for this subclass

        """
        This function ...
        :return: 
        """

        return self.contents

# -----------------------------------------------------------------

class NamedImageList(NamedList):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(NamedImageList, self).__init__()

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path):

        """
        This function ...
        :param path: 
        :return: 
        """

        new = cls()
        for path, name in fs.files_in_path(path, returns=["path", "name"], extension="fits"): new.append(Image.from_file(path), name)
        return new

    # -----------------------------------------------------------------

    @property
    def images(self): # an alias for the contents for this subclass

        """
        This function ...
        :return: 
        """

        return self.contents

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

        return self.keys

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

        # Check
        if name in self.images: raise ValueError("Already an image with the name '" + name + "'")

        # Call the append function of the base class
        super(NamedImageList, self).append(name, image)

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

def convert_to_same_unit(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return: 
    """

    # Inform the user
    log.info("Converting frames to the same unit ...")

    # Check if the unit is defined
    if "unit" in kwargs: unit = kwargs.pop("unit")
    else: unit = frames[0].unit

    # Debugging
    log.debug("Converting to unit '" + str(unit) + "' ...")

    # Initialize list for converted frames
    new_frames = []

    # Convert all
    for frame in frames: new_frames.append(frame.converted_to(unit, **kwargs))

    # Return the new set of frames
    return new_frames

# -----------------------------------------------------------------

def rebin_to_highest_pixelscale(*frames):

    """
    This function ...
    :param frames: 
    :return: 
    """

    # Inform the user
    log.info("Rebinning frames to the coordinate system with the highest pixelscale ...")

    highest_pixelscale = None
    highest_pixelscale_wcs = None

    # Loop over the frames
    for frame in frames:

        wcs = frame.wcs
        if highest_pixelscale is None or wcs.average_pixelscale > highest_pixelscale:

            highest_pixelscale = wcs.average_pixelscale
            highest_pixelscale_wcs = wcs

    # Initialize list for rebinned frames
    new_frames = []

    # Rebin
    for frame in frames:

        if frame.wcs == highest_pixelscale_wcs: new_frames.append(frame.copy())
        else:
            if frame.unit.is_per_angular_area: rebinned = frame.rebinned(highest_pixelscale_wcs)
            else:
                rebinned = frame.converted_to_corresponding_angular_area_unit()
                rebinned.rebin(highest_pixelscale_wcs)
            new_frames.append(rebinned)

    # Return the rebinned frames
    return new_frames

# -----------------------------------------------------------------

def convolve_to_highest_fwhm(*frames):

    """
    This function ...
    :param frames: 
    :return: 
    """

    aniano = AnianoKernels()

    # Inform the user
    log.info("Convolving frames to the resolution of the frame with the highest FWHM ...")

    highest_fwhm = None
    highest_fwhm_filter = None

    # Loop over the frames
    for frame in frames:

        if highest_fwhm is None or frame.fwhm > highest_fwhm:

            highest_fwhm = frame.fwhm
            highest_fwhm_filter = frame.filter

    # Initialize list for convolved frames
    new_frames = []

    # Convolve
    for frame in frames:

        if frame.filter == highest_fwhm_filter: new_frames.append(frame.copy())
        else: new_frames.append(frame.convolved(aniano.get_kernel(frame.filter, highest_fwhm_filter)))

    # Return the convolved frames
    return new_frames

# -----------------------------------------------------------------
