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
import numpy as np

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.filter.filter import parse_filter
from ..region.list import SkyRegionList
from ...core.units.parsing import parse_unit as u
from .frame import Frame
from .image import Image
from ..basics.coordinatesystem import CoordinateSystem
from ..tools import coordinates
from ..basics.coordinate import SkyCoordinate
from ..basics.stretch import SkyStretch
from ..region.rectangle import SkyRectangleRegion
from ...core.basics.containers import NamedList, FilterBasedList
from ...core.tools import filesystem as fs
from ..convolution.aniano import AnianoKernels
from ..convolution.matching import MatchingKernels
from ..convolution.kernels import get_fwhm, get_average_variable_fwhm, has_variable_fwhm, has_average_variable_fwhm
from ...core.tools import sequences, types
from ...core.launch.pts import execute_pts_remote
from ...core.remote.remote import load_remote
from ...core.tools import introspection
from ...core.tools import time
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class CoordinateSystemList(FilterBasedList):

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

    def append(self, frame_or_wcs, fltr=None):

        """
        This function ...
        :param frame_or_wcs: 
        :param fltr: 
        :return: 
        """

        # Get WCS and filter
        if isinstance(frame_or_wcs, Frame):
            wcs = frame_or_wcs.wcs
            fltr = frame_or_wcs.filter
        elif isinstance(frame_or_wcs, CoordinateSystem):
            wcs = frame_or_wcs
            if fltr is None: raise ValueError("Filter must be specified")
        else: raise ValueError("Invalid input")

        # Check the key
        if fltr in self.filters: raise ValueError("Already a coordinate system for the '" + str(fltr) + "' filter")

        # Call the append function of the base class
        super(CoordinateSystemList, self).append(fltr, wcs)

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

class NamedCoordinateSystemList(NamedList):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs: 
        """

        # Call the constructor of the base class
        super(NamedCoordinateSystemList, self).__init__()

        # Add coordinate systems
        for name in kwargs: self.append(name, kwargs[name])

    # -----------------------------------------------------------------

    def append(self, frame_or_wcs, name=None):

        """
        This function ...
        :param frame_or_wcs:
        :param name: 
        :return: 
        """

        # Get WCS and name
        if isinstance(frame_or_wcs, Frame):
            wcs = frame_or_wcs.wcs
            if name is None: name = frame_or_wcs.name
        elif isinstance(frame_or_wcs, CoordinateSystem): wcs = frame_or_wcs
        else: raise ValueError("Invalid input")

        # Check whether name is defined
        if name is None: raise ValueError("Name not specified")

        # Call the append function of the base class
        super(NamedCoordinateSystemList, self).append(name, wcs)

    # -----------------------------------------------------------------

    @property
    def min_ra(self):

        """
        This function ...
        :return:
        """

        min_ra = None

        # Loop over the coordinate systems
        for name in self.names:

            wcs = self[name]
            wcs_min_ra = wcs.min_ra

            if min_ra is None or min_ra > wcs_min_ra: min_ra = wcs_min_ra

        # Return
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
        for name in self.names:

            wcs = self[name]
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
        for name in self.names:

            wcs = self[name]

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
        for name in self.names:

            wcs = self[name]
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
        for name in self.names:

            wcs = self[name]
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
        for name in self.names:

            wcs = self[name]

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

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for name in self.names:

            wcs = self[name]
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the coordinate systems
        for name in self.names:

            wcs = self[name]
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

# -----------------------------------------------------------------

class FrameList(FilterBasedList):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        THe constructor ...
        :param args:
        :param kwargs:
        """

        # LAZY?
        lazy = kwargs.pop("lazy", False)
        if lazy: lazy_evaluator = Frame.from_file
        else: lazy_evaluator = None

        # Call the constructor of the base class
        super(FrameList, self).__init__(lazy=lazy, lazy_evaluator=lazy_evaluator)

        # Add frames
        for frame in args: self.append(frame)
        for filter_name in kwargs: self.append(kwargs[filter_name], fltr=parse_filter(filter_name))

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, recursive=False, contains=None, not_contains=None):

        """
        This function ...
        :param path:
        :param recursive:
        :param contains:
        :param not_contains:
        :return: 
        """

        new = cls()
        for path in fs.files_in_path(path, extension="fits", recursive=recursive, contains=contains, not_contains=not_contains): new.append(Frame.from_file(path))
        return new

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, *paths, **kwargs):

        """
        This function ...
        :param paths:
        :param kwargs::
        :return:
        """

        # Initialize
        lazy = kwargs.pop("lazy", False)
        new = cls(lazy=lazy)

        if lazy and len(paths) > 0: raise ValueError("Cannot give positional arguments when lazy is True")

        # Add args
        if not lazy:
            for path in paths: new.append_from_file(path)

        # Add kwargs
        for filter_name in kwargs: new.append(kwargs[filter_name], fltr=filter_name)

        return new

    # -----------------------------------------------------------------

    def write_to_directory(self, path, update_path=False, replace=True):

        """
        This function ...
        :param path:
        :param update_path:
        :param replace:
        :return:
        """

        # Inform the user
        log.info("Writing all frames to the '" + path + "' directory ...")

        # Loop over the frames
        for fltr in self.filters:

            # Get the frame
            frame = self[fltr]

            # Determine the filepath
            filepath = fs.join(path, str(fltr) + ".fits")

            # Check
            if fs.is_file(filepath) and not replace: raise IOError("File '" + filepath + "' already exists")

            # Save the frame
            frame.saveto(filepath, update_path=update_path)

    # -----------------------------------------------------------------

    @property
    def filters(self):
        return self.keys

    # -----------------------------------------------------------------

    @property
    def frames(self): # an alias for the contents for this subclass
        return self.contents

    # -----------------------------------------------------------------

    def append(self, frame, fltr=None):

        """
        This function ...
        :param frame: 
        :param fltr:
        :return: 
        """

        # Get filter
        if fltr is None:
            if self.lazy: raise ValueError("This is a lazy frame list. The filter has to be passed explicitely")
            else: fltr = frame.filter

        # Check if not None
        if fltr is None: raise ValueError("Filter cannot be determined")

        # Parse filter
        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        # Check keys
        if fltr in self.frames: raise ValueError("Already a frame for the '" + str(fltr) + "' filter")

        # Call the function of the base class
        super(FrameList, self).append(fltr, frame)

    # -----------------------------------------------------------------

    def append_from_file(self, path, fltr=None):

        """
        This function ...
        :param path:
        :param fltr:
        :return:
        """

        if self.lazy: self.append(path, fltr=fltr)
        else: self.append(Frame.from_file(path), fltr=fltr)

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
        return self.bounding_box.center

    # -----------------------------------------------------------------

    @property
    def coordinate_systems(self):
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
        return self.max_ra.to("deg").value

    # -----------------------------------------------------------------

    @property
    def ra_center(self):
        return 0.5 * (self.min_ra + self.max_ra)

    # -----------------------------------------------------------------

    @property
    def ra_center_deg(self):
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
        dec_center = 0.5 * (self.min_dec + self.max_dec)
        return dec_center

    # -----------------------------------------------------------------

    @property
    def dec_center_deg(self):
        return self.dec_center.to("deg").value

    # -----------------------------------------------------------------

    @property
    def max_dec_deg(self):
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
        log.debug("Converting to unit " + tostr(unit, add_physical_type=True) + " ...")

        # Initialize list for converted frames
        new_frames = FrameList()

        # Convert all
        for fltr in self.frames:
            frame = self.frames[fltr]
            new_frames.append(frame.converted_to(unit, **kwargs))

        # Return the new set of frames
        return new_frames

    # -----------------------------------------------------------------

    def convolve_to_highest_fwhm(self, remote=None, ignore=None):

        """
        This function ...
        :param remote:
        :param ignore:
        :return:
        """

        new_frames = convolve_to_highest_fwhm(*self.values, names=self.filter_names, remote=remote, ignore=ignore)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def rebin_to_highest_pixelscale(self, remote=None, rebin_remote_threshold=None, in_place=False, ignore=None):

        """
        This function ...
        :param remote:
        :param rebin_remote_threshold:
        :param in_place:
        :param ignore:
        :return: 
        """

        # In place
        if in_place: rebin_to_highest_pixelscale(*self.values, names=self.filter_names, remote=remote, rebin_remote_threshold=rebin_remote_threshold, in_place=True, ignore=ignore)

        # Replace
        else:

            new_frames = rebin_to_highest_pixelscale(*self.values, names=self.filter_names, remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)
            self.remove_all()
            for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convolve_and_rebin(self, remote=None, rebin_remote_threshold=None):

        """
        This function ...
        :param remote:
        :param rebin_remote_threshold:
        :return:
        """

        new_frames = convolve_and_rebin(*self.values, names=self.filter_names, remote=remote, rebin_remote_threshold=rebin_remote_threshold)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convert_to_same_unit(self, unit=None, **kwargs):

        """
        This function ...
        :param unit:
        :param kwargs:
        :return: 
        """
        
        new_frames = convert_to_same_unit(*self.values, unit=unit, names=self.filter_names, **kwargs)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convolve_rebin_and_convert(self, unit=None, **kwargs):

        """
        This function ...
        :param unit: 
        :param kwargs: 
        :return: 
        """

        new_frames = convolve_rebin_and_convert(*self.values, unit=unit, names=self.filter_names, **kwargs)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def uniformize(self, unit=None, **kwargs):

        """
        This function is an alias for convolve_rebin_and_convert
        :param unit:
        :param kwargs:
        :return: 
        """

        return self.convolve_rebin_and_convert(unit=unit, **kwargs)

    # -----------------------------------------------------------------

    @property
    def uniform_properties(self):

        """
        This function ...
        :return:
        """

        unit, wcs, pixelscale, psf_filter, fwhm, distance = check_uniformity(*self.values)
        return unit, wcs, pixelscale, psf_filter, fwhm, distance

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return check_wcs(*self.values)

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].wcs = value

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        return check_pixelscale(*self.values)

    # -----------------------------------------------------------------

    @pixelscale.setter
    def pixelscale(self, value):

        """
        This function ...
        :return:
        """

        for key in self.keys: self[key].pixelscale = value

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return check_unit(*self.values)

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].unit = value

    # -----------------------------------------------------------------

    @property
    def psf_filter(self):

        """
        This function ...
        :return:
        """

        return check_psf_filter(*self.values)

    # -----------------------------------------------------------------

    @psf_filter.setter
    def psf_filter(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].psf_filter = value

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        return check_fwhm(*self.values)

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].fwhm = value

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return check_distance(*self.values)

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """
        #print(self.keys)
        for key in self.keys: self[key].distance = value

    # -----------------------------------------------------------------

    def normalize(self, to=1.0):

        """
        This function ...
        :return:
        """

        # Normalize each frame
        for key in self.keys: self[key].normalize(to=to)

# -----------------------------------------------------------------

class NamedFrameList(NamedList):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Lazy?
        lazy = kwargs.pop("lazy", False)
        if lazy: lazy_evaluator = Frame.from_file
        else: lazy_evaluator = None

        # Call the constructor of the base class
        super(NamedFrameList, self).__init__(lazy=lazy, lazy_evaluator=lazy_evaluator)

        # Add
        for name in kwargs: self.append(kwargs[name], name)

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
    def from_directory(cls, path, contains=None, plane=None):

        """
        This function ...
        :param path: 
        :param contains:
        :param plane:
        :return: 
        """

        new = cls()
        for path, name in fs.files_in_path(path, returns=["path", "name"], extension="fits", contains=contains):
            frame = Frame.from_file(path, plane=plane)
            new.append(frame, name=name)
        return new

    # -----------------------------------------------------------------

    def write_to_directory(self, path, update_path=False, replace=True):

        """
        This function ...
        :param path:
        :param update_path:
        :param replace:
        :return:
        """

        # Inform the user
        log.info("Writing all frames to the '" + path + "' directory ...")

        # Loop over the frames
        for name in self.names:

            # Get the frame
            frame = self[name]

            # Determine the path
            filepath = fs.join(path, name + ".fits")

            # Check
            if fs.is_file(filepath) and not replace: raise IOError("File '" + filepath + "' already exists")

            # Save the frame
            frame.saveto(filepath, update_path=update_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, **paths):

        """
        This function ...
        :param paths: 
        :return: 
        """

        # Initialize
        lazy = paths.pop("lazy", False)
        new = cls(lazy=paths.pop("lazy", False))

        # Add
        for name in paths: new.append_from_file(paths[name], name)

        # Return the frame list
        return new

    # -----------------------------------------------------------------

    @classmethod
    def from_dictionary(cls, dictionary):

        """
        This function ...
        :param dictionary:
        :param lazy:
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

        if name is None:
            if self.lazy: raise ValueError("This is a lazy frame list, the name should be specified")
            else: name = frame.name
        if name is None: raise ValueError("Frame does not have a name")

        # Call the append function of the base class
        super(NamedFrameList, self).append(name, frame)

    # -----------------------------------------------------------------

    def append_from_file(self, path, name=None):

        """
        This function ...
        :param path:
        :param name:
        :return:
        """

        if self.lazy: self.append(path, name=name)
        else: self.append(Frame.from_file(path), name=name)

    # -----------------------------------------------------------------

    def convolve_to_name(self, name, remote=None):

        """
        This function ...
        :param name:
        :param remote:
        :return:
        """

        # Get FWHM and PSF filter
        fwhm = self[name].fwhm
        psf_filter = self[name].psf_filter

        # Convolve and replace
        new_frames = convolve_to_fwhm(*self.values, names=self.names, fwhm=fwhm, filter=psf_filter, remote=remote)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convolve_to_filter(self, fltr, fwhm=None, remote=None):

        """
        This function ...
        :param fltr:
        :param fwhm:
        :param remote:
        :return:
        """

        # Convolve and replace
        new_frames = convolve_to_fwhm(*self.values, names=self.names, fwhm=fwhm, filter=fltr, remote=remote)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def rebin_to_name(self, name, remote=None, rebin_remote_threshold=None, in_place=False):

        """
        This function ...
        :param name:
        :param remote:
        :param rebin_remote_threshold:
        :param in_place:
        :return:
        """

        # Get pixelscale and wcs
        pixelscale = self[name].average_pixelscale
        wcs = self[name].wcs

        # Internal frames are not used elsewhere
        if in_place: rebin_to_pixelscale(*self.values, names=self.names, pixelscale=pixelscale, wcs=wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold, in_place=True)

        # Internal frames can be used elsewhere, and we don't want them rebinned
        else:

            # Rebin and replace
            new_frames = rebin_to_pixelscale(*self.values, names=self.names, pixelscale=pixelscale, wcs=wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold)
            self.remove_all()
            for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def rebin_to_wcs(self, wcs, remote=None, rebin_remote_threshold=None, in_place=False, ignore=None):

        """
        This function ...
        :param wcs:
        :param remote:
        :param rebin_remote_threshold:
        :param in_place:
        :param ignore:
        :return:
        """

        # Get pixelscale
        pixelscale = wcs.average_pixelscale

        # In place
        if in_place: rebin_to_pixelscale(*self.values, names=self.names, pixelscale=pixelscale, wcs=wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold, in_place=True, ignore=ignore)

        # Replace
        else:

            # Rebin and replace
            new_frames = rebin_to_pixelscale(*self.values, names=self.names, pixelscale=pixelscale, wcs=wcs, remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)
            self.remove_all()
            for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    @property
    def highest_fwhm_name(self):

        """
        This function ...
        :return:
        """

        return get_highest_fwhm_name(*self.values, names=self.names)

    # -----------------------------------------------------------------

    @property
    def highest_pixelscale_name(self):

        """
        This function ...
        :return:
        """

        return get_highest_pixelscale_name(*self.values, names=self.names)

    # -----------------------------------------------------------------

    def highest_fwhm_name_below(self, fwhm):

        """
        Thisf unction ...
        :param fwhm:
        :return:
        """

        return get_highest_fwhm_name(*self.values, names=self.names, below=fwhm)

    # -----------------------------------------------------------------

    def highest_pixelscale_name_below(self, pixelscale):

        """
        This function ...
        :param pixelscale:
        :return:
        """

        return get_highest_pixelscale_name(*self.values, names=self.names, below=pixelscale)

    # -----------------------------------------------------------------

    def highest_pixelscale_name_above_npixels(self, npixels):

        """
        This function ...
        :param npixels:
        :return:
        """

        return get_highest_pixelscale_name(*self.values, names=self.names, below_npixels=npixels)

    # -----------------------------------------------------------------

    def convolve_to_highest_fwhm(self, remote=None, ignore=None):

        """
        This function ...
        :param remote:
        :param ignore:
        :return:
        """

        new_frames = convolve_to_highest_fwhm(*self.values, names=self.names, remote=remote, ignore=ignore)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def rebin_to_highest_pixelscale(self, remote=None, rebin_remote_threshold=None, in_place=False, ignore=None):

        """
        This function ...
        :param remote:
        :param rebin_remote_threshold:
        :param in_place:
        :param ignore:
        :return:
        """

        # In place
        if in_place: rebin_to_highest_pixelscale(*self.values, names=self.names, remote=remote, rebin_remote_threshold=rebin_remote_threshold, in_place=True, ignore=ignore)

        # REplace
        else:

            new_frames = rebin_to_highest_pixelscale(*self.values, names=self.names, remote=remote, rebin_remote_threshold=rebin_remote_threshold, ignore=ignore)
            self.remove_all()
            for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convolve_and_rebin(self, remote=None, unitless=None):

        """
        This function ...
        :param remote:
        :param unitless:
        :return:
        """

        new_frames = convolve_and_rebin(*self.values, names=self.names, remote=remote, unitless=unitless)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convert_to_same_unit(self, unit=None, **kwargs):

        """
        This function ...
        :param unit:
        :param kwargs:
        :return:
        """

        new_frames = convert_to_same_unit(*self.values, unit=unit, names=self.names, **kwargs)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def convolve_rebin_and_convert(self, unit=None, **kwargs):

        """
        This function ...
        :param unit:
        :param kwargs:
        :return:
        """

        new_frames = convolve_rebin_and_convert(*self.values, unit=unit, names=self.names, **kwargs)
        self.remove_all()
        for frame in new_frames: self.append(frame)

    # -----------------------------------------------------------------

    def uniformize(self, unit=None, **kwargs):

        """
        This function is an alias for convolve_rebin_and_convert
        :param unit:
        :param kwargs:
        :return:
        """

        return self.convolve_rebin_and_convert(unit=unit, **kwargs)

    # -----------------------------------------------------------------

    @property
    def uniform_properties(self):

        """
        This function ...
        :return:
        """

        return self.get_uniform_properties()

    # -----------------------------------------------------------------

    def get_uniform_properties(self, strict=True):

        """
        This function ...
        :param strict:
        :return:
        """

        return check_uniformity(*self.values, strict=strict)

    # -----------------------------------------------------------------

    def set_uniform_properties(self):

        """
        This function ...
        :return:
        """

        # Get properties that are uniform (or only defined for one)
        unit, wcs, pixelscale, psf_filter, fwhm, distance = self.get_uniform_properties(strict=False)

        #print("unit", unit)
        #print("wcs", wcs)
        #print("pixelscale", pixelscale)
        #print("psf_filter", psf_filter)
        #print("fwhm", fwhm)
        #print("distance", distance)

        # Set
        self.unit = unit
        self.wcs = wcs
        self.pixelscale = pixelscale
        self.psf_filter = psf_filter
        self.fwhm = fwhm
        self.distance = distance

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return check_wcs(*self.values)

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].wcs = value

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        return check_pixelscale(*self.values)

    # -----------------------------------------------------------------

    @pixelscale.setter
    def pixelscale(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].pixelscale = value

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        return check_unit(*self.values)

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].unit = value

    # -----------------------------------------------------------------

    @property
    def psf_filter(self):

        """
        This function ...
        :return:
        """

        return check_psf_filter(*self.values)

    # -----------------------------------------------------------------

    @psf_filter.setter
    def psf_filter(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].psf_filter = value

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        return check_fwhm(*self.values)

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].fwhm = value

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return check_distance(*self.values)

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].distance = value

    # -----------------------------------------------------------------

    def normalize(self, to=1.0):

        """
        This function ...
        :return:
        """

        # Normalize each frame
        for key in self.keys: self[key].normalize(to=to)

    # -----------------------------------------------------------------

    def replace_nans(self, value):

        """
        THis function ...
        :param value:
        :return:
        """

        for key in self.keys: self[key].replace_nans(value)

    # -----------------------------------------------------------------

    def show_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        # Print wcss:
        print("")
        for name in self.names:
            print(name, self[name].wcs)
            print("")

# -----------------------------------------------------------------

class ImageList(FilterBasedList):
        
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
        super(ImageList, self).__init__()

        # Add frames
        for image in args: self.append(image)
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
        for path in fs.files_in_path(path): new.append(Image.from_file(path))
        return new

    # -----------------------------------------------------------------

    def append(self, image, fltr=None):

        """
        This function ...
        :param image:
        :param fltr:
        :return: 
        """

        fltr = fltr if fltr is not None else image.filter

        # Check keys
        if fltr in self.images: raise ValueError("Already an image for the '" + str(fltr) + "' filter")

        # Call the function of the base class
        super(ImageList, self).append(fltr, image)

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

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(NamedImageList, self).__init__()

        # Add
        for name in kwargs: self.append(kwargs[name], name)

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

def convolve_rebin_and_convert(*frames, **kwargs):

    """
    This function ...
    :param frames: 
    :param kwargs: 
    :return: 
    """

    # First rebin
    frames = rebin_to_highest_pixelscale(*frames, **kwargs)

    # Then convolve
    frames = convolve_to_highest_fwhm(*frames, **kwargs)

    # Then convert
    frames = convert_to_same_unit(*frames, **kwargs)

    # Return the frames
    return frames

# -----------------------------------------------------------------

def convolve_and_rebin(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # First rebin
    frames = rebin_to_highest_pixelscale(*frames, **kwargs)

    # Then convolve
    frames = convolve_to_highest_fwhm(*frames, **kwargs)

    # Return the frames
    return frames

# -----------------------------------------------------------------

def convert_to_same_unit(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return: 
    """

    from .datacube import DataCube

    # Get frame names
    names = kwargs.pop("names", None)

    # Get options for conversion
    distance = kwargs.pop("distance", None)
    density = kwargs.pop("density", False)
    brightness = kwargs.pop("brightness", False)
    density_strict = kwargs.pop("density_strict", False)
    brightness_strict = kwargs.pop("brightness_strict", False)
    wavelength = kwargs.pop("wavelength", None)

    # Inform the user
    log.info("Converting images to the same unit ...")

    # Check if the unit is defined
    if "unit" in kwargs and kwargs["unit"] is not None:
        unit = kwargs.pop("unit")
        #if types.is_string_type(unit): unit = u(unit, **kwargs) # not necessary: converted_to() of frame takes **kwargs
    else: unit = frames[0].unit
    if unit is None: raise ValueError("Target unit is None")

    # Remove unit from kwargs (if None stays in), otherwise frame.converted_to call will crash
    if "unit" in kwargs: kwargs.pop("unit")

    # Debugging
    log.debug("Converting images to unit " + tostr(unit, add_physical_type=True) + " ...")

    # Initialize list for converted frames
    new_frames = []

    # Convert all
    #index = 0
    for index, frame in enumerate(frames):

        # Get frame name
        #name = names[index] if names is not None else ""
        print_name = "'" + names[index] + "' " if names is not None else ""

        # Get type
        if isinstance(frame, Frame): image_type = "frame"
        elif isinstance(frame, DataCube): image_type = "datacube"
        elif isinstance(frame, Image): image_type = "image"
        else: raise ValueError("Invalid argument of type '" + str(type(frame)) + "'")

        # Check unit
        if frame.unit == unit:

            # Debugging
            log.debug("Frame " + print_name + "already has the target unit of '" + tostr(unit, add_physical_type=True) + "' and will not be converted")

            # Create copy
            converted = frame.copy()

        # Convert
        elif frame.unit is not None:

            # Debugging
            log.debug("Converting " + image_type + " " + print_name + "with unit " + tostr(frame.unit, add_physical_type=True) + " to " + tostr(unit, add_physical_type=True) + " ...")

            # Create converted version
            if image_type == "frame": converted = frame.converted_to(unit, distance=distance, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict, wavelength=wavelength)
            elif image_type == "datacube" or image_type == "image":
                if wavelength is not None: raise ValueError("Wavelength cannot be specified when datacubes/images are passed")
                converted = frame.converted_to(unit, distance=distance, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict, silent=True)
            else: raise ValueError("Invalid argument of type '" + str(type(frame)) + "'")

        # Unit is None
        else: raise ValueError("Unit of frame " + print_name + "is not defined")

        # Set name
        if names is not None: converted.name = names[index]

        # Add to the list of new frames
        new_frames.append(converted)

        # Increment index
        #index += 1

    # Return the new set of frames
    return new_frames

# -----------------------------------------------------------------

def get_highest_pixelscale_name(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get frame names
    names = kwargs.pop("names")

    below = kwargs.pop("below", None)
    below_npixels = kwargs.pop("below_npixels", None)

    highest_pixelscale = None
    highest_pixelscale_wcs = None
    highest_pixelscale_name = None

    # Loop over the frames
    for frame, name in zip(frames, names):

        wcs = frame.wcs

        if wcs is None: raise ValueError("Coordinate system of the " + name + " image is not defined")

        # SKIP?
        if below is not None and wcs.average_pixelscale > below: continue
        if below_npixels is not None and min(frame.xsize, frame.ysize) > below_npixels: continue

        if highest_pixelscale is None or wcs.average_pixelscale > highest_pixelscale:

            highest_pixelscale = wcs.average_pixelscale
            highest_pixelscale_wcs = wcs
            highest_pixelscale_name = name

    # Debugging
    log.debug("The frame with the highest FWHM is the '" + highest_pixelscale_name + "' frame ...")

    # Return the name
    return highest_pixelscale_name

# -----------------------------------------------------------------

def rebin_to_highest_pixelscale(*frames, **kwargs):

    """
    This function ...
    :param frames: 
    :param kwargs:
    :return: 
    """

    # Get frame names
    names = kwargs.pop("names", None)

    # Get the remote
    remote = kwargs.pop("remote", None)
    rebin_remote_threshold = kwargs.pop("rebin_remote_threshold", None)

    # In place?
    in_place = kwargs.pop("in_place", False)

    # Which are unitless?
    unitless = kwargs.pop("unitless", None)

    # Ignore?
    ignore = kwargs.pop("ignore", None)
    no_pixelscale = kwargs.pop("no_pixelscale", "error")

    # Get distance
    distance = kwargs.pop("distance", None)

    # Check
    if len(frames) == 1:

        # Success
        log.success("Only one frame: not rebinning")

        frame = frames[0]
        frame.name = names[0]
        return [frame]

    # Inform the user
    log.info("Rebinning frames to the coordinate system with the highest pixelscale ...")

    all_pixelscales = []
    highest_pixelscale = None
    highest_pixelscale_wcs = None
    highest_pixelscale_index = None
    xsizes = []
    ysizes = []

    _check_shapes = False

    # Loop over the frames
    for index, frame in enumerate(frames):

        # Ignore?
        if ignore is not None and names[index] in ignore: continue

        #wcs = frame.wcs
        #pixelscale = frame.average_pixelscale
        #if wcs is None:
        #if pixelscale is None:
        
        # Angular pixelscale
        #if frame.has_angular_pixelscale: pixelscale = frame.average_pixelscale

        # Physical pixelscale
        #elif frame.has_physical_pixelscale: pixelscale = frame.average_angular_pixelscale

        # Set frame distance if given
        if distance is not None: frame.distance = distance

        if frame.has_pixelscale: pixelscale = frame.average_angular_pixelscale

        # Pixelscale not defined
        else:

            if no_pixelscale == "error":
                if names is not None: raise ValueError("Pixelscale of the " + names[index] + " image is not defined")
                else: raise ValueError("Pixelscale of the image is not defined")
            elif no_pixelscale == "skip": continue
            elif no_pixelscale == "return": return frames
            elif no_pixelscale == "shape":
                _check_shapes = True
                xsizes.append(frame.xsize)
                ysizes.append(frame.ysize)
                continue
            else: raise ValueError("Invalid value for 'no_pixelscale'")

        #print(highest_pixelscale)
        #print(wcs.average_pixelscale)

        #all_pixelscales.append(wcs.average_pixelscale)
        all_pixelscales.append(pixelscale)
        xsizes.append(frame.xsize)
        ysizes.append(frame.ysize)

        #print(highest_pixelscale, pixelscale)
        #if highest_pixelscale is None or wcs.average_pixelscale > highest_pixelscale:
        if highest_pixelscale is None or pixelscale > highest_pixelscale:

            #highest_pixelscale = wcs.average_pixelscale
            highest_pixelscale = pixelscale
            highest_pixelscale_wcs = frame.wcs
            highest_pixelscale_index = index

    # Debugging
    log.debug("All pixelscales: " + tostr(all_pixelscales))

    # Check shapes?
    if _check_shapes:
        if sequences.all_equal(xsizes) and sequences.all_equal(ysizes): return frames
        else: raise ValueError("Some pixelscales are undefined, and not all frames have the same shape")

    # Debugging
    if names is not None: log.debug("The frame with the highest pixelscale is the '" + names[highest_pixelscale_index] + "' frame ...")
    log.debug("The highest pixelscale is " + tostr(highest_pixelscale))

    # Rebin
    return rebin_to_pixelscale(*frames, names=names, pixelscale=highest_pixelscale, wcs=highest_pixelscale_wcs,
                               remote=remote, rebin_remote_threshold=rebin_remote_threshold, in_place=in_place,
                               unitless=unitless, ignore=ignore)

# -----------------------------------------------------------------

def rebin_to_median_pixelscale(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get frame names
    names = kwargs.pop("names", None)

    # Get the remote
    remote = kwargs.pop("remote", None)

    # In-place
    in_place = kwargs.pop("in_place", False)

    # Check
    if len(frames) == 1:

        # Success
        log.success("Only one frame: not rebinning")

        frame = frames[0]
        frame.name = names[0]
        return [frame]

    # Inform the user
    log.info("Rebinning frames to the coordinate system with the median pixelscale ...")

    # Determine which frame contains the median pixelscale
    pixelscales = [frame.average_pixelscale.to("arcsec").value for frame in frames]
    median_index = np.argsort(pixelscales)[len(pixelscales) // 2]

    # Get pixelscale and WCS
    pixelscale = frames[median_index].pixelscale
    wcs = frames[median_index].wcs

    # Debugging
    if names is not None: log.debug("The frame with the median pixelscale is the '" + names[median_index] + "' frame ...")

    # Rebin
    return rebin_to_pixelscale(*frames, names=names, pixelscale=pixelscale, wcs=wcs, remote=remote, in_place=in_place)

# -----------------------------------------------------------------

def rebin_to_pixelscale(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # REMOTE?
    if "remote" in kwargs and kwargs["remote"] is not None:

        # Threshold is defined
        if "rebin_remote_threshold" in kwargs and kwargs["rebin_remote_threshold"] is not None:
            return rebin_to_pixelscale_local(*frames, **kwargs)

        # Threshold is not defined
        else:

            # Check that in_place is not enabled
            if kwargs.get("in_place", False): raise ValueError("Cannot enable 'in_place' for remote rebinning")

            # Rebin
            return rebin_to_pixelscale_remote(*frames, **kwargs) # all remote

    # LOCAL?
    else:

        # Threshold is defined: not possible because remote is not defined
        if "rebin_remote_threshold" in kwargs and kwargs["rebin_remote_threshold"] is not None:
            #raise ValueError("Cannot specify 'rebin_remote_threshold' if 'remote' is not defined")
            log.warning("'rebin_remote_threshold' is defined but 'remote' is not specified: rebinning locally ...")

        # Rebin the frames
        return rebin_to_pixelscale_local(*frames, **kwargs)

# -----------------------------------------------------------------

def rebin_to_pixelscale_remote(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get input
    names = kwargs.pop("names")
    #pixelscale = kwargs.pop("pixelscale")
    wcs = kwargs.pop("wcs")
    ignore = kwargs.pop("ignore", None)

    # Get remote
    remote = kwargs.pop("remote")
    remote = load_remote(remote)

    # Create temporary directory
    dirname = time.unique_name("rebin")
    temp_path = fs.create_directory_in(introspection.pts_temp_dir, dirname)

    # Save the frames
    for index, frame in enumerate(frames):
        name = names[index]
        if ignore is not None and name in ignore: continue # ignore?
        filepath = fs.join(temp_path, name + ".fits")
        frame.saveto(filepath)

    # Upload temporary directory to the remote
    remote_temp_path = remote.upload_directory_to(temp_path, remote.pts_temp_path, compress=True, show_output=log.is_debug)

    # Run PTS remotely
    output = execute_pts_remote(remote, "rebin", cwd=remote_temp_path, show_output=log.is_debug, wcs=wcs, backup=False)

    # Download the rebinned frames
    temp_path = remote.download_directory_to(remote_temp_path, introspection.pts_temp_dir)

    # Initialize list for rebinned frames
    new_frames = []

    # Load the frames
    for index, name in enumerate(names):

        # Ignored frame?
        if ignore is not None and name in ignore:

            new = frames[index].copy()
            new_frames.append(new)

        else:

            # Determine the filepath
            filepath = fs.join(temp_path, name + ".fits")

            # Load the frame
            frame = Frame.from_file(filepath)
            new_frames.append(frame)

    # Return the list of rebinned frames
    return new_frames

# -----------------------------------------------------------------

def any_frame_above_threshold(frames, threshold):

    """
    This function ...
    :param frames:
    :param threshold:
    :return:
    """

    for frame in frames:
        if frame.data_size > threshold: return True
    return False

# -----------------------------------------------------------------

def rebin_to_pixelscale_local(*frames, **kwargs):

    """
    THis function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get input
    names = kwargs.pop("names")
    highest_pixelscale = kwargs.pop("pixelscale")
    highest_pixelscale_wcs = kwargs.pop("wcs")

    # FOR SOME FRAMES TO BE REBINNED REMOTELY
    session = None
    remote = kwargs.pop("remote", None)
    rebin_remote_threshold = kwargs.pop("rebin_remote_threshold", None)

    # IN PLACE?
    in_place = kwargs.pop("in_place", False)

    # Which are unitless
    unitless = kwargs.pop("unitless", None)

    # Ignore?
    ignore = kwargs.pop("ignore", None)

    if rebin_remote_threshold is not None and remote is None:
        log.warning("'rebin_remote_threshold' is defined but 'remote' is not specified: rebinning locally ...")
        rebin_remote_threshold = None

    # Remote rebinning check
    if rebin_remote_threshold is not None:
        #if remote is None: raise ValueError("Cannot specify 'rebin_remote_threshold' when 'remote' is not specified")
        # Make remote session, ONLY IF IT WILL BE NECESSARY
        if any_frame_above_threshold(frames, rebin_remote_threshold):

            # Make remote
            remote = load_remote(remote)

            # START SESSION
            new_connection = False
            session = remote.start_python_session(attached=True, new_connection_for_attached=new_connection)

    # If session is not None, there are frames for which remote rebinning will be used
    if session is not None and in_place:
        raise ValueError("in_place cannot be enabled when there are frames to be rebinned remotely")

    # Initialize list for rebinned frames
    if in_place: new_frames = None
    else: new_frames = []

    # Rebin
    index = 0
    for frame in frames:

        # Determine frame name
        name = names[index] if names is not None else ""
        print_name = "'" + names[index] + "' " if names is not None else ""

        # Ignore?
        if ignore is not None and name in ignore:

            # Debugging
            if names is not None: log.debug("Frame " + print_name + "will be ignored for rebinning")

            # Not in place: create copy
            if not in_place:

                # Create new and set the name
                new = frame.copy()
                if names is not None: new.name = names[index]

                # Add
                new_frames.append(new)

        # If the current frame is the frame with the highest pixelscale
        elif frame.wcs == highest_pixelscale_wcs:

            # Debugging
            log.debug("Frame " + print_name + "has highest pixelscale of '" + tostr(highest_pixelscale) + "' and is not rebinned")

            # Not in place, create copy
            if not in_place:

                # Create new and set the name
                new = frame.copy()
                if names is not None: new.name = names[index]

                # Add
                new_frames.append(new)

        # The frame has a lower pixelscale, has to be rebinned
        else:

            # Debugging
            log.debug("Frame " + print_name + "will be rebinned ...")

            # Unitless frame
            if unitless is not None: unitless_frame = name in unitless
            else: unitless_frame = False

            # In place?
            if in_place: rebin_frame(name, frame, highest_pixelscale_wcs, rebin_remote_threshold=rebin_remote_threshold,
                                     session=session, in_place=True, unitless=unitless_frame)

            # New frames
            else:

                # Create rebinned frame
                rebinned = rebin_frame(name, frame, highest_pixelscale_wcs, rebin_remote_threshold=rebin_remote_threshold,
                                       session=session, unitless=unitless_frame)

                # Set the name
                if names is not None: rebinned.name = names[index]

                # Add the rebinned frame
                new_frames.append(rebinned)

        # Increment the index for the frames
        index += 1

    # Return the rebinned frames
    if not in_place: return new_frames

# -----------------------------------------------------------------

def rebin_frame(name, frame, wcs, rebin_remote_threshold=None, session=None, in_place=False, unitless=False):

    """
    This function ...
    :param name:
    :param frame:
    :param wcs:
    :param rebin_remote_threshold:
    :param session:
    :param in_place:
    :param unitless:
    :return:
    """

    # Debugging
    if frame.unit is not None: log.debug("Rebinning frame '" + name + "' with unit " + tostr(frame.unit, add_physical_type=True) + " ...")
    else: log.debug("Rebinning frame '" + name + "' without unit ...")

    # CONVERT TO PER ANGULAR OR INTRINSIC AREA, IF UNIT IS DEFINED
    if frame.unit is not None and frame.is_photometric:

        # Convert to the corresponding brightness unit
        #original_unit = frame.convert_to_corresponding_brightness_unit()
        original_unit = frame.unit
        conversion_factor = frame.convert_to_corresponding_angular_or_intrinsic_area_unit()
        correction_factor = None

        # Converted?
        if original_unit != frame.unit: log.debug("Unit has been converted to " + tostr(frame.unit, add_physical_type=True) + " prior to rebinning (conversion factor of " + str(conversion_factor) + ")")

    # UNIT IS DEFINED, NOT PHOTOMETRIC
    elif frame.unit is not None:

        conversion_factor = None
        original_unit = None

        #print(wcs.pixelarea, frame.pixel_solid_angle, frame.pixel_area)
        if wcs.is_celestial: correction_factor = (wcs.pixelarea / frame.pixel_solid_angle).to("").value
        elif wcs.is_physical: correction_factor = (wcs.pixelarea / frame.pixel_area).to("").value
        else: raise ValueError("Uknown coordinate system type")

        #print(type(correction_factor))
        #exit()

        # Debugging
        log.debug("Frame '" + name + "' will be corrected after rebinning with a factor of " + str(correction_factor) + " to account for the changing pixelscale ...")

    # UNIT IS NOT DEFINED
    else:
        if not unitless: log.warning("The unit of the '" + name + "' frame is not defined: make sure it is per unit of angular or intrinsic area (or only the relative variation is important)")
        original_unit = None
        conversion_factor = None
        correction_factor = None

    # REBIN remotely
    if rebin_remote_threshold is not None and frame.data_size > rebin_remote_threshold:

        # In place is not possible
        if in_place: raise ValueError("Cannot enable 'in_place' when frames are remotely rebinned")

        # Debugging
        log.debug("Rebinning frame remotely ...")

        # Remote rebinning
        from .remote import RemoteFrame
        remoteframe = RemoteFrame.from_local(frame, session)
        remoteframe.rebin(wcs, convert=False)  # should already be converted, or correction factor to be applied later
        rebinned = remoteframe.to_local()

    # REBIN locally
    else:

        # Debugging
        log.debug("Rebinning frame locally ...")

        #print(frame.name, frame.wcs)
        #print(wcs)

        if in_place:
            frame.rebin(wcs, convert=False) # should already be converted, or correction factor to be applied later
            rebinned = None
        else: rebinned = frame.rebinned(wcs, convert=False) # should already be converted, or correction factor to be applied later

    # IF THERE WAS AN ORIGINAL UNIT
    if original_unit is not None:

        # Debugging
        log.debug("Converting the unit back to the original unit " + tostr(original_unit, add_physical_type=True))

        # Convert the original back to the original unit
        frame.convert_to(original_unit)

        # CONVERT THE RESULT AS WELL IF NOT 'IN PLACE'
        if not in_place: rebinned.convert_to(original_unit)

    # CORRECT
    elif correction_factor is not None: rebinned *= correction_factor

    # Return rebinned frame
    if not in_place: return rebinned

# -----------------------------------------------------------------

# def rebin_frame_old(name, frame, wcs, rebin_remote_threshold=None, session=None, in_place=False):
#
#     """
#     THis function ...
#     :param name:
#     :param frame:
#     :param wcs:
#     :param rebin_remote_threshold:
#     :param session:
#     :param in_place:
#     :return:
#     """
#
#     # Is per pixelsize
#     if frame.unit is not None and frame.unit.is_per_pixelsize:
#
#         # Debugging
#         log.debug("Frame " + name + "is expressed in units per angular or intrinsic area (pixelsize squared)")
#
#         # Debugging
#         log.debug("Rebinning frame " + name + "with unit " + tostr(frame.unit, add_physical_type=True) + " ...")
#
#         # REBIN remotely
#         if rebin_remote_threshold is not None and frame.data_size > rebin_remote_threshold:
#
#             # In place is not possible
#             if in_place: raise ValueError("Cannot enable 'in_place' when frames are remotely rebinned")
#
#             from .remote import RemoteFrame
#             remoteframe = RemoteFrame.from_local(frame, session)
#             remoteframe.rebin(wcs)
#             rebinned = remoteframe.to_local()
#
#         # REBIN locally
#         else:
#
#             # IN PLACE
#             if in_place:
#
#                 frame.rebin(wcs)
#                 rebinned = None
#
#             # NOT IN PLACE
#             else: rebinned = frame.rebinned(wcs)
#
#     # Not per pixelsize
#     else:
#
#         # Debugging
#         log.debug("Frame " + name + "is not expressed in units per angular or intrinsic area (pixelsize squared)")
#
#         # Debugging
#         # log.debug("Converting frame " + name + "with unit " + str(frame.unit) + " to " + str(frame.corresponding_angular_area_unit) + " prior to rebinning ...")
#         # old_unit = frame.unit
#         # rebinned = frame.converted_to_corresponding_angular_area_unit(**kwargs)
#         # rebinned.rebin(highest_pixelscale_wcs)
#         # Convert back to old unit
#         # rebinned.convert_to(old_unit)
#
#         # print(rebinned)
#
#         # NEW WAY:
#
#         # Converting unit is not necessary if we calculate the ratio of both pixel areas
#         ratio = wcs.pixelarea / frame.wcs.pixelarea
#
#         # Debugging
#         log.debug("Rebinning frame " + name + "and multiplying with a factor of " + str(ratio) + " to correct for the changing pixelscale ...")
#
#         # REBIN
#         # print("threshold", rebin_remote_threshold)
#         # print("FILESIZE", frame.file_size)
#         # print("DATASIZE", frame.data_size)
#         # if rebin_remote_threshold is not None and frame.file_size > rebin_remote_threshold:
#         if rebin_remote_threshold is not None and frame.data_size > rebin_remote_threshold:
#
#             # In place is not possible
#             if in_place: raise ValueError("Cannot enable 'in_place' when frames are remotely rebinned")
#
#             from .remote import RemoteFrame
#             remoteframe = RemoteFrame.from_local(frame, session)
#             remoteframe.rebin(wcs)
#             rebinned = remoteframe.to_local()
#
#             # Multiply with ratio
#             if ratio != 1.0: rebinned *= float(ratio)
#
#         # LOCAL
#         else:
#
#             # Rebin and multiply
#             try:
#
#                 # IN PLACE
#                 if in_place:
#
#                     # Rebin in place
#                     frame.rebin(wcs)
#                     rebinned = None
#
#                     # Multiply with ratio
#                     if ratio != 1.0: frame *= float(ratio)
#
#                 # NOT IN PLACE
#                 else:
#
#                     # Create rebinned frame
#                     rebinned = frame.rebinned(wcs)
#
#                     # Multiply with ratio
#                     if ratio != 1.0: rebinned *= float(ratio)
#
#             # Something went wrong
#             except ValueError as e:
#
#                 print("")
#                 print("INPUT WCS:", frame.wcs)
#                 print("")
#                 print("OUTPUT WCS:", wcs)
#                 print("")
#                 raise RuntimeError("Rebinning the " + name + " image failed: " + str(e))
#
#         # Multiply with ratio
#         #if ratio != 1.0: rebinned *= float(ratio)
#
#     # Return the rebinned frame
#     if not in_place: return rebinned

# -----------------------------------------------------------------

def get_highest_fwhm_name(*frames, **kwargs):

    """
    THis function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get frame names
    names = kwargs.pop("names")

    below = kwargs.pop("below", None)

    highest_fwhm = None
    highest_fwhm_filter = None
    highest_fwhm_name = None

    # Loop over the frames
    for frame, name in zip(frames, names):

        # Search and set frame FWHM
        frame_fwhm = frame.fwhm
        if frame_fwhm is None: frame_fwhm = get_fwhm(frame.filter)
        frame.fwhm = frame_fwhm

        # Check again
        if frame.fwhm is None: raise ValueError("FWHM of the " + name + " image cannot be determined")

        # SKIP
        if below is not None and frame.fwhm > below: continue

        if highest_fwhm is None or frame.fwhm > highest_fwhm:

            highest_fwhm = frame.fwhm
            highest_fwhm_filter = frame.psf_filter
            highest_fwhm_name = name

    # Debugging
    log.debug("The frame with the highest FWHM is the '" + highest_fwhm_name + "' frame ...")

    # Return the name
    return highest_fwhm_name

# -----------------------------------------------------------------

def convolve_to_highest_fwhm(*frames, **kwargs):

    """
    This function ...
    :param frames: 
    :return: 
    """

    # Get frame names
    names = kwargs.pop("names", None)

    # Get remote
    remote = kwargs.pop("remote", None)

    # Get ignore
    ignore = kwargs.pop("ignore", None)
    no_fwhm = kwargs.pop("no_fwhm", "error")

    # Check
    if len(frames) == 1:

        # Success
        log.success("Only one frame: not convolving")

        frame = frames[0]
        frame.name = names[0]
        return [frame]

    # Inform the user
    log.info("Convolving frames to the resolution of the frame with the highest FWHM ...")

    highest_fwhm = None
    highest_fwhm_filter = None
    highest_fwhm_index = None

    all_fwhms = []

    # Loop over the frames
    for index, frame in enumerate(frames):

        # Ignore?
        if ignore and names[index] in ignore: continue

        # Search and set frame FWHM
        frame_fwhm = frame.fwhm
        if frame_fwhm is None:

            if frame.psf_filter is None:

                if no_fwhm == "error":
                    if names is not None: raise ValueError("Neither FWHM nor filter of the " + names[index] + " image is defined")
                    else: raise ValueError("Neither FWHM nor filter of the frame is defined")
                elif no_fwhm == "skip": continue
                elif no_fwhm == "return": return frames
                else: raise ValueError("Invalid value for 'no_fwhm'")

            #frame_fwhm = get_fwhm(frame.psf_filter)
            if has_variable_fwhm(frame.psf_filter):
                log.warning("Using average value for the FWHM for " + str(frame.psf_filter) + " images from Clark et al. (2017) to proceed ...")
                frame_fwhm = get_average_variable_fwhm(frame.psf_filter)
            else: frame_fwhm = get_fwhm(frame.psf_filter)

        frame.fwhm = frame_fwhm

        # Check again
        if frame.fwhm is None:

            if no_fwhm == "error":
                if names is not None: raise ValueError("FWHM of the " + names[index] + " image cannot be determined")
                else: raise ValueError("FWHM of the image cannot be determined")
            elif no_fwhm == "skip": continue
            elif no_fwhm == "return": return frames
            else: raise ValueError("Invalid value for 'no_fwhm'")

        # Add to list
        else: all_fwhms.append(frame_fwhm)

        if highest_fwhm is None or frame.fwhm > highest_fwhm:

            highest_fwhm = frame.fwhm
            highest_fwhm_filter = frame.psf_filter
            highest_fwhm_index = index

    # Show all FWHMs
    log.debug("All FWHMs: " + tostr(all_fwhms))

    # Debugging
    if names is not None: log.debug("The frame with the highest FWHM is the '" + names[highest_fwhm_index] + "' frame ...")

    # Convolve
    return convolve_to_fwhm(*frames, names=names, fwhm=highest_fwhm, filter=highest_fwhm_filter, remote=remote, ignore=ignore)

# -----------------------------------------------------------------

def convolve_to_fwhm(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    if "remote" in kwargs and kwargs["remote"] is not None: return convolve_to_fwhm_remote(*frames, **kwargs)
    else: return convolve_to_fwhm_local(*frames, **kwargs)

# -----------------------------------------------------------------

def convolve_to_fwhm_remote(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get input
    names = kwargs.pop("names")
    fwhm = kwargs.pop("fwhm")
    fltr = kwargs.pop("filter")
    ignore = kwargs.pop("ignore", None)

    # Get remote
    remote = kwargs.pop("remote")
    remote = load_remote(remote)

    # Create temporary directory
    dirname = time.unique_name("convolve")
    temp_path = fs.create_directory_in(introspection.pts_temp_dir, dirname)

    # Save the frames
    for index, frame in enumerate(frames):
        name = names[index]
        if ignore is not None and name in ignore: continue
        filepath = fs.join(temp_path, name + ".fits")
        frame.saveto(filepath)

    # Upload temporary directory to the remote
    remote_temp_path = remote.upload_directory_to(temp_path, remote.pts_temp_path, compress=True, show_output=log.is_debug)

    # Run PTS remotely
    output = execute_pts_remote(remote, "convolve", cwd=remote_temp_path, show_output=log.is_debug, filter=fltr, fwhm=fwhm, backup=False)

    # Download the rebinned frames
    temp_path = remote.download_directory_to(remote_temp_path, introspection.pts_temp_dir)

    # Initialize list for rebinned frames
    new_frames = []

    # Load the frames
    for index, name in enumerate(names):

        # Ignored frame?
        if ignore is not None and name in ignore:

            # Add copy
            new = frames[index].copy()
            new_frames.append(new)

        # Not ignored
        else:

            # Determine the filepath
            filepath = fs.join(temp_path, name + ".fits")

            # Load the frame
            frame = Frame.from_file(filepath)
            new_frames.append(frame)

    # Return the list of rebinned frames
    return new_frames

# -----------------------------------------------------------------

def convolve_to_fwhm_local(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get input
    names = kwargs.pop("names")
    highest_fwhm = kwargs.pop("fwhm")
    highest_fwhm_filter = kwargs.pop("filter")
    ignore = kwargs.pop("ignore", None)

    # Get kernel services
    aniano = AnianoKernels()
    matching = MatchingKernels()

    # Initialize list for convolved frames
    new_frames = []

    # Convolve
    index = 0
    for frame in frames:

        # Get frame name
        name = names[index] if names is not None else ""
        print_name = "'" + names[index] + "' " if names is not None else ""

        #print(print_name, frame.psf_filter, frame.fwhm, highest_fwhm_filter, highest_fwhm)
        frame_psf_filter_defined = frame.psf_filter is not None
        highest_fwhm_filter_defined = highest_fwhm_filter is not None
        same_psf_filter = frame_psf_filter_defined and highest_fwhm_filter_defined and frame.psf_filter == highest_fwhm_filter
        frame_fwhm_defined = frame.fwhm is not None
        highest_fwhm_defined = highest_fwhm is not None
        same_fwhm = frame_fwhm_defined and highest_fwhm_defined and frame.fwhm == highest_fwhm
        both_filters_defined = frame_psf_filter_defined and highest_fwhm_filter_defined
        some_filters_undefined = not both_filters_defined

        # Ignore?
        if ignore is not None and name in ignore:

            # Debugging
            log.debug("Frame " + print_name + "will be ignored for convolution")

            # Make new frame and set the name
            new = frame.copy()
            if names is not None: new.name = names[index]

            # Add a copy of the frame
            new_frames.append(new)

        # Doesn't need convolution
        #elif frame.psf_filter == highest_fwhm_filter:
        # elif same_psf_filter:
        #
        #     # Debugging
        #     log.debug("Frame " + print_name + "has highest FWHM of " + str(highest_fwhm) + " and will not be convolved")
        #
        #     # Make new frame and set the name
        #     new = frame.copy()
        #     if names is not None: new.name = names[index]
        #
        #     # Add a copy of the frame
        #     new_frames.append(new)

        # Same FWHM, and potentially the same filter?
        elif same_fwhm and (some_filters_undefined or same_psf_filter):

            # Debugging
            log.debug("Frame " + print_name + "has highest FWHM of " + str(highest_fwhm) + " and will not be convolved")

            # Make new frame and set the name
            new = frame.copy()
            if names is not None: new.name = names[index]

            # Add a copy of the frame
            new_frames.append(new)

        # Convolve
        else:

            # Debugging
            log.debug("Frame " + print_name + "will be convolved to a PSF with FWHM = " + str(highest_fwhm) + " ...")

            # Get the kernel, either from aniano or from matching kernels
            #print(frame.psf_filter, highest_fwhm_filter)
            #print(frame.filter, frame.psf_filter)
            #fwhm1 = get_fwhm(frame.psf_filter)
            #fwhm2 = get_fwhm(highest_fwhm_filter)
            #print(fwhm1, fwhm2, frame.fwhm, highest_fwhm)

            if has_variable_fwhm(frame.psf_filter): from_standard_resolution = True # assume OK
            else: from_standard_resolution = frame.fwhm == get_fwhm(frame.psf_filter)

            if has_variable_fwhm(highest_fwhm_filter): to_standard_resolution = True # assume OK
            else: to_standard_resolution = highest_fwhm == get_fwhm(highest_fwhm_filter)

            #print(from_standard_resolution, to_standard_resolution)
            #print(highest_fwhm, highest_fwhm_filter, get_fwhm(highest_fwhm_filter))
            if not from_standard_resolution and frame.has_psf_filter: log.warning("PSF filter is defined in frame (" + tostr(frame.psf_filter) + ") but it does not seem to be correct comparing to the frame's FWHM (" + tostr(frame.fwhm) + ")")
            all_standard = from_standard_resolution and to_standard_resolution

            if all_standard and aniano.has_kernel_for_filters(frame.psf_filter, highest_fwhm_filter): kernel = aniano.get_kernel(frame.psf_filter, highest_fwhm_filter, from_fwhm=frame.fwhm, to_fwhm=highest_fwhm)
            else:

                # Get from and to filter
                from_filter = frame.psf_filter
                to_filter = highest_fwhm_filter

                # Get from and to FWHM
                if frame.fwhm is not None: from_fwhm = frame.fwhm
                elif has_variable_fwhm(from_filter):
                    if has_average_variable_fwhm(from_filter): from_fwhm = get_average_variable_fwhm(from_filter)
                    else: raise ValueError("FWHM for the " + tostr(from_filter) + " frame is undefined")
                else: from_fwhm = get_fwhm(from_filter)
                to_fwhm = highest_fwhm

                # Generate the kernel
                kernel = matching.get_kernel(from_filter, to_filter, frame.angular_pixelscale, from_fwhm=from_fwhm, to_fwhm=to_fwhm)
                #kernel.saveto("kernel.fits")

            # Convolve with the kernel
            convolved = frame.convolved(kernel)

            # Set the PSF filter to be sure
            convolved.psf_filter = highest_fwhm_filter

            # Set the name
            if names is not None: convolved.name = names[index]

            # Add to the list
            new_frames.append(convolved)

        # Increment the index
        index += 1

    # Return the convolved frames
    return new_frames

# -----------------------------------------------------------------

def check_uniformity(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get unit
    unit = check_unit(*frames, **kwargs)

    # Get the wcs
    wcs = check_wcs(*frames, **kwargs)

    # Get the pixelscale
    pixelscale = check_pixelscale(*frames, **kwargs)

    # Get the PSF filter
    psf_filter = check_psf_filter(*frames, **kwargs)

    # Get the FWHM
    fwhm = check_fwhm(*frames, **kwargs)

    # Get the distance
    distance = check_distance(*frames, **kwargs)

    # Return the common properties
    return unit, wcs, pixelscale, psf_filter, fwhm, distance

# -----------------------------------------------------------------

def check_unit(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param strict:
    :return:
    """

    # Check units
    units = [frame.unit for frame in frames]
    if kwargs.pop("strict", True) and not sequences.all_equal(units, ignore_none=True):
        # raise ValueError("Frames have to be in the same unit")
        log.error("Frames have to be in the same unit")
        log.error("Units:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.unit))
        print("")
        exit()
    unit = sequences.find_first_not_none(units)
    return unit

# -----------------------------------------------------------------

def check_wcs(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param strict:
    :return:
    """

    # Get WCS
    wcss = [frame.wcs for frame in frames]
    # print(wcss)
    if kwargs.pop("strict", True) and not sequences.all_equal(wcss, ignore_none=True):
        # raise ValueError("Frames have to be transformed to same pixel grid")
        log.error("Frames have to be transformed to the same pixel grid")
        log.error("Coordinate systems:")
        print("")
        for frame in frames:
            print(" - " + frame.name + ": " + str(frame.wcs))
            print("")
        #print("")
        exit()
    wcs = sequences.find_first_not_none(wcss)
    return wcs

# -----------------------------------------------------------------

def check_pixelscale(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param strict:
    :return:
    """

    # Get pixelscale
    pixelscales = [frame.average_pixelscale for frame in frames]
    if kwargs.pop("strict", True) and not sequences.all_close(pixelscales, ignore_none=True):
        # raise ValueError("Frames must have the same pixelscale")
        log.error("Frames must have the same pixelscale")
        log.error("Pixelscales:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.average_pixelscale))
        print("")
        exit()
    pixelscale = sequences.find_first_not_none(pixelscales)
    return pixelscale

# -----------------------------------------------------------------

def check_psf_filter(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return:
    """

    # Get PSF filter
    psf_filters = [frame.psf_filter for frame in frames]
    if kwargs.pop("strict", True) and not sequences.all_equal(psf_filters):
        # raise ValueError("Frames have to be convolved to the same resolution")
        log.error("Frames have to be convolved to the same resolution")
        log.error("PSF filters:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.psf_filter))
        print("")
        exit()
    psf_filter = sequences.find_first_not_none(psf_filters, return_none=True)
    return psf_filter

# -----------------------------------------------------------------

def check_fwhm(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param strict:
    :return:
    """

    # Get FWHM
    fwhms = [frame.fwhm for frame in frames]
    if sequences.all_none(fwhms): return None

    # Check FWHM
    difference = abs(fwhms[0] - fwhms[1])
    rel_difference = difference / fwhms[0]
    if rel_difference < 0.1 and kwargs.pop("strict", True) and not sequences.all_close(fwhms, ignore_none=True):
        log.warning("The FWHM difference is" + str(rel_difference * 100) + "% from the highest FWHM frame")
        log.warning("FWHMs:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.fwhm))
        print("")
    elif rel_difference > 0.1 and kwargs.pop("strict", True) and not sequences.all_close(fwhms, ignore_none=True):
        log.error("Frames have to have the same FWHM")
        log.error("FWHMs:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.fwhm))
        print("")
        exit()
    fwhm = sequences.find_first_not_none(fwhms)
    return fwhm

# -----------------------------------------------------------------

def check_distance(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param strict:
    :return:
    """

    # Get distance
    distances = [frame.distance for frame in frames]
    # print(distances)
    if sequences.all_none(distances):
        log.warning("Distances of the frames are undefined")
        return None
    elif kwargs.pop("strict", True) and not sequences.all_close(distances, ignore_none=True):
        # raise ValueError("Frames have to have the same distance to the object")
        log.error("Frames have to have the same distance to the object")
        log.error("Distances:")
        print("")
        for frame in frames: print(" - " + frame.name + ": " + str(frame.distance))
        print("")
        exit()
    distance = sequences.find_first_not_none(distances)
    return distance

# -----------------------------------------------------------------

def uniformize(*frames, **kwargs):

    """
    This function ...
    :param frames:
    :param kwargs:
    :return: 
    """

    # Add flags
    if "no_fwhm" not in kwargs: kwargs["no_fwhm"] = "return"
    if "no_pixelscale" not in kwargs: kwargs["no_pixelscale"] = "shape"
    #kwargs["distance"] = self.galaxy_distance

    # Get rebin for which frames
    do_rebin = True
    rebin_indices = None
    if "rebin" in kwargs:
        rebin = kwargs.pop("rebin")
        if types.is_integer_sequence_or_tuple(rebin): rebin_indices = rebin
        elif types.is_boolean_type(rebin): do_rebin = rebin
        else: raise ValueError("Invalid argument for 'rebin'")

    # Get convolve for which frames
    do_convolve = True
    convolve_indices = None
    if "convolve" in kwargs:
        convolve = kwargs.pop("convolve")
        if types.is_integer_sequence_or_tuple(convolve): convolve_indices = convolve
        elif types.is_boolean_type(convolve): do_convolve = convolve
        else: raise ValueError("Invalid argument for 'convolve'")

    # Get convert for which frames
    do_convert = True
    convert_indices = None
    if "convert" in kwargs:
        convert = kwargs.pop("convert")
        if types.is_integer_sequence_or_tuple(convert): convert_indices = convert
        elif types.is_boolean_type(convert): do_convert = convert
        else: raise ValueError("Invalid argument for 'convert'")

    # Rebin?
    if do_rebin:

        log.debug("Images will be rebinned to the highest pixelscale ...")

        if rebin_indices is not None:

            to_rebin_frames = sequences.subset(frames, rebin_indices)
            rebinned_frames = rebin_to_highest_pixelscale(*to_rebin_frames, **kwargs)
            sequences.put(rebinned_frames, frames, rebin_indices)

        else: frames = rebin_to_highest_pixelscale(*frames, **kwargs)

    # Convolve?
    if do_convolve:

        log.debug("Images will be convolved to the highest FWHM ...")

        if convolve_indices is not None:

            to_convolve_frames = sequences.subset(frames, convolve_indices)
            convolved_frames = convolve_to_highest_fwhm(*to_convolve_frames, **kwargs)
            sequences.put(convolved_frames, frames, convolve_indices)

        else: frames = convolve_to_highest_fwhm(*frames, **kwargs)

    # Convert?
    if do_convert:

        log.debug("Images will be converted to the same unit ...")

        if convert_indices is not None:

            to_convert_frames = sequences.subset(frames, convert_indices)
            converted_frames = convert_to_same_unit(*to_convert_frames, **kwargs)
            sequences.put(converted_frames, frames, convert_indices)

        else: frames = convert_to_same_unit(*frames, **kwargs)

    # Return the frames
    return frames

# -----------------------------------------------------------------
