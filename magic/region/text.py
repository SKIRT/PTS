#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.text Contains the TextRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion
from ..basics.mask import Mask
from .region import add_info, make_text_template, coordsys_name_mapping

# -----------------------------------------------------------------

class TextRegion(Region):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set the attributes
        self.center = center
        self.text = text

        # Call the constructor of the base class
        super(TextRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelTextRegion(TextRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param x:
        :param y:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, PixelCoordinate): raise ValueError("Center must be pixel coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        x = self.center.x
        y = self.center.y
        text = self.text

        string = prefix + make_text_template(fmt).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class SkyTextRegion(TextRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param ra:
        :param dec:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, SkyCoordinate): raise ValueError("Center must be sky coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        # convert coordsys string to coordsys object
        if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
        else: frame = None  # for pixel/image/physical frames

        x = float(self.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(self.center.transform_to(frame).spherical.lat.to('deg').value)
        text = self.text

        string = prefix + make_text_template(fmt).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class PhysicalTextRegion(TextRegion, PhysicalCoordinate):

    """
    This class ...
    """

    def __init__(self, center, text, **kwargs):

        """
        This function ...
        :param center:
        :param text:
        :param kwargs:
        """

        # Check the center coordinate
        if not isinstance(center, PhysicalCoordinate): raise ValueError("Center must be physical coordinate")

        # Call the constructor of TextRegion class
        TextRegion.__init__(self, center, text, **kwargs)

# -----------------------------------------------------------------
