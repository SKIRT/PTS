#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.region Contains the (abstract) Region class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import string
from abc import ABCMeta

# Import astronomical modules
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from ...core.tools import types

# -----------------------------------------------------------------

def make_point_template(fmt):

    """
    This function ...
    :param fmt:
    :return:
    """

    return 'point({x:' + fmt + '},{y:' + fmt + '})'

# -----------------------------------------------------------------

def make_line_template(fmt):

    """
    Thisf ucntion ...
    :return:
    """

    return 'line({x1:' + fmt + '},{y1:' + fmt + '},{x2:' + fmt + '},{y2:' + fmt + '}) # line=0 0'

# -----------------------------------------------------------------

def make_vector_template(fmt, radunitstr):

    """
    This function ...
    :param fmtr:
    :param radunitstr:
    :return:
    """

    return '# vector({x:' + fmt + '},{y:' + fmt + '},{l:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '}) vector=1'

# -----------------------------------------------------------------

def make_circle_template(fmt, radunitstr):

    """
    This function ...
    :param fmt:
    :param radunitstr:
    :return:
    """

    return 'circle({x:' + fmt + '},{y:' + fmt + '},{r:' + fmt + '}' + radunitstr + ')'

# -----------------------------------------------------------------

def make_ellipse_template(fmt, radunitstr, hmsdms=False):

    """
    This functio n...
    :param fmtr:
    :param radunitstr:
    :param hmsdms:
    :return:
    """

    if hmsdms: return 'ellipse({x},{y},{r1:' + fmt + '}' + radunitstr + ',{r2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})'
    else: return 'ellipse({x:' + fmt + '},{y:' + fmt + '},{r1:' + fmt + '}' + radunitstr + ',{r2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})'

# -----------------------------------------------------------------

def make_rectangle_template(fmt, radunitstr):

    """
    This function ...
    :param fmt:
    :param radunitstr:
    :return:
    """

    return 'box({x:' + fmt + '},{y:' + fmt + '},{d1:' + fmt + '}' + radunitstr + ',{d2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})'

# -----------------------------------------------------------------

def make_polygon_template():

    """
    Thisn function ...
    :return:
    """

    return 'polygon({c})'

# -----------------------------------------------------------------

def make_text_template(fmt):

    """
    This function ...
    :param fmt:
    :return:
    """

    return '# text({x:' + fmt + '},{y:' + fmt + '}) text="{text:}"'

# -----------------------------------------------------------------

def make_composite_template(fmt):

    """
    This function ...
    :param fmt:
    :return:
    """

    return '# composite({x:' + fmt + '},{y:' + fmt + '},{ang:' + fmt + '}) || composite=1'

# -----------------------------------------------------------------

def add_info(string, reg):

    """
    This function ...
    :param string:
    :param reg:
    :return:
    """

    start_chars = " #" if not string.startswith("#") else " "

    if reg.has_info: string += start_chars
    if reg.has_label: string += " text={" + reg.label + "}"
    if reg.has_meta:
        if "text" in reg.meta: string += " text={" + reg.meta["text"] + "}"
        string += " " + " ".join(key + "=" + value for key, value in reg.meta.items() if types.is_string_type(value) and key != "text")
    if reg.has_appearance: string += " " + " ".join(key + "=" + value for key, value in reg.appearance.items())
    return string

# -----------------------------------------------------------------

coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image', 'ecliptic']
coordinate_systems += ['wcs{0}'.format(letter) for letter in string.ascii_lowercase]

coordsys_name_mapping = dict(zip(frame_transform_graph.get_names(), frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO

# -----------------------------------------------------------------

class Region(object):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set the label
        self.label = kwargs.pop("label", None)

        # Set the 'include' flag
        self.include = kwargs.pop("include", True)

        # Set the appearance info
        self.appearance = kwargs.pop("appearance", dict())

        # Set the meta information
        self.meta = kwargs.pop("meta", dict())

    # -----------------------------------------------------------------

    @property
    def has_label(self):

        """
        This function ...
        :return:
        """

        return self.label is not None

    # -----------------------------------------------------------------

    @property
    def has_appearance(self):

        """
        This function ...
        :return:
        """

        return len(self.appearance) > 0

    # -----------------------------------------------------------------

    @property
    def has_meta(self):

        """
        This function ...
        :return:
        """

        return len(self.meta) > 0

    # -----------------------------------------------------------------

    @property
    def has_info(self):

        """
        This property ...
        :return:
        """

        return self.has_label or self.has_appearance or self.has_meta

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

class PixelRegion(Region):

    """
    This class ...
    """

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        copy = self.copy()
        copy += other
        return copy

    # -----------------------------------------------------------------

    def __iadd__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        self.center += other
        return self

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        copy = self.copy()
        copy -= other
        return copy

    # -----------------------------------------------------------------

    def __isub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        self.center -= other
        return self

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This property ...
        :return:
        """

        return self.axis1_min

    # -----------------------------------------------------------------

    @property
    def x_min_pixel(self):

        """
        This function ...
        :return:
        """

        return int(round(self.x_min))

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return self.axis1_max

    # -----------------------------------------------------------------

    @property
    def x_max_pixel(self):

        """
        This function ...
        :return:
        """

        return int(round(self.x_max))

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return self.axis2_min

    # -----------------------------------------------------------------

    @property
    def y_min_pixel(self):

        """
        This function ...
        :return:
        """

        return int(round(self.y_min))

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return self.axis2_max

    # -----------------------------------------------------------------

    @property
    def y_max_pixel(self):

        """
        This function ...
        :return:
        """

        return int(round(self.y_max))

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import PixelRectangleRegion

        # Create the rectangle region and return it
        return PixelRectangleRegion(self.center, self.unrotated_radius)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        coordsys = 'image'
        output = '# Region file format: DS9 PTS/magic/region\n'
        output += '{}\n'.format(coordsys)
        output += str(self)

        # Write
        with open(path, 'w') as fh: fh.write(output)

# -----------------------------------------------------------------

class SkyRegion(Region):

    """
    This class ...
    """

    @property
    def ra_min(self):
        return self.axis1_min

    # -----------------------------------------------------------------

    @property
    def ra_max(self):
        return self.axis1_max

    # -----------------------------------------------------------------

    @property
    def dec_min(self):
        return self.axis2_min

    # -----------------------------------------------------------------

    @property
    def dec_max(self):
        return self.axis2_max

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import SkyRectangleRegion

        # Create the rectangle region and return it
        return SkyRectangleRegion(self.center, self.unrotated_radius)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        coordsys = 'fk5'
        output = '# Region file format: DS9 PTS/magic/region\n'
        output += '{}\n'.format(coordsys)
        output += str(self)

        # Write
        with open(path, 'w') as fh: fh.write(output)

# -----------------------------------------------------------------

class PhysicalRegion(Region):

    """
    This class ...
    """

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Import the relevant PTS classes and modules
        from .rectangle import PhysicalRectangleRegion

        # Create the rectangle region and return it
        return PhysicalRectangleRegion(self.center, self.radius)

# -----------------------------------------------------------------
