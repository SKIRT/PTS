#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.list Contains the RegionList class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
import string
import itertools

# Import astronomical modules
from astropy.coordinates import Angle, frame_transform_graph, UnitSphericalRepresentation, BaseCoordinateFrame
from astropy.units import Unit, Quantity, dimensionless_unscaled

# Import the relevant PTS classes and modules
from ..basics.vector import Extent
from ..basics.mask import Mask

from ..basics.coordinate import Coordinate, PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from .point import PixelPointRegion, SkyPointRegion, PhysicalPointRegion
from .circle import PixelCircleRegion, SkyCircleRegion, PhysicalCircleRegion
from .ellipse import PixelEllipseRegion, SkyEllipseRegion, PhysicalEllipseRegion
from .rectangle import PixelRectangleRegion, SkyRectangleRegion, PhysicalRectangleRegion

# -----------------------------------------------------------------

region_type_or_coordsys_re = re.compile("^#? *(-?)([a-zA-Z0-9]+)")

coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image', 'ecliptic']
coordinate_systems += ['wcs{0}'.format(letter) for letter in string.ascii_lowercase]

coordsys_name_mapping = dict(zip(frame_transform_graph.get_names(), frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO

# -----------------------------------------------------------------

hour_or_deg = 'hour_or_deg'
coordinate_units = {'fk5': (hour_or_deg, Unit("deg")),
                    'fk4': (hour_or_deg, Unit("deg")),
                    'icrs': (hour_or_deg, Unit("deg")),
                    'geocentrictrueecliptic': (Unit("deg"), Unit("deg")),
                    'galactic': (Unit("deg"), Unit("deg")),
                    'physical': (dimensionless_unscaled, dimensionless_unscaled),
                    'image': (dimensionless_unscaled, dimensionless_unscaled),
                    'wcs': (dimensionless_unscaled, dimensionless_unscaled),
                    }
for letter in string.ascii_lowercase:
    coordinate_units['wcs{0}'.format(letter)] = (dimensionless_unscaled, dimensionless_unscaled)

# -----------------------------------------------------------------

unit_mapping = {
    '"': Unit("arcsec"),
    "'": Unit("arcmin"),
    'r': Unit("rad"),
    'i': dimensionless_unscaled,
}

# -----------------------------------------------------------------

paren = re.compile("[()]")

def strip_paren(string_rep):
    return paren.sub("", string_rep)

# -----------------------------------------------------------------

def parse_coordinate(string_rep, unit):

    """
    Parse a single coordinate
    """
    # Any ds9 coordinate representation (sexagesimal or degrees)

    if 'd' in string_rep or 'h' in string_rep:
        return Angle(string_rep)
    elif unit is 'hour_or_deg':
        if ':' in string_rep:
            spl = tuple([float(x) for x in string_rep.split(":")])
            return Angle(spl, "hourangle")
        else:
            ang = float(string_rep)
            return Angle(ang, "deg")
    elif unit.is_equivalent("deg"):
        # return coordinates.Angle(string_rep, unit=unit)
        if ':' in string_rep:
            ang = tuple([float(x) for x in string_rep.split(":")])
        else:
            ang = float(string_rep)
        return Angle(ang, "deg")
    else: return Quantity(float(string_rep), unit)

# -----------------------------------------------------------------

def parse_angular_length_quantity(string_rep):

    """
    Given a string that is either a number or a number and a unit, return a
    Quantity of that string.  e.g.:
        23.9 -> 23.9*u.deg
        50" -> 50*u.arcsec
    """
    has_unit = string_rep[-1] not in string.digits
    if has_unit:
        unit = unit_mapping[string_rep[-1]]
        return Quantity(float(string_rep[:-1]), unit=unit)
    else: return Quantity(float(string_rep), unit="deg")

# -----------------------------------------------------------------

# these are the same function, just different names
radius = parse_angular_length_quantity
width = parse_angular_length_quantity
height = parse_angular_length_quantity
angle = parse_angular_length_quantity

# For the sake of readability in describing the spec, parse_coordinate etc. are renamed here
coordinate = parse_coordinate
language_spec = {'point': (coordinate, coordinate),
                 'circle': (coordinate, coordinate, radius),
                 # This is a special case to deal with n elliptical annuli
                 'ellipse': itertools.chain((coordinate, coordinate), itertools.cycle((radius,))),
                 'box': (coordinate, coordinate, width, height, angle),
                 'polygon': itertools.cycle((coordinate,)),
                 }

# -----------------------------------------------------------------

def type_parser(string_rep, specification, coordsys):

    """
    For a given region line in which the type has already been determined,
    parse the coordinate definition
    Parameters
    ----------
    string_rep : str
        The string containing the coordinates.  For example, if your region is
        `circle(1,2,3)` this string would be `(1,2,3)`
    specification : iterable
        An iterable of coordinate specifications.  For example, for a circle,
        this would be a list of (coordinate, coordinate, radius).  Each
        individual specification should be a function that takes a string and
        returns the appropriate astropy object.  See ``language_spec`` for the
        definition of the grammar used here.
    coordsys : str
        The string name of the global coordinate system
    Returns
    -------
    coord_list : list
        The list of astropy coordinates and/or quantities representing radius,
        width, etc. for the region
    """
    coord_list = []
    splitter = re.compile("[, ]")
    for ii, (element, element_parser) in enumerate(zip(splitter.split(string_rep), specification)):
        if element_parser is coordinate:
            unit = coordinate_units[coordsys][ii % 2]
            coord_list.append(element_parser(element, unit))
        else:
            coord_list.append(element_parser(element))

    return coord_list

# -----------------------------------------------------------------

# match an x=y pair (where y can be any set of characters) that may or may not
# be followed by another one
meta_token = re.compile("([a-zA-Z]+)(=)([^= ]+) ?")

def meta_parser(meta_str):

    """Parse the metadata for a single ds9 region string.
    The metadata is everything after the close-paren of the region coordinate specification.
    All metadata is specified as key=value pairs separated by whitespace, but
    sometimes the values can also be whitespace separated.
    """
    meta_token_split = [x for x in meta_token.split(meta_str.strip()) if x]
    equals_inds = [i for i, x in enumerate(meta_token_split) if x is '=']
    result = {meta_token_split[ii - 1]:
                  " ".join(meta_token_split[ii + 1:jj - 1 if jj is not None else None])
              for ii, jj in zip(equals_inds, equals_inds[1:] + [None])}

    return result

# -----------------------------------------------------------------

def line_parser(line, coordsys=None):

    """
    Parse a single ds9 region line into a string
    Parameters
    ----------
    line : str
        A single ds9 region contained in a string
    coordsys : str
        The global coordinate system name declared at the top of the ds9 file
    Returns
    -------
    (region_type, parsed_return, parsed_meta, composite, include)
    region_type : str
    coord_list : list of coordinate objects
    meta : metadata dict
    composite : bool
        indicates whether region is a composite region
    include : bool
        Whether the region is included (False -> excluded)
    """
    region_type_search = region_type_or_coordsys_re.search(line)
    if region_type_search:
        include = region_type_search.groups()[0]
        region_type = region_type_search.groups()[1]
    else:
        return

    if region_type in coordinate_systems:
        return region_type  # outer loop has to do something with the coordinate system information
    elif region_type in language_spec:
        if coordsys is None:
            raise ValueError("No coordinate system specified and a region has been found.")

        if "||" in line: composite = True
        else: composite = False

        # end_of_region_name is the coordinate of the end of the region's name, e.g.:
        # circle would be 6 because circle is 6 characters
        end_of_region_name = region_type_search.span()[1]
        # coordinate of the # symbol or end of the line (-1) if not found
        hash_or_end = line.find("#")
        coords_etc = strip_paren(line[end_of_region_name:hash_or_end].strip(" |"))
        meta_str = line[hash_or_end:]

        parsed_meta = meta_parser(meta_str)

        if coordsys in coordsys_name_mapping:
            parsed = type_parser(coords_etc, language_spec[region_type],
                                 coordsys_name_mapping[coordsys])

            # Reset iterator for ellipse annulus
            if region_type == 'ellipse':
                language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            parsed_angles = [
                (x, y) for x, y in zip(parsed[:-1:2], parsed[1::2])
                if isinstance(x, Angle) and isinstance(x, Angle)
                ]
            frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])

            lon, lat = zip(*parsed_angles)
            lon, lat = Quantity(lon), Quantity(lat)
            sphcoords = UnitSphericalRepresentation(lon, lat)
            coords = frame(sphcoords)

            return region_type, [coords] + parsed[len(coords) * 2:], parsed_meta, composite, include

        else:

            parsed = type_parser(coords_etc, language_spec[region_type], coordsys)
            if region_type == 'polygon':
                # have to special-case polygon in the phys coord case b/c can't typecheck when iterating as in sky coord case
                #coord = PixCoord(parsed[0::2], parsed[1::2])
                coord = PixelCoordinate(parsed[0::2], parsed[1::2])
                parsed_return = [coord]
            else:
                parsed = [_.value for _ in parsed]
                #coord = PixCoord(parsed[0], parsed[1])
                coord = PixelCoordinate(parsed[0], parsed[1])
                parsed_return = [coord] + parsed[2:]

            # Reset iterator for ellipse annulus
            if region_type == 'ellipse':
                language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            return region_type, parsed_return, parsed_meta, composite, include

# -----------------------------------------------------------------

def ds9_string_to_region_list(region_string):

    """Parse a DS9 region string.
    Parameters
    ----------
    region_string : str
        DS9 region string
    Returns
    -------
    list of (region type, coord_list, meta, composite, include) tuples
    region_type : str
    coord_list : list of coordinate objects
    meta : metadata dict
    composite : bool
        indicates whether region is a composite region
    include : bool
        Whether the region is included (False -> excluded)
    """
    coordsys = None
    regions = []
    composite_region = None

    # ds9 regions can be split on \n or ;
    lines = []
    for line_ in region_string.split('\n'):
        for line in line_.split(";"):
            lines.append(line)
            parsed = line_parser(line, coordsys)
            if parsed in coordinate_systems:
                coordsys = parsed
            elif parsed:
                region_type, coordlist, meta, composite, include = parsed
                meta['include'] = include
                #log.debug("Region type = {0}".format(region_type))
                if composite and composite_region is None:
                    composite_region = [(region_type, coordlist)]
                elif composite:
                    composite_region.append((region_type, coordlist))
                elif composite_region is not None:
                    composite_region.append((region_type, coordlist))
                    regions.append(composite_region)
                    composite_region = None
                else:
                    regions.append((region_type, coordlist, meta))

    return regions

# -----------------------------------------------------------------

def ds9_region_list_to_objects(region_list):

    """
    Given a list of parsed region tuples, product a list of astropy objects.
    TODO: show example what a "region list" is.
    Parameters
    ----------
    region_list : list
        List of TODO???
    Returns
    -------
    regions : list
        List of `regions.Region` objects
    """

    viz_keywords = ['color', 'dashed', 'width', 'point', 'font', 'text']

    regions = RegionList()

    for region_type, coord_list, meta in region_list:

        # TODO: refactor, possible on the basis of # of parameters + sometimes handle corner cases

        # CIRCLES
        if region_type == 'circle':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                #reg = circle.CircleSkyRegion(coord_list[0], coord_list[1])
            #elif isinstance(coord_list[0], PixCoord):
            elif isinstance(coord_list[0], PixelCoordinate):
                #reg = circle.CirclePixelRegion(coord_list[0], coord_list[1])
            else: raise ValueError("No central coordinate")

        # ELLIPSES
        elif region_type == 'ellipse':
            # Do not read elliptical annuli for now
            if len(coord_list) > 4:
                continue
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = ellipse.EllipseSkyRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            #elif isinstance(coord_list[0], PixCoord):
            elif isinstance(coord_list[0], PixelCoordinate):
                reg = ellipse.EllipsePixelRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            else: raise ValueError("No central coordinate")

        # POLYGONS
        elif region_type == 'polygon':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = polygon.PolygonSkyRegion(coord_list[0])
            #elif isinstance(coord_list[0], PixCoord):
            elif isinstance(coord_list[0], PixelCoordinate):
                reg = polygon.PolygonPixelRegion(coord_list[0])
            else: raise ValueError("No central coordinate")

        # RECTANGLES
        elif region_type == 'rectangle':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = rectangle.RectangleSkyRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            #elif isinstance(coord_list[0], PixCoord):
            elif isinstance(coord_list[0], PixelCoordinate):
                reg = rectangle.RectanglePixelRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            else: raise ValueError("No central coordinate")

        # POINTS
        elif region_type == 'point':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                #reg = point.PointSkyRegion(coord_list[0])
            #elif isinstance(coord_list[0], PixCoord):
            elif isinstance(coord_list[0], PixelCoordinate):
                #reg = point.PointPixelRegion(coord_list[0])
                reg = coord_list[0]
            else: raise ValueError("No central coordinate")
        else: continue

        reg.vizmeta = {key: meta[key] for key in meta.keys() if key in viz_keywords}
        reg.meta = {key: meta[key] for key in meta.keys() if key not in viz_keywords}
        #output_list.append(reg)
        regions.append(reg)

    # Return the region list
    return regions

# -----------------------------------------------------------------

class RegionList(list):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, region_string):

        """
        This function ...
        :param region_string:
        :return:
        """

        # Get list
        region_list = ds9_string_to_region_list(region_string)

        # Get region objects
        regions = ds9_region_list_to_objects(region_list)

        # Return
        return regions

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the region string
        with open(path) as fh: region_string = fh.read()

        # Open from string
        return cls.from_string(region_string)

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class (list)
        super(RegionList, self).__init__()

    # -----------------------------------------------------------------

    def append(self, region):

        """
        This function ...
        :param region:
        :return:
        """

        # Check if the region is indeed a derived 'Region' object
        if not isinstance(region, Region): raise ValueError("Not a region object")

        # Otherwise, add the region
        super(RegionList, self).append(region)

    # -----------------------------------------------------------------

    def points(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, PointRegion)]

    # -----------------------------------------------------------------

    def circles(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, CircleRegion)]

    # -----------------------------------------------------------------

    def ellipses(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, EllipseRegion)]

    # -----------------------------------------------------------------

    def rectangles(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, RectangleRegion)]

# -----------------------------------------------------------------

class PixelRegionList(RegionList):

    """
    This class ...
    """

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :param only:
        :param ignore:
        :param color:
        :param ignore_color:
        :return:
        """

        # Create a new region
        #region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        #try:
        #    _region = pyregion.open(path)
        #    if not _region.check_imagecoord(): raise IOError("Region is not in image coordinates")
        #except ValueError: # If a ValueError comes out, assume the region file is empty (no shapes)
        #    _region = []

        #region_list = ds9_string_to_region_list(region_string)
        #regions = ds9_region_list_to_objects(region_list)
        #return regions

        with open(path) as fh:
            region_string = fh.read()

        #return ds9_string_to_objects(region_string)
        # Implementation of ds9_string_to_objects
        region_list = ds9_string_to_region_list(region_string)
        regions = ds9_region_list_to_objects(region_list)
        return regions

        # Loop over all shapes in the region
        for shape in _region:

            # Meta information
            meta = {}
            if "text" in shape.attr[1]: meta["text"] = shape.attr[1]["text"]
            if "color" in shape.attr[1]: meta["color"] = shape.attr[1]["color"]
            if "point" in shape.attr[1]: meta["point"] = shape.attr[1]["point"]

            # Check the color of the shape
            if color is not None and shape.attr[1]["color"] != color: continue
            if ignore_color is not None and shape.attr[1]["color"] == ignore_color: continue

            # If the shape is a point -> Position
            if shape.name == "point":

                if only is not None and "point" not in only: continue
                if ignore is not None and "point" in ignore: continue

                # Get the position
                x = shape.coord_list[0]
                y = shape.coord_list[1]
                coordinate = Coordinate(x, y, meta=meta)

                # Add the position to the region
                new_shape = coordinate

            # If the shape is a line -> Line
            elif shape.name == "line" or shape.name == "vector":

                if only is not None and "line" not in only: continue
                if ignore is not None and "line" in ignore: continue

                # Get the position of the two points
                x_1 = shape.coord_list[0]
                y_1 = shape.coord_list[1]
                x_2 = shape.coord_list[2]
                y_2 = shape.coord_list[3]

                # Create the positions
                position_1 = Coordinate(x_1, y_1)
                position_2 = Coordinate(x_2, y_2)

                # Create the line
                line = Line(position_1, position_2, meta=meta)
                new_shape = line

            # If the shape is a circle -> Circle
            elif shape.name == "circle":

                if only is not None and "circle" not in only: continue
                if ignore is not None and "circle" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the radius
                radius = shape.coord_list[2]

                # Create a circle
                circle = Circle(center, radius, meta=meta)
                new_shape = circle

            # If the shape is an ellipse -> Ellipse
            elif shape.name == "ellipse":

                if only is not None and "ellipse" not in only: continue
                if ignore is not None and "ellipse" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the radius
                x_radius = shape.coord_list[2]
                y_radius = shape.coord_list[3]
                radius = Extent(x_radius, y_radius)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create an ellipse
                ellipse = Ellipse(center, radius, angle, meta=meta)
                new_shape = ellipse

            # If the shape is a rectangle -> Rectangle
            elif shape.name == "box":

                if only is not None and "box" not in only: continue
                if ignore is not None and "box" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the width and height
                width = shape.coord_list[2]
                height = shape.coord_list[3]

                # Create radius
                radius = Extent(0.5 * width, 0.5 * height)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create a Rectangle
                rectangle = Rectangle(center, radius, angle, meta=meta)
                new_shape = rectangle

            # If the shape is a polygon -> Polygon
            elif shape.name == "polygon":

                if only is not None and "polygon" not in only: continue
                if ignore is not None and "polygon" in ignore: continue

                # Get the number of points in the polygon
                number_of_points = 0.5 * len(shape.coord_list)
                assert int(number_of_points) == number_of_points
                number_of_points = int(number_of_points)

                # Create a new Polygon
                polygon = Polygon(meta=meta)

                # Get the position of the different points
                for i in range(number_of_points):

                    # Create a new Position
                    x = shape.coord_list[2*i]
                    y = shape.coord_list[2*i + 1]
                    position = Coordinate(x, y)

                    # Add the point to the polygon
                    polygon.add_point(position)

                # Add the polygon to the region
                new_shape = polygon

            # Unrecognized shape
            else: raise ValueError("Unrecognized shape: " + shape.name)

            # If this shape should be excluded
            if shape.exclude:

                previous_shape = region[len(region)-1]
                region[len(region)-1] = Composite(previous_shape, new_shape, meta=previous_shape.meta)

            # Add the shape to the region
            else: region.append(new_shape)

        # Return the new region
        return region

    # -----------------------------------------------------------------

    def append(self, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        # Check whether the argument is a valid shape
        if not (shape.__class__.__name__ == "Coordinate" or shape.__class__.__name__ == "Line"
                or shape.__class__.__name__ == "Circle" or shape.__class__.__name__ == "Ellipse"
                or shape.__class__.__name__ == "Rectangle" or shape.__class__.__name__ == "Polygon"
                or shape.__class__.__name__ == "Composite"):
            raise ValueError("Shape must of of type Coordinate, Line, Circle, Ellipse, Rectangle, Polygon or Composite")

        # Otherwise, add the shape
        super(Region, self).append(shape)

    # -----------------------------------------------------------------

    def cropped(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # Create a new region
        region = Region()

        # Loop over all shapes
        for shape in self:

            # Ignore shapes that fall outside of the crop region
            if shape.x < x_min or shape.x > x_max or shape.y < y_min or shape.y > y_max: continue

            # Create a new shape that has its coordinates shifted
            new_shape = shape - Extent(x_min, y_min)

            # Add the new shape to the new region
            region.append(new_shape)

        # Return the new region
        return region

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = Region()
        for shape in self: new_region.append(shape * value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = Region()
        for shape in self: new_region.append(shape / value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Create a new SkyRegion
        region = SkyRegionList()

        # Add the shapes to the sky region
        for shape in self: region.append(shape.to_sky(wcs))

        # Return the region
        return region

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        # Create empty mask
        mask = Mask.empty(x_size, y_size)

        # Add each shape to the mask
        for shape in self: mask += Mask.from_shape(shape, x_size, y_size)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def to_mpl_patches(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the patches
        patches = []

        # Add the patches from the shapes
        for shape in self: patches.append(shape.to_mpl_patch())

        # Return the list of patches
        return patches

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Print the coordinate system
        print("image", file=f)

        # Loop over all shapes, get string and print it to the region file
        for shape in self: print(shape.to_region_string(coordinate_system=False), file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def homogenized(self):

        """
        This function returns a copy of the region where Composite objects have been dissolved into their individual components
        :return:
        """

        import copy

        new = Region()

        for shape in self:

            if isinstance(shape, Composite):

                copy_base = copy.deepcopy(shape.base)
                copy_exclude = copy.deepcopy(shape.exclude)
                copy_base.meta = shape.meta
                copy_exclude.meta = shape.meta

                new.append(copy_base)
                new.append(copy_exclude)

            else: new.append(copy.deepcopy(shape))

        # Return the new region
        return new

# -----------------------------------------------------------------

class SkyRegionList(RegionList):

    """
    This class ...
    """



# -----------------------------------------------------------------
