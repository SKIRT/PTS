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
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from .point import PointRegion, PixelPointRegion, SkyPointRegion, PhysicalPointRegion
from .line import LineRegion, PixelLineRegion, SkyLineRegion, PhysicalLineRegion
from .vector import VectorRegion, PixelVectorRegion, SkyVectorRegion, PhysicalVectorRegion
from .circle import CircleRegion, PixelCircleRegion, SkyCircleRegion, PhysicalCircleRegion
from .ellipse import EllipseRegion, PixelEllipseRegion, SkyEllipseRegion, PhysicalEllipseRegion
from .rectangle import RectangleRegion, PixelRectangleRegion, SkyRectangleRegion, PhysicalRectangleRegion
from .polygon import PolygonRegion, PixelPolygonRegion, SkyPolygonRegion, PhysicalPolygonRegion
from .text import TextRegion, PixelTextRegion, SkyTextRegion, PhysicalTextRegion
from .composite import CompositeRegion, PixelCompositeRegion, SkyCompositeRegion, PhysicalCompositeRegion

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
length = parse_angular_length_quantity
width = parse_angular_length_quantity
height = parse_angular_length_quantity
angle = parse_angular_length_quantity

# For the sake of readability in describing the spec, parse_coordinate etc. are renamed here
coordinate = parse_coordinate
language_spec = {'point': (coordinate, coordinate),
                 'line': (coordinate, coordinate, coordinate, coordinate),
                 'vector': (coordinate, coordinate, length, angle),
                 'circle': (coordinate, coordinate, radius),
                 'ellipse': itertools.chain((coordinate, coordinate), itertools.cycle((radius,))), # This is a special case to deal with n elliptical annuli
                 'box': (coordinate, coordinate, width, height, angle),
                 'polygon': itertools.cycle((coordinate,)),
                 'text': (coordinate, coordinate),
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
        else: coord_list.append(element_parser(element))

    return coord_list

# -----------------------------------------------------------------

# match an x=y pair (where y can be any set of characters) that may or may not
# be followed by another one
meta_token = re.compile("([a-zA-Z]+)(=)([^= ]+) ?")

# -----------------------------------------------------------------

def meta_parser(meta_str):

    """
    Parse the metadata for a single ds9 region string.
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

    startswithhash = False
    if line.startswith("#"):
        startswithhash = True
        line = line[1:].strip()

    region_type_search = region_type_or_coordsys_re.search(line)

    #print(region_type_search)

    if region_type_search:
        include = region_type_search.groups()[0]
        region_type = region_type_search.groups()[1]
        #print("INCLUDE", include)
        #print("REGION TYPE", region_type)
    else: return

    if include == "": include = True
    elif include == "-": include = False
    else: raise RuntimeError("Encountered a problem (include = '" + include + "')")

    if region_type in coordinate_systems: return region_type  # outer loop has to do something with the coordinate system information
    elif region_type == "composite":

        #print(line)

        # end_of_region_name is the coordinate of the end of the region's name, e.g.:
        # circle would be 6 because circle is 6 characters
        end_of_region_name = region_type_search.span()[1]

        # coordinate of the # symbol or end of the line (-1) if not found
        if startswithhash: hash_or_end = line.find("||")
        else: hash_or_end = line.find("#")
        meta_str = line[hash_or_end:]

        #print(hash_or_end)
        #print(meta_str)

        parsed_meta = meta_parser(meta_str)

        #print(parsed_meta)

        return region_type, parsed_meta, include

    elif region_type in language_spec:

        if coordsys is None: raise ValueError("No coordinate system specified and a region has been found.")

        if "||" in line: composite = True
        else: composite = False

        # end_of_region_name is the coordinate of the end of the region's name, e.g.:
        # circle would be 6 because circle is 6 characters
        end_of_region_name = region_type_search.span()[1]

        # coordinate of the # symbol or end of the line (-1) if not found
        if startswithhash:
            hash_or_end = line.find(") ")
            coords_etc = line.split("(")[1].split(")")[0]
            meta_str = line[hash_or_end:]
        else:
            hash_or_end = line.find("#")
            coords_etc = strip_paren(line[end_of_region_name:hash_or_end].strip(" |"))
            meta_str = line[hash_or_end:]

        parsed_meta = meta_parser(meta_str)

        if coordsys in coordsys_name_mapping:

            #print(line)
            #print(coords_etc)
            #print(language_spec[region_type])
            #print(coordsys_name_mapping[coordsys])

            parsed = type_parser(coords_etc, language_spec[region_type], coordsys_name_mapping[coordsys])

            # Reset iterator for ellipse annulus
            if region_type == 'ellipse':
                language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            parsed_angles = [(x, y) for x, y in zip(parsed[:-1:2], parsed[1::2]) if isinstance(x, Angle) and isinstance(x, Angle)]
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
            if region_type == 'ellipse': language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            return region_type, parsed_return, parsed_meta, composite, include

# -----------------------------------------------------------------

def add_ds9_regions_from_string(region_string, regions, only=None, ignore=None, color=None, ignore_color=None):

    """
    Parse a DS9 region string.
    Parameters
    ----------
    region_string : str
        DS9 region string
    regions:
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
    composite_region = None
    composite_meta = None
    composite_include = None

    # ds9 regions can be split on \n or ;
    lines = []
    for line_ in region_string.split('\n'):
        for line in line_.split(";"):

            lines.append(line)

            parsed = line_parser(line, coordsys)

            if parsed is None: continue

            # If the line specifies the coordinate system for the region
            elif parsed in coordinate_systems: coordsys = parsed

            elif parsed[0] == "composite":

                composite_meta = parsed[1]
                composite_include = parsed[2]

            # Else: a line with one or more regions
            elif parsed:

                region_type, coordlist, meta, composite, include = parsed
                #meta['include'] = include
                #log.debug("Region type = {0}".format(region_type))

                # If the parsed region is part of a composite and it is the first
                if composite and composite_region is None: composite_region = [(region_type, coordlist, include)]

                # If the parsed region is part of a composite and it is not the first
                elif composite: composite_region.append((region_type, coordlist, include))

                # Finish composite
                elif composite_region is not None:

                    composite_region.append((region_type, coordlist, include))

                    # MAKE COMPOSITE REGION
                    specs = (composite_region, composite_meta, composite_include)
                    region = make_composite_region(specs)

                    composite_region = None
                    composite_meta = None
                    composite_include = None

                    regions.append(region)

                else:

                    # MAKE ORDINARY REGION

                    specs = (region_type, coordlist, meta, include)
                    region = make_regular_region(specs)

                    #regions.append((region_type, coordlist, meta))

                    regions.append(region)

# -----------------------------------------------------------------

def composite_to_string(composite, ds9_strings, frame, radunit, fmt):

    """
    This function ...
    :param composite:
    :param ds9_strings:
    :param frame:
    :param radunit:
    :param fmt:
    :return:
    """

    prefix = "" if composite.include else "-"

    if isinstance(composite, SkyCompositeRegion):

        # Get properties
        x = float(composite.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(composite.center.transform_to(frame).spherical.lat.to('deg').value)
        ang = composite.angle

        # Create string
        composite_string = prefix + ds9_strings['composite'].format(**locals())

    elif isinstance(composite, PixelCompositeRegion):

        # Get properties
        x = composite.center.x
        y = composite.center.y
        ang = composite.angle

        # Create string
        composite_string = prefix + ds9_strings['composite'].format(**locals())

    # Invalid value for the composite region
    else: raise ValueError("Invalid value for 'composite'")

    # Add the strings for the composite elements
    output = composite_string + "\n" + " ||\n".join([regular_to_string(element, ds9_strings, frame, radunit, fmt) for element in composite.elements])

    # Return the string
    return output

# -----------------------------------------------------------------

def regular_to_string(reg, ds9_strings, frame, radunit, fmt):

    """
    This function ...
    :param reg:
    :param ds9_strings:
    :param frame:
    :param radunit:
    :param fmt:
    :return:
    """

    if reg.include: prefix = ""
    else: prefix = "-"

    # Point region in sky coordinates
    if isinstance(reg, SkyPointRegion):

        x = float(reg.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.transform_to(frame).spherical.lat.to('deg').value)

        return prefix + ds9_strings['point'].format(**locals())

    # Point region in pixel coordinates
    elif isinstance(reg, PixelPointRegion):

        x = reg.x
        y = reg.y

        return prefix + ds9_strings['point'].format(**locals())

    elif isinstance(reg, SkyLineRegion):

        x1 = float(reg.start.transform_to(frame).spherical.lon.to('deg').value)
        y1 = float(reg.start.transform_to(frame).spherical.lat.to('deg').value)

        x2 = float(reg.end.transform_to(frame).spherical.lon.to('deg').value)
        y2 = float(reg.end.transform_to(frame).spherical.lat.to('deg').value)

        # TO HAVE IT IN hh:mm:ss NOTATION:

        #print(vars(reg.start.transform_to(frame)))
        skycoordinate = SkyCoordinate(reg.start.transform_to(frame).spherical.lon, reg.start.transform_to(frame).spherical.lat, frame=frame, representation="spherical")
        str1 = skycoordinate.to_string('hmsdms').replace("d", ":").replace("h", ":").replace("m", ":").replace("s ", ",")[:-1]

        skycoordinate = SkyCoordinate(reg.end.transform_to(frame).spherical.lon, reg.end.transform_to(frame).spherical.lat, frame=frame, representation="spherical")
        str2 = skycoordinate.to_string('hmsdms').replace("d", ":").replace("h", ":").replace("m", ":").replace("s ", ",")[:-1]

        return prefix + ds9_strings['line'].format(**locals())

    elif isinstance(reg, PixelLineRegion):

        x1 = reg.start.x
        y1 = reg.end.y

        x2 = reg.end.x
        y2 = reg.end.y

        return prefix + ds9_strings['line'].format(**locals())

    elif isinstance(reg, SkyVectorRegion):

        x = float(reg.start.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.start.transform_to(frame).spherical.lat.to('deg').value)
        l = float(reg.length.to(radunit).value)
        ang = float(reg.angle.to('deg').value)

        return prefix + ds9_strings['vector'].format(**locals())

    elif isinstance(reg, PixelVectorRegion):

        x = reg.start.x
        y = reg.start.y
        l = reg.length
        ang = reg.angle

        return prefix + ds9_strings['vector'].format(**locals())

    elif isinstance(reg, SkyCircleRegion):

        x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
        r = float(reg.radius.to(radunit).value)

        return prefix + ds9_strings['circle'].format(**locals())

    elif isinstance(reg, PixelCircleRegion):

        x = reg.center.x
        y = reg.center.y
        r = reg.radius

        return prefix + ds9_strings['circle'].format(**locals())

    elif isinstance(reg, SkyEllipseRegion):

        x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
        r1 = float(reg.major.to(radunit).value)
        r2 = float(reg.minor.to(radunit).value)
        ang = float(reg.angle.to('deg').value)

        return prefix + ds9_strings['ellipse'].format(**locals())

    elif isinstance(reg, PixelEllipseRegion):

        x = reg.center.x
        y = reg.center.y
        r1 = reg.major
        r2 = reg.minor
        ang = reg.angle

        return prefix + ds9_strings['ellipse'].format(**locals())

    elif isinstance(reg, SkyRectangleRegion):

        x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
        d1 = 2.0 * float(reg.radius.ra.to(radunit).value)
        d2 = 2.0 * float(reg.radius.dec.to(radunit).value)
        ang = float(reg.angle.to('deg').value)

        return prefix + ds9_strings['rectangle'].format(**locals())

    elif isinstance(reg, PixelRectangleRegion):

        x = reg.center.x
        y = reg.center.y
        d1 = 2.0 * reg.radius.x
        d2 = 2.0 * reg.radius.y
        ang = reg.angle

        return prefix + ds9_strings['rectangle'].format(**locals())

    elif isinstance(reg, SkyPolygonRegion):

        v = reg.points.transform_to(frame)
        coords = [(x.to('deg').value, y.to('deg').value) for x, y in
                  zip(v.spherical.lon, v.spherical.lat)]
        val = "{:" + fmt + "}"
        temp = [val.format(x) for _ in coords for x in _]
        c = ",".join(temp)

        return prefix + ds9_strings['polygon'].format(**locals())

    elif isinstance(reg, PixelPolygonRegion):

        v = reg.points
        coords = [(x, y) for x, y in zip(v.x, v.y)]
        val = "{:" + fmt + "}"
        temp = [val.format(x) for _ in coords for x in _]
        c = ",".join(temp)

        return prefix + ds9_strings['polygon'].format(**locals())

    elif isinstance(reg, SkyTextRegion):

        x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
        y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
        text = reg.text

        return prefix + ds9_strings['text'].format(**locals())

    elif isinstance(reg, PixelTextRegion):

        x = reg.center.x
        y = reg.center.y
        text = reg.text

        return prefix + ds9_strings['text'].format(**locals())

# -----------------------------------------------------------------

def ds9_objects_to_string(regions, coordsys='fk5', fmt='.4f', radunit='deg'):

    """
    Convert list of regions ato ds9 region strings.
    :param regions:
    :param coordsys:
    :param fmt:
    :param radunit:
    """

    if radunit == 'arcsec':
        if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
        else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
    else: radunitstr = ''

    ds9_strings = {
        'point': 'point({x:' + fmt + '},{y:' + fmt + '})',
        #'line': 'line({x1:' + fmt + '},{y1:' + fmt + '},{x2:' + fmt + '},{y2:' + fmt + '}) # line=0 0',
        'line': 'line({str1:},{str2:}) # line=0 0',
        'vector': '# vector({x:' + fmt + '},{y:' + fmt + '},{l:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '}) vector=1',
        'circle': 'circle({x:' + fmt + '},{y:' + fmt + '},{r:' + fmt + '}' + radunitstr + ')',
        'ellipse': 'ellipse({x:' + fmt + '},{y:' + fmt + '},{r1:' + fmt + '}' + radunitstr + ',{r2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})',
        'rectangle': 'box({x:' + fmt + '},{y:' + fmt + '},{d1:' + fmt + '}' + radunitstr + ',{d2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})',
        'polygon': 'polygon({c})',
        'text': '# text({x:' + fmt + '},{y:' + fmt + '}) text={text:}',
        'composite': '# composite({x:' + fmt + '},{y:' + fmt + '},{ang:' + fmt + '}) || composite=1'
    }

    output = '# Region file format: DS9 PTS/magic/region\n'
    output += '{}\n'.format(coordsys)

    # convert coordsys string to coordsys object
    if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
    else: frame = None # for pixel/image/physical frames

    # Loop over the regions
    for reg in regions:

        if isinstance(reg, CompositeRegion): output += composite_to_string(reg, ds9_strings, frame, radunit, fmt) + "\n"
        else: output += regular_to_string(reg, ds9_strings, frame, radunit, fmt) + "\n"

    # Return the output string
    return output

# -----------------------------------------------------------------

def make_composite_region(specs):

    """
    This function ...
    :param specs:
    :return:
    """

    regions = []

    composite_specs, meta, include = specs

    appearance = {key: meta[key] for key in meta.keys() if key in viz_keywords}
    meta = {key: meta[key] for key in meta.keys() if key not in viz_keywords}

    # Create regions from the specs
    for spec in composite_specs:

        spec = (spec[0], spec[1], {}, spec[2])
        reg = make_regular_region(spec)
        regions.append(reg)

    # Pixel region
    if isinstance(regions[0], PixelRegion):

        # Check if other are also pixel regions: is done by PixelCompositeRegion class

        # Add all pixel regions as a composite
        region = PixelCompositeRegion(*regions, meta=meta, appearance=appearance, include=include)

    elif isinstance(regions[0], SkyRegion):

        # Add all sky regions as a composite
        region = SkyCompositeRegion(*regions, meta=meta, appearance=appearance, include=include)

    # Hmm
    else: raise ValueError("Something went wrong: encountered" + repr(regions[0]) + " of type " + str(type(regions[0])))

    # Return the region
    return region

# -----------------------------------------------------------------

viz_keywords = ['color', 'dashed', 'width', 'point', 'font', 'text']

# -----------------------------------------------------------------

def make_regular_region(specs):

    """
    This function ...
    :param specs:
    :return:
    """

    region_type, coord_list, meta, include = specs

    appearance = {key: meta[key] for key in meta.keys() if key in viz_keywords}
    meta = {key: meta[key] for key in meta.keys() if key not in viz_keywords}

    # POINTS
    if region_type == 'point':

        if isinstance(coord_list[0], BaseCoordinateFrame):

            ra = coord_list[0].ra
            dec = coord_list[0].dec

            # Create the region
            reg = SkyPointRegion(ra, dec, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            # Create the region
            reg = PixelPointRegion(coord_list[0].x, coord_list[0].y, meta=meta, appearance=appearance, include=include)

        else: raise ValueError("No central coordinate")

    # LINES
    elif region_type == "line":

        if isinstance(coord_list[0], BaseCoordinateFrame):

            coord = coord_list[0]

            coord_1 = SkyCoordinate(coord[0].ra, coord[0].dec)
            coord_2 = SkyCoordinate(coord[1].ra, coord[1].dec)

            # Create the line
            reg = SkyLineRegion(coord_1, coord_2, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            #coord = coord_list[]

            #print(coord_list)

            coord_1 = coord_list[0]
            coord_2 = coord_list[1]

            # Create the line
            reg = PixelLineRegion(coord_1, coord_2, meta=meta, appearance=appearance, include=include)

    # VECTORS
    elif region_type == "vector":

        # Sky coordinates
        if isinstance(coord_list[0], BaseCoordinateFrame):

            # Get the start coordinate
            coord = coord_list[0]
            start = SkyCoordinate(coord.ra, coord.dec)

            # Get the length
            length = coord_list[1]

            # Get the angle
            angle = coord_list[2]
            angle = Angle(angle.to("deg"), "deg")

            # Create the vector
            reg = SkyVectorRegion(start, length, angle, meta=meta, appearance=appearance, include=include)

        # Pixel coordinates
        elif isinstance(coord_list[0], PixelCoordinate):

            # Get the start coordinate
            start = coord_list[0]

            # Get the length
            length = coord_list[1]

            # Get the angle
            angle = coord_list[2]
            angle = Angle(angle.to("deg"), "deg")

            # Create the vector
            reg = PixelVectorRegion(start, length, angle, meta=meta, appearance=appearance, include=include)

    # CIRCLES
    elif region_type == 'circle':

        if isinstance(coord_list[0], BaseCoordinateFrame):

            ra = coord_list[0].ra
            dec = coord_list[0].dec

            # Get the center coordinate
            center = SkyCoordinate(ra, dec)

            # Get the radius
            radius = coord_list[1]

            # Create a circle
            reg = SkyCircleRegion(center, radius, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            center = coord_list[0]
            radius = coord_list[1]

            # Create the circle
            reg = PixelCircleRegion(center, radius, meta=meta, appearance=appearance, include=include)

        # No central coordinate for this circle
        else: raise ValueError("No central coordinate")

    # ELLIPSES
    elif region_type == 'ellipse':

        # Do not read elliptical annuli for now
        #if len(coord_list) > 4: continue
        if isinstance(coord_list[0], BaseCoordinateFrame):

            ra = coord_list[0].ra
            dec = coord_list[0].dec

            # Get the center coordinate
            center = SkyCoordinate(ra, dec)

            # Get the radius
            radius = SkyStretch(coord_list[1], coord_list[2])

            # Get the angle
            angle = coord_list[3]
            angle = Angle(angle.to("deg"), "deg")

            # Create an ellipse
            reg = SkyEllipseRegion(center, radius, angle, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            center = coord_list[0]

            radius = PixelStretch(coord_list[1], coord_list[2])

            # Get the angle
            angle = coord_list[3]
            angle = Angle(angle.to("deg"), "deg")

            # Create the region
            reg = PixelEllipseRegion(center, radius, angle, meta=meta, appearance=appearance, include=include)

        # No central coordinate found for this ellipse
        else: raise ValueError("No central coordinate")

    # RECTANGLES
    elif region_type == 'rectangle':

        if isinstance(coord_list[0], BaseCoordinateFrame):

            ra = coord_list[0].ra
            dec = coord_list[0].dec

            # Get the center coordinate
            center = SkyCoordinate(ra, dec)

            # Get the radius
            radius = SkyStretch(0.5 * coord_list[1], 0.5 * coord_list[2])

            # Get the angle
            angle = coord_list[3]
            angle = Angle(angle.to("deg"), "deg")

            # Create the region
            reg = SkyRectangleRegion(center, radius, angle, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            center = coord_list[0]

            radius = PixelStretch(0.5 * coord_list[1], 0.5 * coord_list[2])

            # Get the angle
            angle = coord_list[3]
            angle = Angle(angle.to("deg"), "deg")

            # Create the region
            reg = PixelRectangleRegion(center, radius, angle, meta=meta, appearance=appearance, include=include)

        # No central coordinate given
        else: raise ValueError("No central coordinate")

    # POLYGONS
    elif region_type == 'polygon':

        if isinstance(coord_list[0], BaseCoordinateFrame):

            reg = SkyPolygonRegion(coord_list[0], meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            reg = PixelPolygonRegion(coord_list[0], meta=meta, appearance=appearance, include=include)

        # No central coordinate given
        else: raise ValueError("No central coordinate")

    # TEXT
    elif region_type == "text":

        text = appearance.pop("text")

        if isinstance(coord_list[0], BaseCoordinateFrame):

            ra = coord_list[0].ra
            dec = coord_list[0].dec

            reg = SkyTextRegion(SkyCoordinate(ra, dec), text, meta=meta, appearance=appearance, include=include)

        elif isinstance(coord_list[0], PixelCoordinate):

            reg = PixelTextRegion(coord_list[0], text, meta=meta, appearance=appearance, include=include)

        # No central coordinate given
        else: raise ValueError("No central coordinate")

    # UNKNOWN
    else: raise ValueError("Region specification not recognized")

    # Return the region
    return reg

# -----------------------------------------------------------------

class RegionList(list):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, region_string, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param region_string:
        :return:
        """

        # Create the region list
        regions = cls()

        # Get list
        add_ds9_regions_from_string(region_string, regions, only, ignore, color, ignore_color)

        # Return the regions
        return regions

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Get the region string
        with open(path) as fh: region_string = fh.read()

        # Open from string
        return cls.from_string(region_string, only, ignore, color, ignore_color)

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class (list)
        super(RegionList, self).__init__()

    # -----------------------------------------------------------------

    def save(self, path, coordsys='fk5'):

        """
        This function ...
        :param path:
        :param coordsys:
        :return:
        """

        # Convert to string
        output = ds9_objects_to_string(self, coordsys)

        # Write
        with open(path, 'w') as fh: fh.write(output)

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

    def pixel_regions(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, PixelRegion)]

    # -----------------------------------------------------------------

    def sky_regions(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, SkyRegion)]

    # -----------------------------------------------------------------

    def physical_regions(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, PhysicalRegion)]

    # -----------------------------------------------------------------

    @property
    def only_pixel(self):

        """
        This function ...
        :return:
        """

        return len(self.sky_regions()) == 0 and len(self.physical_regions()) == 0

    # -----------------------------------------------------------------

    @property
    def only_sky(self):

        """
        This function ...
        :return:
        """

        return len(self.pixel_regions()) == 0 and len(self.physical_regions()) == 0

    # -----------------------------------------------------------------

    @property
    def only_physical(self):

        """
        This function ...
        :return:
        """

        return len(self.pixel_regions()) == 0 and len(self.sky_regions()) == 0

    # -----------------------------------------------------------------

    def points(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, PointRegion)]

    # -----------------------------------------------------------------

    def lines(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, LineRegion)]

    # -----------------------------------------------------------------

    def vectors(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, VectorRegion)]

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

    def polygons(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, PolygonRegion)]

    # -----------------------------------------------------------------

    def texts(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, TextRegion)]

    # -----------------------------------------------------------------

    def composites(self):

        """
        This function ...
        :return:
        """

        return [region for region in self if isinstance(region, CompositeRegion)]

# -----------------------------------------------------------------

class PixelRegionList(RegionList):

    """
    This class ...
    """

    @classmethod
    def from_string(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(PixelRegionList, cls).from_string(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_pixel: raise ValueError("String does not only contain pixel regions")

        # Return
        return regions

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(PixelRegionList, cls).from_file(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_pixel: raise ValueError("String does not only contain pixel regions")

        # Return
        return regions

    # -----------------------------------------------------------------

    def append(self, region):

        """
        This function ...
        :param region:
        :return:
        """

        # Check whether the argument is a valid shape
        if not isinstance(region, PixelRegion): raise ValueError("Region is not a pixel region")

        # Otherwise, add the shape
        super(PixelRegionList, self).append(region)

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

    @classmethod
    def from_string(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(SkyRegionList, cls).from_string(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_sky: raise ValueError("String does not only contain sky regions")

        # Return the region list
        return regions

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(SkyRegionList, cls).from_file(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_sky: raise ValueError("String does not only contain sky regions")

        # Return the region list
        return regions

# -----------------------------------------------------------------

class PhysicalRegionList(RegionList):

    """
    This class ...
    """

    @classmethod
    def from_string(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(PhysicalRegionList, cls).from_file(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_physical: raise ValueError("String does not only contain physical regions")

        # Return the region list
        return regions

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :return:
        """

        # Call the function of the base class
        regions = super(PhysicalRegionList, cls).from_file(path, only, ignore, color, ignore_color)

        # Check
        if not regions.only_physical: raise ValueError("String does not only contain physical regions")

        # Return the region list
        return regions

# -----------------------------------------------------------------
