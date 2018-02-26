#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.extended Contains the ExtendedSourceCatalog class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ..basics.coordinate import SkyCoordinate
from ...core.basics.log import log
from ..core.extendedsource import ExtendedSource
from ..region.list import SkyRegionList
from ..region.list import PixelRegionList
from ..region.ellipse import SkyEllipseRegion
from ..basics.stretch import SkyStretch
from ..region.point import SkyPointRegion
from ...core.units.parsing import parse_quantity

# -----------------------------------------------------------------

class ExtendedSourceCatalog(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Name"] = (str, None, "Name of the galaxy")
    _column_info["RA"] = (float, "deg", "Right ascension")
    _column_info["DEC"] = (float, "deg", "Declination")
    _column_info["Redshift"] = (float, None, "Redshift")
    _column_info["Type"] = (str, None, "Galaxy type")
    _column_info["Names"] = (str, None, "Alternative names")
    _column_info["Distance"] = (float, "Mpc", "distance")
    _column_info["Incl"] = (float, "deg", "inclination")
    _column_info["D25"] = (float, "arcmin", "D25")
    _column_info["Major"] = (float, "arcmin", "Major axis length")
    _column_info["Minor"] = (float, "arcmin", "Minor axis length")
    _column_info["Posangle"] = (float, "deg", "Position angle")
    _column_info["Principal"] = (bool, None, "Is principal galaxy")
    _column_info["Companions"] = (str, None, "Companion galaxies")
    _column_info["Parent"] = (str, None, "Parent galaxy (is companion galaxy)")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExtendedSourceCatalog, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_entry(self, name, ra, dec, z, galtype, alternative_names, distance, inclination, d25, major, minor,
                  posangle, principal, companions, parent):

        """
        This function ...
        :param name:
        :param ra:
        :param dec:
        :param z:
        :param galtype:
        :param alternative_names:
        :param distance:
        :param inclination:
        :param d25:
        :param major:
        :param minor:
        :param posangle:
        :param principal:
        :param companions:
        :param parent:
        :return:
        """

        values = [name, ra, dec, z, galtype, alternative_names, distance, inclination, d25, major, minor, posangle,
                  principal, companions, parent]
        self.add_row(values)

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index], dec=self["DEC"][index], unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def get_name(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Name", index)

    # -----------------------------------------------------------------

    def get_redshift(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Redshift", index)

    # -----------------------------------------------------------------

    def get_type(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Type", index)

    # -----------------------------------------------------------------

    def get_names(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        names = self.get_quantity("Names", index)
        if names is not None: names = names.split(",")
        return names

    # -----------------------------------------------------------------

    def get_distance(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Distance", index)

    # -----------------------------------------------------------------

    def get_inclination(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Incl", index)

    # -----------------------------------------------------------------

    def get_d25(self, index):

        """
        This function ...
        :return:
        """

        return self.get_quantity("D25", index)

    # -----------------------------------------------------------------

    def get_major(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Major", index)

    # -----------------------------------------------------------------

    def get_minor(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Minor", index)

    # -----------------------------------------------------------------

    def get_position_angle(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Posangle", index)

    # -----------------------------------------------------------------

    def is_principal(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Principal", index)

    # -----------------------------------------------------------------

    def get_companions(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        companions = self.get_quantity("Companions", index)
        if companions is not None: companions = companions.split(",")
        return companions

    # -----------------------------------------------------------------

    def get_parent(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Parent", index)

    # -----------------------------------------------------------------

    def create_source(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Debugging
        log.debug("Creating an extended source for entry " + str(index) + " in the catalog ...")

        # Get the row
        row = self.get_row(index)

        # Get properties
        name = row["Name"]
        redshift = row["Redshift"]
        galaxy_type = row["Type"]
        position = self.get_position(index)
        names = row["Names"].split(",") if row["Names"] is not None else None
        distance = row["Distance"]
        inclination = row["Incl"]
        d25 = row["D25"]
        major = row["Major"]
        minor = row["Minor"]
        position_angle = row["Posangle"]
        principal = row["Principal"]
        companions = row["Companions"].split(",") if row["Companions"] is not None else None
        parent = row["Parent"]

        # Create a new ExtendedSource instance
        source = ExtendedSource(index=index, name=name, position=position, redshift=redshift, galaxy_type=galaxy_type,
                                names=names, distance=distance, inclination=inclination, d25=d25, major=major,
                                minor=minor, position_angle=position_angle, principal=principal, companions=companions,
                                parent=parent)

        # Return the source
        return source

    # -----------------------------------------------------------------

    def create_point_region(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        position = self.get_position(index)

        # Create a coordinate for the center and add it to the region
        meta = {"point": "x"}

        # Create and return the region
        region = SkyPointRegion(position.ra, position.dec, meta=meta)
        return region

    # -----------------------------------------------------------------

    def create_pixel_point_region(self, index, wcs):

        """
        This function ... 
        :param index:
        :param wcs:
        :return: 
        """

        # Sky point region
        region = self.create_point_region(index)

        # Convert and return
        return region.to_pixel(wcs)

    # -----------------------------------------------------------------

    def create_region(self, index, default_radius=parse_quantity("20 arcsec")):

        """
        This function ...
        :param index: 
        :param default_radius:
        :return: 
        """

        # Debugging
        log.debug("Creating region for extended source " + str(index+1) + " ...")

        # Get the row
        row = self.get_row(index)

        position = self.get_position(index)

        position_angle = row["Posangle"]

        if position_angle is not None:
            angle = Angle(position_angle.to("deg").value, "deg")
        else: angle = Angle(0.0, "deg")

        major = row["Major"]
        minor = row["Minor"]

        if major is None:

            color = "red"
            ra_radius = default_radius
            dec_radius = default_radius

        elif minor is None or position_angle is None:

            color = "green"
            ra_radius = 0.5 * major.to("arcsec")
            dec_radius = ra_radius

        else:

            color = "green"
            ra_radius = 0.5 * major.to("arcsec")
            dec_radius = 0.5 * minor.to("arcsec")

        # Get name
        name = row["Name"]
        principal = row["Principal"]

        text = name
        if principal: text += " (principal)"

        # Create ellipse shape
        center = position
        radius = SkyStretch(ra_radius, dec_radius)
        shape = SkyEllipseRegion(center, radius, angle)

        # Set meta information
        meta = {"text": text, "color": color, "index": index}
        shape.meta = meta

        # Return the region
        return shape

    # -----------------------------------------------------------------

    def create_pixel_region(self, index, wcs, default_radius=20.):

        """
        This function ...
        :param index: 
        :param wcs:
        :param default_radius:
        :return: 
        """

        # Calculate radius in sky coordinates (degrees)
        default_sky_radius = default_radius * wcs.average_pixelscale

        # Create sky region and return it
        sky_region = self.create_region(index, default_sky_radius)

        # Convert to pixel and return
        return sky_region.to_pixel(wcs)

    # -----------------------------------------------------------------

    def create_region_list(self, check_in_wcs=None, default_radius=parse_quantity("20 arcsec"), add_point=False):

        """
        This function ...
        :param check_in_wcs:
        :param default_radius:
        :param add_point:
        :return: 
        """

        # Inform the user
        log.info("Creating region list from extended source catalog ...")

        # Initialize the region list
        regions = SkyRegionList()

        # Loop over the entries in the catalog
        for index in range(len(self)):

            # Get position
            position = self.get_position(index)

            # If the source falls outside of the frame, skip it
            if check_in_wcs is not None and not check_in_wcs.contains(position): continue

            # Create the region
            region = self.create_region(index, default_radius=default_radius)

            # Add the region
            regions.append(region)

            # If add point
            if add_point:

                region = self.create_point_region(index)
                regions.append(region)

        # return the region list
        return regions

    # -----------------------------------------------------------------

    def create_regions(self, check_in_wcs=None, default_radius=parse_quantity("20 arcsec"), add_point=False):

        """
        THis function ...
        :param check_in_wcs: 
        :param default_radius:
        :param add_point:
        :return: 
        """

        return self.create_region_list(check_in_wcs=check_in_wcs, default_radius=default_radius, add_point=add_point)

    # -----------------------------------------------------------------

    def create_pixel_region_list(self, wcs, check_in_wcs=False, default_radius=20.0, add_point=False):

        """
        This function ...
        :param wcs: 
        :param check_in_wcs: 
        :param default_radius:
        :param add_point:
        :return: 
        """

        # Initialize the region list
        regions = PixelRegionList()

        # Loop over the entries in the catalog
        for index in range(len(self)):

            # Get position
            position = self.get_position(index)

            # If the source falls outside of the frame, skip it
            if check_in_wcs and not wcs.contains(position): continue

            # Create the pixel region
            region = self.create_pixel_region(index, wcs, default_radius=default_radius)

            # Add the region
            regions.append(region)

            # If add point
            if add_point:
                region = self.create_pixel_point_region(index, wcs)
                regions.append(region)

        # return the region list
        return regions

    # -----------------------------------------------------------------

    def create_pixel_regions(self, wcs, check_in_wcs=False, default_radius=20.0, add_point=False):

        """
        This function ...
        :param wcs: 
        :param check_in_wcs: 
        :param default_radius:
        :param add_point:
        :return: 
        """

        return self.create_pixel_region_list(wcs, check_in_wcs=check_in_wcs, default_radius=default_radius, add_point=add_point)

# -----------------------------------------------------------------
