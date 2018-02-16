#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.point Contains the PointSourceCatalog class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.basics.log import log
from ..core.pointsource import PointSource
from ..basics.coordinate import SkyCoordinate
from ..region.point import SkyPointRegion
from ...core.units.parsing import parse_quantity
from ..region.list import PixelRegionList, SkyRegionList
from ..region.circle import SkyCircleRegion

# -----------------------------------------------------------------

class PointSourceCatalog(SmartTable):

    """
    This class ...
    """

    # Add columns
    _column_info = OrderedDict()
    _column_info["Catalog"] = (str, None, "Original catalog")
    _column_info["ID"] = ( str, None, "ID in the catalog")
    _column_info["RA"] = (float, "deg", "Right ascension")
    _column_info["DEC"] = (float, "deg", "Declination")
    _column_info["RA error"] = (float, "mas", "Error on right ascension")
    _column_info["DEC error"] = (float, "mas", "Error on declination")
    _column_info["Confidence"] = (int, None, "Confidence level")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSourceCatalog, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_coordinate(self, coordinate, catalog=None, id=None, ra_error=None, dec_error=None):

        """
        This function ...
        :param coordinate:
        :param catalog:
        :param id:
        :param ra_error:
        :param dec_error:
        :return:
        """

        self.add_entry(catalog, id, coordinate.ra, coordinate.dec, ra_error, dec_error, confidence_level=None)

    # -----------------------------------------------------------------

    def add_entry(self, catalog, id, ra, dec, ra_error, dec_error, confidence_level):

        """
        This function ...
        :param catalog:
        :param id:
        :param ra:
        :param dec:
        :param ra_error:
        :param dec_error:
        :param confidence_level
        :return:
        """

        # Add
        values = [catalog, id, ra, dec, ra_error, dec_error, confidence_level]
        self.add_row(values)

    # -----------------------------------------------------------------

    def coordinates(self):

        """
        This function ...
        :return:
        """

        coordinates = []

        unit = self.column_unit("RA")
        assert unit == self.column_unit("DEC")

        for index in range(len(self)):

            ra = self["RA"][index]
            dec = self["DEC"][index]

            coordinate = SkyCoordinate(ra=ra, dec=dec, unit=unit)
            coordinates.append(coordinate)

        # Return the coordinates
        return coordinates

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index], dec=self["DEC"][index], unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def get_catalog(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Catalog", index)

    # -----------------------------------------------------------------

    def get_id(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("ID", index)

    # -----------------------------------------------------------------

    def get_ra_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("RA error", index)

    # -----------------------------------------------------------------

    def get_dec_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("DEC error", index)

    # -----------------------------------------------------------------

    def get_confidence(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Confidence", index)

    # -----------------------------------------------------------------

    def create_source(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Debugging
        log.debug("Creating a point source for entry " + str(index) + " in the catalog ...")

        #print(index, len(self))

        # Get the row
        row = self.get_row(index)

        #print(row)

        # Get properties
        catalog = row["Catalog"]
        id = row["ID"]
        position = self.get_position(index)
        ra_error = row["RA error"]
        dec_error = row["DEC error"]
        confidence = row["Confidence"]

        # Check for which bands magnitudes are defined
        magnitudes = dict()
        magnitude_errors = {}
        for name in row:
            if "magnitude" in name:
                band = name.split(" magnitude")[0]
                magnitudes[band] = row[name]
                magnitude_errors[band] = row[name + " error"]

        # Create a new PointSource instance
        source = PointSource(index=index, position=position, catalog=catalog, id=id, ra_error=ra_error,
                             dec_error=dec_error, confidence=confidence, magnitudes=magnitudes,
                             magnitude_errors=magnitude_errors)

        # Return the source
        return source

    # -----------------------------------------------------------------

    def create_point_region(self, index, color="white"):

        """
        This function ...
        :param index: 
        :param color:
        :return: 
        """

        # Create meta information for the position
        meta = {"point": "x", "color": color}

        position = self.get_position(index)

        # Create the position and add it to the region
        region = SkyPointRegion(position.ra, position.dec, meta=meta)
        return region

    # -----------------------------------------------------------------

    def create_pixel_point_region(self, index, wcs, color="white"):

        """
        This function ...
        :param index: 
        :param wcs: 
        :param color:
        :return: 
        """

        # SKy region
        region = self.create_point_region(index, color=color)

        # Convert and return
        return region.to_pixel(wcs)

    # -----------------------------------------------------------------

    def create_region(self, index, radius=parse_quantity("10 arcsec"), color="white"):

        """
        This function ...
        :param index: 
        :param radius:
        :param color:
        :return: 
        """

        # Debugging
        log.debug("Creating region for point source " + str(index + 1) + " ...")

        # Get the row and position
        #row = self.get_row(index)
        position = self.get_position(index)

        # Determine the color, based on the detection level
        #if source.has_model: color = "blue"
        #elif source.has_detection: color = "green"
        #else: color = "red"

        # Convert the source index to a string
        text = str(index)

        # Create meta information
        meta = {"color": color, "text": text, "index": index}

        # Create hape
        shape = SkyCircleRegion(position, radius, meta=meta)

        #return the shape
        return shape

    # -----------------------------------------------------------------

    def create_pixel_region(self, index, wcs, radius=10., color="white"):

        """
        THis function ...
        :param index: 
        :param wcs: 
        :param radius:
        :param color:
        :return: 
        """

        # Calculate radius in sky coordinates (degrees)
        sky_radius = radius * wcs.average_pixelscale

        # Create sky region and return it
        sky_region = self.create_region(index, sky_radius, color=color)

        # Convert to pixel and return
        return sky_region.to_pixel(wcs)

    # -----------------------------------------------------------------

    def create_region_list(self, check_in_wcs=None, radius=parse_quantity("10 arcsec"), add_point=False, color="white"):

        """
        This function ...
        :return: 
        """

        # Initialize region list
        regions = SkyRegionList()

        # Loop over the sources in the catalog
        for index in range(len(self)):

            # Get position
            position = self.get_position(index)

            # If the stars falls outside of the frame, skip
            if check_in_wcs is not None and not check_in_wcs.contains(position): continue

            # Create region
            region = self.create_region(index, radius=radius, color=color)

            # Add the region
            regions.append(region)

            if add_point:

                region = self.create_point_region(index)
                regions.append(region)

        # Return the region list
        return regions

    # -----------------------------------------------------------------

    def create_regions(self, check_in_wcs=None, radius=parse_quantity("10. arcsec"), add_point=False, color="white"):

        """
        THis function ...
        :return: 
        """

        return self.create_region_list(check_in_wcs=check_in_wcs, radius=radius, add_point=add_point, color=color)

    # -----------------------------------------------------------------

    def create_pixel_region_list(self, wcs, check_in_wcs=False, radius=10., add_point=False, color="white"):

        """
        This function ...
        :param wcs:
        :param check_in_wcs:
        :param radius:
        :param add_point:
        :param color:
        :return: 
        """

        # Initializ region list
        regions = PixelRegionList()

        # Calculate radius in sky coordinates (degrees)
        sky_radius = radius * wcs.average_pixelscale

        # Loop over the sources in the catalog
        for index in range(len(self)):

            # Get position
            position = self.get_position(index)

            # If the star fall outside of the frame, skip
            if check_in_wcs and not wcs.contains(position): continue

            # Create pixel region
            region = self.create_pixel_region(index, wcs, radius=sky_radius, color=color)

            # Add the region
            regions.append(region)

            # Add point
            if add_point:

                region = self.create_point_region(index, color=color)
                regions.append(region)

        # Return the region list
        return regions

    # -----------------------------------------------------------------

    def create_pixel_regions(self, wcs, check_in_wcs=False, radius=10., add_point=False, color="white"):

        """
        This function ...
        :param wcs:
        :param check_in_wcs:
        :param radius:
        :param add_point:
        :param color:
        :return: 
        """

        return self.create_pixel_region_list(wcs, check_in_wcs=check_in_wcs, radius=radius, add_point=add_point, color=color)

# -----------------------------------------------------------------
