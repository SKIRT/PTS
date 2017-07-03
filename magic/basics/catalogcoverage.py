#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalogcoverage Contains the CatalogCoverage class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.vector import Position, Extent
from ..basics.coordinate import SkyCoordinate
from ..region.rectangle import PixelRectangleRegion, SkyRectangleRegion
from ...core.tools import tables, introspection
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class CatalogCoverage(object):

    """
    This class ...
    """

    def __init__(self, name):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the user catalogs directory
        catalogs_user_path = os.path.join(introspection.pts_user_dir, "magic", "catalogs")

        # Determine the path to the galaxy directory
        galaxy_path = os.path.join(catalogs_user_path, name)

        # Load the 'ranges' table
        self.ranges_table_path = os.path.join(galaxy_path, "ranges.dat")

        if fs.is_file(self.ranges_table_path):

            self.table = tables.from_file(self.ranges_table_path)

            # Get the ranges of the first entry
            galaxy_ra = self.table["Central right ascension"][0]
            galaxy_dec = self.table["Central declination"][0]
            galaxy_ra_span = self.table["Right ascension span"][0]
            galaxy_dec_span = self.table["Declination span"][0]

            self.first_box = PixelRectangleRegion(Position(galaxy_ra, galaxy_dec), Extent(galaxy_ra_span, galaxy_dec_span))

        else:

            data = [[],[],[],[]]
            names = ["Central right ascension", "Central declination", "Right ascension span", "Declination span"]
            self.table = tables.new(data, names)

    # -----------------------------------------------------------------

    def matches(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        # If the center position does not fall within the ranges of the first entry
        # OPTIMIZATION !!!! ASSUMES GALAXY IS ALWAYS MORE OR LESS IN THE CENTER OF THE FRAME FOR ALL IMAGES
        if self.first_box.contains(box.center): return True

    # -----------------------------------------------------------------

    def covers_position(self, position):

        """
        This function ...
        :param position: is Position(x=ra(in deg), y=dec(in deg))
        :return:
        """

        # Loop over all boxes
        for box in self.boxes:

            # If at least one of the boxes contains this position, return True
            if box.contains(position): return True

    # -----------------------------------------------------------------

    def covers(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        # List of boxes from ranges in the table
        self.boxes = []

        ra = box.center.ra.to("deg").value
        dec = box.center.dec.to("deg").value
        ra_span = 2. * box.radius.x.to("deg").value
        dec_span = 2. * box.radius.y.to("deg").value

        # Loop over all entries in the table
        for i in range(len(self.table)):

            entry_ra = self.table["Central right ascension"][i]
            entry_dec = self.table["Central declination"][i]
            entry_ra_span = self.table["Right ascension span"][i]
            entry_dec_span = self.table["Declination span"][i]

            # Check whether the box properties are the same (neglecting small roundoff errors)
            same_ra = np.isclose(entry_ra, ra)
            same_dec = np.isclose(entry_dec, dec)
            same_ra_span = np.isclose(entry_ra_span, ra_span)
            same_dec_span = np.isclose(entry_dec_span, dec_span)

            ## First check for identical entry
            identical_box_already_processed = same_ra and same_dec and same_ra_span and same_dec_span

            # If a frame with identical ranges has already been processed before
            if identical_box_already_processed: return True

            entry_center = SkyCoordinate(entry_ra, entry_dec)
            entry_radius = Extent(0.5 * entry_ra_span, 0.5 * entry_dec_span)

            entry_box = SkyRectangleRegion(entry_center, entry_radius)

            self.boxes.append(entry_box)

        # Loop over corners of the box
        for corner in box.corners:

            # Loop over all boxes
            for box in self.boxes:

                # If at least one of the boxes contains this corner, this corner is OK
                if box.contains(corner): break

            # If a break is not encountered: this corner is not covered by any of the boxes
            else: return False

        # If a break IS encountered for each corner, each corner is covered by one or more boxes, so this box
        # does not extent the range of these other boxes
        return True

    # -----------------------------------------------------------------

    def add_box(self, box):

        """
        This function ...
        :param box:
        :return:
        """

        ra = box.center.x
        dec = box.center.y
        ra_span = 2. * box.radius.x
        dec_span = 2. * box.radius.y

        self.table.add_row([ra, dec, ra_span, dec_span])

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Save table over the old one
        tables.write(self.table, self.ranges_table_path)

# -----------------------------------------------------------------
