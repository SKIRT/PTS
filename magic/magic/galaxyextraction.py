#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import numpy as np
from config import Config
import inspect

# Import Astromagic modules
from ..core import regions
from ..core.masks import Mask
from ..core.galaxy import Galaxy
from ..core.vector import Position

# Import astromagic modules
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier

# *****************************************************************

class GalaxyExtractor(object):

    """
    This class
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "galaxyextractor.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

        # Initialize an empty list for the galaxies
        self.galaxies = []

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the list of galaxies
        self.galaxies = []

    # *****************************************************************

    def run(self, frame):

        """
        This function ...
        """

        # Get list of galaxies
        self.fetch_galaxies(frame)

        # Find the sources
        self.find_sources(frame)

        # If requested, remove
        if self.config.remove: self.remove_galaxies(frame)

    # *****************************************************************

    def fetch_galaxies(self, frame):

        """
        This function ...
        :param image:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy positions from an online catalog")

        # Get the range of right ascension and declination of the image
        center, ra_span, dec_span = frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["galaxies", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=["VII/237"])
        table = result[0]

        # Loop over all galaxies in the table
        for entry in table:

            # Get the PGC number
            pgc_id = entry["PGC"]

            # Get the right ascension and the declination
            center = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # Get the names given to this galaxy
            names = entry["ANames"].split() if entry["ANames"] else None

            # Get the size of the galaxy
            ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
            diameter = np.power(10.0, entry["logD25"])*0.1*u.arcmin if entry["logD25"] else None

            # Get the size of major and minor axes
            major = diameter
            minor = diameter / ratio if diameter is not None and ratio is not None else None

            # Get the position angle of the galaxy
            pa = entry["PA"]*u.deg if entry["PA"] else None

            # Add a new galaxy to the list
            self.galaxies.append(Galaxy(pgc_id=pgc_id, center=center, names=names, major=major, minor=minor, pa=pa))

        # Define a function that returns the length of the major axis of the galaxy
        def major_axis(galaxy):

            if galaxy.major is None: return 0.0
            else: return galaxy.major

        # Indicate which galaxy is the principal galaxy
        galaxy = max(self.galaxies, key=major_axis)
        galaxy.principal = True

        # Inform the user
        log.debug("Number of galaxies: " + str(len(self.galaxies)))

    # *****************************************************************

    def find_sources(self, frame):

        """
        This function ...
        :param image:
        :return:
        """

        # Inform the user
        log.info("Looking for sources near the galaxy positions")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Do not perform the procedure for the principal galaxy
            if galaxy.principal: continue

            # Find a source
            galaxy.find_source(frame, self.config.detection)

        # Inform the user
        log.debug("Success ratio: {0:.2f}%".format(self.have_source/len(self.galaxies)*100.0))

    # *****************************************************************

    def remove_galaxies(self, frame):

        """
        This function ...
        :param image:
        :param galaxies:
        :param config:
        :return:
        """

        # Inform the user
        log.info("Removing the galaxies from the frame (except for the principal galaxy)")

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # Remove the galaxy from the frame
            if galaxy.source is not None: galaxy.remove(frame, self.config.removal)

    # *****************************************************************

    def create_region(self, frame):

        """
        This function ...
        :param image:
        :param stars:
        :param config:
        :return:
        """

        ra_list = []
        dec_list = []
        height_list = []
        width_list = []
        angle_list = []

        # Loop over all galaxies
        for galaxy in self.galaxies:

            ra_list.append(galaxy.center.ra.value)
            dec_list.append(galaxy.center.dec.value)

            if galaxy.major is not None:

                width = galaxy.major.to("arcsec").value

                if galaxy.minor is not None: height = galaxy.minor.to("arcsec").value
                else: height = width

            else:

                width = self.config.region.default_radius * frame.pixelscale.value
                height = self.config.region.default_radius * frame.pixelscale.value

            angle = galaxy.pa.value if galaxy.pa is not None else 0.0

            height_list.append(height)
            width_list.append(width)
            angle_list.append(angle)

        # Create a region and return it
        return regions.ellipses(ra_list, dec_list, height_list, width_list, angle_list)

    # *****************************************************************

    def create_mask(self, frame, principal=False):

        """
        This function ...
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(frame))

        # Loop over all galaxies
        for galaxy in self.galaxies:

            # If no source was found for the galaxy, skip it
            if not galaxy.has_source: continue

            # Add this galaxy to the mask
            mask[galaxy.source.cutout.y_min:galaxy.source.cutout.y_max, galaxy.source.cutout.x_min:galaxy.source.cutout.x_max] = galaxy.source.mask

        # Return the mask
        return mask

    # *****************************************************************

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for galaxy in self.galaxies: count += galaxy.has_source
        return count

    # *****************************************************************

    def position_list(self, frame):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the galaxy positions
        position_list = []

        # Loop over the galaxies
        for galaxy in self.galaxies:

            # Calculate the pixel coordinate in the frame and add it to the list
            x_center, y_center = galaxy.center.to_pixel(frame.wcs)
            position_list.append(Position(x=x_center, y=y_center))

        # Return the list
        return position_list

    # *****************************************************************

    def table(self):

        """
        This function ...
        :return:
        """

        pass

# *****************************************************************
