#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import os.path
import numpy as np
from config import Config

from frames import Frame
import catalogs
import regions
import analysis

# Import astromagic modules
from astropy import log

# *****************************************************************

class StarRemover(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        if config is None:

            # Load the default configurations for the star remover
            config_path = os.path.join(os.getcwd(), "config", "starremover.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

        # Set attributes to None
        self.number = None
        self.region = None

    # *****************************************************************

    def run(self, image):

        fetcher = StarFetcher(self.config)

        region = fetcher.region

        # Get the name of the active frame
        frame_name = image.frames.get_selected(require_single=True)

        sources = analysis.find_sources_in_region(image.frames[frame_name].data, region, self.config.model, self.config.detection_method, plot=self.config.plot, plot_custom=self.config.plot_custom)
        log.info("Number of sources = " + str(len(sources)))

        # Remove duplicates
        unique_sources = analysis.remove_duplicate_sources(sources)
        log.info("Number of unique sources = " + str(len(unique_sources)))

        self.number = len(unique_sources)

        if self.number > 0: self.region = regions.ellipses_from_coordinates(unique_sources)

# *****************************************************************

class StarFetcher(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        if config is None:

            # Load the default configurations for the star fetcher
            config_path = os.path.join(os.getcwd(), "config", "starfetcher.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

        # Set attributes to None
        self.number = None
        self.region = None

    # *****************************************************************

    def run(self, image):

        # Get the range of RA and DEC of this image
        ra_center, dec_center, size_ra_deg, size_dec_deg = image._get_coordinate_range()

        # Search for stars
        radius_in_arcsec = self.config.radius * image.pixelscale
        box = (ra_center, dec_center, size_ra_deg, size_dec_deg)
        if self.config.column_filters is None: region = catalogs.fetch_objects_in_box(box, self.config.catalog, ["stars"], radius_in_arcsec)
        else: region = catalogs.fetch_objects_in_box(box, self.config.catalog, ["optical", "stars"], radius_in_arcsec, column_filters={"Vmag":"<20"})

        if self.config.galaxy_name is not None:

            # Search for the position of the specified galaxy in the image
            galaxy_region = catalogs.fetch_object_by_name(self.config.galaxy_name, self.config.radius)
            original_len = len(region)
            region = regions.subtract(region, galaxy_region, 3.0, image.header)

            if len(region) == original_len - 1: log.info("Removed star position that matches the galaxy center")
            elif len(region) < original_len - 1: log.warning("Removed multiple star positions")

        # Set the number of unique stars
        self.number = len(region)

        # Set the region
        self.region = region

    # *****************************************************************