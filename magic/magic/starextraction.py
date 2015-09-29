#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
import os.path
import inspect
from config import Config

# Import Astromagic modules
from ..core.star import Star

# Import astronomical modules
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier

# *****************************************************************

class StarExtractor(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "starextractor.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

        # Set the log level
        log.setLevel(self.config.log_level)

        # Initialize an empty list for the stars
        self.stars = []

    # *****************************************************************

    def run(self, image):

        """
        This function ...
        """

        # Get list of stars
        self.fetch_stars(image)

        # Find sources
        self.find_sources(image)

        # Fit ...
        self.fit_stars()

        # If requested, add a frame to the image with the stars
        #if self.config.make_star_frame: self.create_frame(image)

        # If requested, remove the stars
        if self.config.remove: self.remove_stars(image)

        # If requested, add the stars region
        if self.config.add_region: self.create_region(image)

    # *****************************************************************

    def fetch_stars(self, image):

        """
        This function ...
        """

        # Inform the user
        log.info("Fetching star positions from an online catalog")

        # Get the range of right ascension and declination of this image
        center, ra_span, dec_span = image.frames.selected(require_single=True).coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=self.config.fetching.catalogs)
        table = result[0]

        # Loop over all stars in the table
        for entry in table:

            # Get the PGC number
            ucac_id = entry["UCAC4"]

            # Get the right ascension and the declination
            position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # If this star does not lie within the frame, skip it
            if not image.frames.selected(require_single=True).contains(position): continue

            # Get the mean error on the right ascension and declination
            position_error = entry["ePos"]*u.mas
            ra_error = entry["e_RAJ2000"]*u.mas
            dec_error = entry["e_DEJ2000"]*u.mas

            # Get the magnitude in different bands
            k_mag = entry["Kmag"]*u.mag
            b_mag = entry["Bmag"]*u.mag
            v_mag = entry["Vmag"]*u.mag
            r_mag = entry["rmag"]*u.mag
            i_mag = entry["imag"]*u.mag

            # Create a star object
            star = Star(ucac_id=ucac_id, position=position, position_error=position_error, ra_error=ra_error,
                        dec_error=dec_error, k_mag=k_mag, b_mag=b_mag, v_mag=v_mag, r_mag=r_mag, i_mag=i_mag)

            # If requested, enable track record
            if self.config.track_record: star.enable_track_record()

            # Add the star to the list of stars
            self.stars.append(star)

        # Inform the user
        log.debug("Number of stars: " + str(len(self.stars)))

    # *****************************************************************

    def find_sources(self, image):

        """
        This function ...
        """

        # Inform the user
        log.info("Looking for sources near the star positions")

        # Get the currently selected frame
        frame_name = image.frames.get_selected(require_single=True)
        frame = image.frames[frame_name]

        # Loop over all stars in the list
        for star in self.stars:

            # Find a source
            star.find_source(frame, self.config.detection)

        # Inform the user
        log.debug("Success ratio: {0:.2f}%".format(self.have_source/len(self.stars)*100.0))

    # *****************************************************************

    def fit_stars(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Fitting analytical profiles to the sources")

        # Loop over all stars in the list
        for star in self.stars:

            # Find a source
            if star.has_source: star.fit_model(self.config.fitting)

        # Inform the user
        log.debug("Success ratio: {0:.2f}%".format(self.have_model/self.have_source*100.0))

    # *****************************************************************

    def remove_stars(self, image):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing the stars from the frame")

        # Get the currently selected frame
        frame_name = image.frames.get_selected(require_single=True)
        frame = image.frames[frame_name]

        # Loop over all stars in the list
        for star in self.stars:

            # Remove the star in the frame
            star.remove(frame)

    # *****************************************************************

    def create_frame(self, image):

        """
        This function ...
        """

        pass

    # *****************************************************************

    def create_region(self, image):

        """
        This function ...
        """

        pass

    # *****************************************************************

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_source
        return count

    # *****************************************************************

    @property
    def have_model(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_model
        return count

    # *****************************************************************

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fwhm of the fitted stars
        fwhm_list = []

        # Loop over all stars
        for star in self.stars:

            # If the star contains a model, add the fwhm of that model to the list
            if star.has_model: fwhm_list.append(star.fwhm)

        # Return the mean of the fwhm for all fitted stars
        return np.mean(fwhm_list)

# *****************************************************************
