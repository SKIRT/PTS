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
from ..core.masks import Mask
from ..core.star import Star
from ..tools import statistics
from ..core.vector import Position

# Import astronomical modules
from astropy import log
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
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

    def clear(self):

        """
        This function ...
        :return:
        """

        # Clear the stars list
        self.stars = []

    # *****************************************************************

    def run(self, frame, galaxyextractor=None):

        """
        This function ...
        """

        # Create a list of stars based on online catalogs
        self.fetch_stars(frame, galaxyextractor)

        # For each star, find a corresponding source in the image
        self.find_sources(frame)

        # Fit analytical models to the stars
        self.fit_stars()

        # If requested, remove the stars
        if self.config.remove: self.remove_stars(frame)

        # If requested, remove saturation in the image
        if self.config.remove_saturation: self.remove_saturation(frame)

    # *****************************************************************

    def fetch_stars(self, frame, galaxyextractor=None):

        """
        This function ...
        """

        # Inform the user
        log.info("Fetching star positions from an online catalog")

        # Get the range of right ascension and declination of this image
        center, ra_span, dec_span = frame.coordinate_range()

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=self.config.fetching.catalogs)
        table = result[0]

        # If a galaxyextractor is passed, get the positions of the galaxies in the frame
        if galaxyextractor is not None: galaxies = galaxyextractor.position_list(frame)
        else: galaxies = []

        # Loop over all stars in the table
        norms = []
        for entry in table:

            # Get the PGC number
            ucac_id = entry["UCAC4"]

            # Get the right ascension and the declination
            position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

            # If this star does not lie within the frame, skip it
            if not frame.contains(position): continue

            # Get the mean error on the right ascension and declination
            position_error = entry["ePos"]*u.mas
            ra_error = entry["e_RAJ2000"]*u.mas
            dec_error = entry["e_DEJ2000"]*u.mas

            # Loop over all galaxy positions in the list
            for galaxy in galaxies:

                # Calculate the distance between the star's position and the galaxy's center
                x_center, y_center = position.to_pixel(frame.wcs)
                difference = galaxy - Position(x=x_center, y=y_center)

                norms.append(difference.norm)

                if difference.norm <= self.config.fetching.min_difference_from_galaxy:

                    # Remove the position of this galaxy from the list
                    galaxies.remove(galaxy)
                    break

            else:

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
        if galaxyextractor is not None: log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(norm) for norm in sorted(norms)[:10]))

        # Inform the user
        log.debug("Number of stars: " + str(len(self.stars)))

    # *****************************************************************

    def find_sources(self, frame):

        """
        This function ...
        """

        # Inform the user
        log.info("Looking for sources near the star positions")

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

    def remove_stars(self, frame):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing the stars from the frame")

        # Get a list of the fwhm of all the stars that were fitted
        fwhm_list = self.fwhm_list

        # Determine the fwhm for the stars that were not fitted
        if self.config.removal.no_source_fwhm == "max": default_fwhm = max(fwhm_list)
        elif self.config.removal.no_source_fwhm == "mean": default_fwhm = np.mean(fwhm_list)
        elif self.config.removal.no_source_fwhm == "median": default_fwhm = np.median(fwhm_list)
        else: raise ValueError("Unkown measure for determining the fwhm of stars without detected source")

        # Inform the user
        log.debug("Default FWHM used when star could not be fitted: {0:.2f} pixels".format(default_fwhm))

        # Loop over all stars in the list
        for star in self.stars:

            # Remove the star in the frame
            star.remove(frame, self.config.removal, default_fwhm)

    # *****************************************************************

    def remove_saturation(self, frame):

        """
        This function ...
        :param image:
        :return:
        """

        # Inform the user
        log.info("Removing saturation from the frame")

        # TODO: allow "central_brightness" as config.saturation.criterion

        # Get a list of the fwhm of all the stars that were fitted
        fwhm_list = self.fwhm_list

        # Determine the fwhm for the stars that were not fitted
        if self.config.saturation.no_source_fwhm == "max": default_fwhm = max(fwhm_list)
        elif self.config.saturation.no_source_fwhm == "mean": default_fwhm = np.mean(fwhm_list)
        elif self.config.saturation.no_source_fwhm == "median": default_fwhm = np.median(fwhm_list)
        else: raise ValueError("Unkown measure for determining the fwhm of stars without detected source")

        # Get a list of the fluxes of the stars
        flux_list = self.flux_list

        # Calculate the minimal flux/central brightness
        minimum = statistics.cutoff(flux_list, self.config.saturation.brightest_method, self.config.saturation.limit)

        # Inform the user
        quantity = "flux" if self.config.saturation.criterion == "flux" else "central brightness"
        log.debug("Minimum value of the " + quantity + " for saturation removal: {0:.2f}".format(minimum))

        # Inform the user
        eligible = len([flux for flux in flux_list if flux >= minimum])
        log.debug("Number of stars eligible for saturation removal: " + str(eligible) + " ({0:.2f}%)".format(eligible/len(self.stars)*100.0))

        # Loop over all stars in the list
        removed = 0
        for star in self.stars:

            # If a source was not found for this star, skip it
            if not star.has_source: continue

            # Calculate the value (flux or brightness) for this star
            value = star.flux

            # Remove the saturation
            if value >= minimum:

                star.remove_saturation(frame, self.config.saturation, default_fwhm)
                if star.has_source: removed += 1

        # Inform the user
        log.debug("Removed saturation in " + str(removed) + " out of " + str(eligible) + " stars ({0:.2f}%)".format(removed/eligible*100.0))

    # *****************************************************************

    def create_frame(self, frame):

        """
        This function ...
        """

        return None

    # *****************************************************************

    def create_region(self, frame):

        """
        This function ...
        """

        return None

    # *****************************************************************

    def create_mask(self, frame):

        """
        This function ...
        :param frame:
        :return:
        """

        # Initialize a mask with the dimensions of the frame
        mask = Mask(np.zeros_like(frame))

        # Loop over all stars
        for star in self.stars:

            # If no source was found for the galaxy, skip it
            if not star.has_source: continue

            # Add this galaxy to the mask
            mask[star.source.cutout.y_min:star.source.cutout.y_max, star.source.cutout.x_min:star.source.cutout.x_max] = star.source.mask



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
    def fwhm_list(self):

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

        # Return the list
        return fwhm_list

    # *****************************************************************

    @property
    def flux_list(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fluxes of the stars
        flux_list = []

        # Loop over all stars
        for star in self.stars:

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            if star.has_source and star.source.is_estimated:

                # Subtract the background if necessary
                if not star.source.is_subtracted: star.source.subtract_background()

                # Add the flux to the list
                flux_list.append(star.flux)

        # Return the list
        return flux_list

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

    def table(self):

        """
        This function ...
        :return:
        """

        ids = []
        ascensions = []
        declinations = []
        sources = []
        models = []
        fwhms = []

        # Loop over all stars
        for star in self.stars:

            ids.append(star.ucac_id)
            ascensions.append(star.position.ra.value)
            declinations.append(star.position.dec.value)
            sources.append(star.has_source)
            models.append(star.has_model)
            if star.has_model: fwhms.append(star.fwhm)
            else: fwhms.append(None)

        # Create and return the table
        return Table([ids, ascensions, declinations, sources, models, fwhms], names=('UCAC-ID', 'RA', 'DEC', 'Source', 'Model', 'FWHM'), meta={'name': 'stars'})

# *****************************************************************
