#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.starextraction Contains the StarExtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import config

# Import astronomical modules
import aplpy
import astropy.io.fits as pyfits
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, MaskedColumn
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.convolution import Gaussian2DKernel

# Import the relevant AstroMagic classes and modules
from .basics import Position, Extent, Mask, Region
from .core import Source
from .sky import Star
from .tools import statistics, fitting, regions

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable
from ..core.tools import tables

# -----------------------------------------------------------------

class StarExtractor(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(StarExtractor, self).__init__(config, "magic")

        # -- Attributes --

        # Initialize an empty list for the stars
        self.stars = []

        # Initialize an empty list to contain the manual sources
        self.manual_sources = []

        # Set the frame to None
        self.frame = None
        self.original_frame = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.input_mask = None

        # The input catalog
        self.input_catalog = None

        # Set the mask to None
        self.mask = None

        # Reference to the galaxy extractor
        self.galaxy_extractor = None

    # -----------------------------------------------------------------

    def run(self, frame, input_mask, galaxyextractor=None, catalog=None):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup(frame, input_mask, galaxyextractor, catalog)

        # 2. Find and remove the stars
        self.find_fit_and_remove_stars()

        # 3. If requested, find and remove saturated stars
        if self.config.find_saturation: self.find_and_remove_saturation()

        # 4. If requested, find and remove apertures
        if self.config.find_apertures: self.find_and_remove_apertures()

        # 5. If specified, remove manually selected stars
        if self.config.manual_region is not None: self.set_and_remove_manual()

        # 6. Writing phase
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame, input_mask, galaxyextractor=None, catalog=None):

        """
        This function ...
        """

        # Call the setup function of the base class
        super(StarExtractor, self).setup()

        # Make a local reference to the frame and 'bad' mask
        self.frame = frame
        self.original_frame = frame.copy()
        self.input_mask = input_mask
        self.input_catalog = catalog

        # Make a local reference to the galaxy extractor (if any)
        self.galaxy_extractor = galaxyextractor

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask(np.zeros_like(self.frame))

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the star extractor")

        # Clear the list of stars
        self.stars = []

        # Clear the list of manual sources
        self.manual_sources = []

        # Clear the frame and the mask
        self.frame = None
        self.mask = None

    # -----------------------------------------------------------------

    def find_fit_and_remove_stars(self):

        """
        This function ...
        :return:
        """

        # If no input catalog was given
        if self.input_catalog is None:

            # Create a list of stars based on online catalogs
            if self.config.fetching.use_catalog_file: self.load_catalog()
            else: self.fetch_stars()

        # If an input catalog was given
        else: self.process_catalog()

        # Set special stars
        if self.config.special_region is not None: self.set_special()

        # Set ignore stars
        if self.config.ignore_region is not None: self.set_ignore()

        # For each star, find a corresponding source in the image
        self.find_sources()

        # Fit analytical models to the stars
        self.fit_stars()

        # If requested, remove the stars
        if self.config.remove: self.remove_stars()

    # -----------------------------------------------------------------

    def find_and_remove_saturation(self):

        """
        This function ...
        :return:
        """

        # If requested, remove saturation in the image
        if self.config.remove_saturation: self.remove_saturation()

    # -----------------------------------------------------------------

    def find_and_remove_apertures(self):

        """
        This function ...
        :return:
        """

        # Find the apertures
        self.find_apertures()

        # If requested, remove apertures
        if self.config.remove_apertures: self.remove_apertures()

    # -----------------------------------------------------------------

    def set_and_remove_manual(self):

        """
        This function ...
        :return:
        """

        # Set manual stars
        self.set_manual()

        # If requested, remove the manually specified stars
        if self.config.remove_manual: self.remove_manual()

    # -----------------------------------------------------------------

    def load_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_input_path(self.config.fetching.catalog_path)

        # Inform the user
        self.log.info("Loading stellar catalog from file " + path)

        # Load the catalog
        self.input_catalog = tables.from_file(path)

        # Process the catalog
        self.process_catalog()

    # -----------------------------------------------------------------

    def process_catalog(self):

        """
        This function ...
        :return:
        """

        # Create the list of stars
        for i in range(len(self.input_catalog)):

            # Get the star properties
            catalog = self.input_catalog["Catalog"][i]
            id = self.input_catalog["Id"][i]
            ra = self.input_catalog["Right ascension"][i]
            dec = self.input_catalog["Declination"][i]
            ra_error = self.input_catalog["Right ascension error"][i] * u.deg
            dec_error = self.input_catalog["Declination error"][i] * u.deg
            on_galaxy = self.input_catalog["On galaxy"][i]
            confidence_level = self.input_catalog["Confidence level"][i]

            magnitudes = {}
            magnitude_errors = {}

            for name in self.input_catalog.colnames:
                if "magnitude" in name:

                    band = name.split(" magnitude")[0]
                    magnitudes[band] = self.input_catalog[name][i] * u.mag
                    magnitude_errors[band] = self.input_catalog[name + " error"][i] * u.mag

            # Create a sky coordinate for the star position
            position = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='fk5')

            # Create a star object
            star = Star(catalog=catalog, id=id, position=position, ra_error=ra_error,
                        dec_error=dec_error, magnitudes=magnitudes, magnitude_errors=magnitude_errors)

            # Set other attributes
            star.on_galaxy = on_galaxy
            star.confidence_level = confidence_level

            # What to do with:
            # Ignored
            # Detected
            # Fitted
            # ?

            # Add the star to the list
            self.stars.append(star)

    # -----------------------------------------------------------------

    def fetch_stars(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Fetching star positions from an online catalog")

        if isinstance(self.config.fetching.catalogs, basestring):
            catalogs = [self.config.fetching.catalogs]
        elif isinstance(self.config.fetching.catalogs, config.Sequence):
            catalogs = self.config.fetching.catalogs
        else: raise ValueError("Invalid option for 'catalogs', should be a string or a list of strings")

        # Get the range of right ascension and declination of this image
        try:
            center, ra_span, dec_span = self.frame.coordinate_range()
        except AssertionError:
            self.log.warning("The coordinate system and pixelscale do not match")
            center, ra_span, dec_span = self.frame.coordinate_range(silent=True)

        # Create a new Vizier object and set the row limit to -1 (unlimited)
        viz = Vizier(keywords=["stars", "optical"])
        viz.ROW_LIMIT = -1

        # Loop over the different catalogs
        for catalog in catalogs:

            # Copy the list of stars already added from other catalogs
            previous_stars = list(self.stars)

            galaxies = list(self.galaxy_extractor.galaxies) if self.galaxy_extractor is not None else []

            # Initialize a list for the stars obtained from this catalog
            stars = []

            # Inform the user
            self.log.debug("Querying " + catalog + " catalog")

            # Query Vizier and obtain the resulting table
            result = viz.query_region(center, width=ra_span, height=dec_span, catalog=catalog)
            table = result[0]

            number_of_stars = 0
            number_of_stars_in_frame = 0
            number_of_new_stars = 0

            # Loop over all stars in the table
            distances = []
            #for entry in table: # slower
            for i in range(len(table)):

                # Initialize an empty dictionary to contain the magnitudes and the magnitude errors
                magnitudes = {}
                mag_errors = {}

                # Get the ID of this star in the catalog
                if catalog == "UCAC4": star_id = table["UCAC4"][i]
                elif catalog == "NOMAD": star_id = table["NOMAD1"][i]
                elif catalog == "II/246": star_id = table["_2MASS"][i]
                else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

                # Get the position of the star as a SkyCoord object
                position = coord.SkyCoord(ra=table["_RAJ2000"][i], dec=table["_DEJ2000"][i], unit=(u.deg, u.deg), frame='fk5')

                number_of_stars += 1

                # If this star does not lie within the frame, skip it
                if not self.frame.contains(position, self.config.transformation_method): continue

                number_of_stars_in_frame += 1

                # Get the mean error on the right ascension and declination
                if catalog == "UCAC4" or catalog == "NOMAD":

                    ra_error = table["e_RAJ2000"][i] * u.mas
                    dec_error = table["e_DEJ2000"][i] * u.mas

                elif catalog == "II/246":

                    error_maj = table["errMaj"][i] * u.arcsec
                    error_min = table["errMin"][i] * u.arcsec
                    error_theta = Angle(table["errPA"][i], u.deg)

                    # Temporary: use only the major axis error (convert the error ellipse into a circle)
                    ra_error = error_maj
                    dec_error = error_maj

                else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

                # Loop over all galaxies
                for galaxy in galaxies:

                    # Calculate the pixel position of the galaxy
                    galaxy_position = galaxy.pixel_position(self.frame.wcs)

                    # Calculate the distance between the star's position and the galaxy's center
                    x_center, y_center = position.to_pixel(self.frame.wcs, mode=self.config.transformation_method)
                    difference = galaxy_position - Position(x=x_center, y=y_center)

                    distance = difference.norm

                    # Add the star-galaxy distance to the list of distances
                    distances.append(distance)

                    # The principal galaxy/galaxies
                    if galaxy.principal:

                        # Check whether the star-galaxy distance is smaller than a certain threshold
                        if distance <= self.config.fetching.min_distance_from_galaxy.principal:
                            break

                    # Companion galaxies
                    elif galaxy.companion:

                        if distance <= self.config.fetching.min_distance_from_galaxy.companion:
                            # Remove the position of this galaxy from the list (one star is already identified with it)
                            galaxies.remove(galaxy)
                            break

                    # All other galaxies in the frame
                    else:

                        if distance <= self.config.fetching.min_distance_from_galaxy.other:
                            # Remove the position of this galaxy from the list (one star is already identified with it)
                            galaxies.remove(galaxy)
                            break

                # If a break is not encountered
                else:

                    # Get the magnitude in different bands
                    for name in table.colnames:

                        # If this column name does not end with "mag", skip it
                        if not name.endswith("mag"): continue

                        # If the column name contains more than one character before "mag", skip it
                        if len(name.split("mag")[0]) > 1: continue

                        # Get the name of the band
                        band = name.split("mag")[0]

                        if not isinstance(table[name][i], np.ma.core.MaskedConstant):

                            # Add the magnitude in this band to the dictionary
                            magnitudes[band] = table[name][i] * u.mag

                            # Check whether an error on the magnitude is present
                            if "e_" + name in table.colnames and not isinstance(table["e_" + name][i], np.ma.core.MaskedConstant):

                                # If so, add it to the mag_errors dictionary
                                mag_errors[band] = table["e_" + name][i] * u.mag

                    # Create a star object
                    star = Star(catalog=catalog, id=star_id, position=position, ra_error=ra_error,
                                dec_error=dec_error, magnitudes=magnitudes, magnitude_errors=mag_errors)

                    # Check whether this star is on top of the galaxy, and label it so (by default, star.on_galaxy is False)
                    if self.galaxy_extractor is not None:

                        star.on_galaxy = self.galaxy_extractor.principal.contains(star.pixel_position(self.frame.wcs, self.config.transformation_method))

                    # If requested, enable track record
                    if self.config.track_record: star.enable_track_record()

                    # If there are already stars in the list, check for correspondences with the current stars
                    for saved_star in previous_stars:

                        # Calculate the distance between the star already in the list and the new star
                        difference = saved_star.pixel_position(self.frame.wcs, self.config.transformation_method) - star.pixel_position(self.frame.wcs, self.config.transformation_method)

                        if difference.norm < 3.0:

                            # Inform the user
                            self.log.debug("Star " + star_id + " could be identified with star " + saved_star.id + " from the " + saved_star.catalog + " catalog")

                            saved_star.confidence_level += 1
                            # Remove stars from previous_stars
                            previous_stars.remove(saved_star)
                            break

                    # If no other stars are in the list yet or no corresponding star was found (no break was
                    # encountered), just add all stars of the current catalog
                    else:

                        number_of_new_stars += 1

                        # Inform the user
                        self.log.debug("Adding star " + star_id + " at " + str(position.to_string("hmsdms")))
                        stars.append(star)

            self.log.debug("Number of stars that were in the catalog: " + str(number_of_stars))
            self.log.debug("Number of stars that fell within the frame: " + str(number_of_stars_in_frame))
            self.log.debug("Number of stars that were only present in this catalog: " + str(number_of_new_stars))

            # Add the new stars to the stars list
            self.stars += stars

        # Inform the user
        if self.galaxy_extractor is not None: self.log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(distance) for distance in sorted(distances)[:10]))

        # Inform the user
        self.log.debug("Number of stars: " + str(len(self.stars)))

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        self.log.info("Looking for sources near the star positions")

        # Loop over all sky objects in the list
        for skyobject in self.stars:

            # If this sky object should be ignored, skip it
            if skyobject.ignore: continue

            # Find a source
            try:
                skyobject.find_source(self.frame, self.config.detection)

            except Exception as e:

                #import traceback

                self.log.error("Error when finding source")
                #print(type(e))
                #print(e)
                #traceback.print_exc()

                if self.config.plot_track_record_if_exception:

                    if skyobject.has_track_record: skyobject.track_record.plot()
                    else: self.log.warning("Track record is not enabled")

                self.log.error("Continuing with next source")

        # Inform the user
        self.log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.stars), self.have_source/len(self.stars)*100.0))

    # -----------------------------------------------------------------

    def fit_stars(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Fitting analytical profiles to the sources")

        # Loop over all stars in the list
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Check if the star has a source (has been detected)
            if not star.has_source and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                center, radius, angle = star.ellipse_parameters(self.frame.wcs, self.frame.pixelscale, self.config.fitting.initial_radius)

                # Create a source object
                source = Source.from_ellipse(self.frame, center, radius, angle, self.config.fitting.background_outer_factor)

            else: source = None

            # Find a source
            if star.has_source or source is not None: star.fit_model(self.config.fitting, source)

        # If requested, perform sigma-clipping to the list of FWHM's to filter out outliers
        if self.config.fitting.sigma_clip_fwhms:

            mean, median, stddev = statistics.sigma_clipped_statistics(self.fwhms, self.config.fitting.fwhm_sigma_level)
            lower = median - self.config.fitting.fwhm_sigma_level * stddev
            upper = median + self.config.fitting.fwhm_sigma_level * stddev

            # Loop over all stars for which a model was found
            for star in self.stars:

                # Ignore stars without model
                if not star.has_model: continue

                # Remove the model if its FWHM is clipped out
                if star.fwhm > upper or star.fwhm < lower: star.model = None

        # Inform the user
        self.log.debug("Found a model for {0} out of {1} stars with source ({2:.2f}%)".format(self.have_model, self.have_source, self.have_model/self.have_source*100.0))

    # -----------------------------------------------------------------

    def remove_stars(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Removing the stars from the frame")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm

        # Inform the user
        self.log.debug("Default FWHM used when star could not be fitted: {0:.2f} pixels".format(default_fwhm))

        if self.config.removal.method == "model":

            # Calculate the relative differences between the ampltides of the fitted models and the corresponding sources
            differences = np.array(self.amplitude_differences) * 100.0

            print(np.mean(differences))
            print(np.median(differences))
            print(np.std(differences))
            print(differences)

        # Loop over all stars in the list
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
            if not self.config.removal.remove_foreground and self.galaxy_extractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

            # Remove the star in the frame
            star.remove(self.frame, self.mask, self.config.removal, default_fwhm)

    # -----------------------------------------------------------------

    def remove_saturation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Removing saturation from the frame")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm

        # Set the number of stars where saturation was removed to zero initially
        removed = 0

        # Detect and remove saturation for all stars
        if self.config.saturation.method == "all":

            # Inform the user on the number of stars that have a source
            self.log.debug("Number of stars with source = " + str(self.have_source))

            # Loop over all stars
            for star in self.stars:

                # If this star should be ignored, skip it
                if star.ignore: continue

                # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
                if not self.config.saturation.remove_foreground and self.galaxy_extractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

                # If a model was not found for this star, skip it unless the remove_if_not_fitted flag is enabled
                if not star.has_model and not self.config.saturation.remove_if_not_fitted: continue
                if star.has_model: assert star.has_source

                # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Find a saturation source and remove it from the frame
                self.config.saturation.centroid_table_path = self.full_output_path(self.config.saturation.centroid_table_path)
                success = star.remove_saturation(self.frame, self.mask, self.config.saturation, default_fwhm)
                if success: star.has_saturation = True
                removed += success

            # Inform the user
            self.log.debug("Removed saturation in " + str(removed) + " out of " + str(self.have_source) + " stars with source ({0:.2f}%)".format(removed/self.have_source*100.0))

        # Detect and remove saturation for the brightest stars
        elif self.config.saturation.method == "brightest":

            # TODO: allow "central_brightness" as config.saturation.criterion

            # Get a list of the fluxes of the stars
            flux_list = self.fluxes

            # Calculate the minimal flux/central brightness
            minimum = statistics.cutoff(flux_list, self.config.saturation.brightest_method, self.config.saturation.limit)

            # Inform the user
            quantity = "flux" if self.config.saturation.criterion == "flux" else "central brightness"
            self.log.debug("Minimum value of the " + quantity + " for saturation removal: {0:.2f}".format(minimum))

            # Inform the user
            eligible = len([flux for flux in flux_list if flux >= minimum])
            self.log.debug("Number of stars eligible for saturation removal: " + str(eligible) + " ({0:.2f}%)".format(eligible/len(self.stars)*100.0))

            # Inform the user on the number of stars that have a source
            self.log.debug("Number of stars with source = " + str(self.have_source))

            # Loop over all stars
            for star in self.stars:

                # If this star should be ignored, skip it
                if star.ignore: continue

                # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
                if not self.config.saturation.remove_foreground and self.galaxy_extractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

                # If a model was not found for this star, skip it unless the remove_if_not_fitted flag is enabled
                if not star.has_model and not self.config.saturation.remove_if_not_fitted: continue
                if star.has_model: assert star.has_source

                # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
                if not star.has_source and not self.config.saturation.remove_if_undetected: continue

                # Calculate the value (flux or brightness) for this star
                try: value = star.flux
                except AttributeError: value = 0.0

                # Remove the saturation if the value is greater than the minimum value or the star has no source and 'remove_if_undetected' is enabled
                if value >= minimum or (self.config.saturation.remove_if_undetected and not star.has_source):

                    # Find a saturation source and remove it from the frame
                    self.config.saturation.centroid_table_path = self.full_output_path(self.config.saturation.centroid_table_path)
                    success = star.remove_saturation(self.frame, self.mask, self.config.saturation, default_fwhm)
                    if success: star.has_saturation = True
                    removed += success

            # Inform the user
            self.log.debug("Removed saturation in {0} out of {1} stars ({2:.2f}%)".format(removed, eligible, removed/eligible*100.0))

        # Unkown saturation
        else: raise ValueError("Unknown method (should be 'brightest' or 'all'")

    # -----------------------------------------------------------------

    def find_apertures(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        self.log.info("Constructing elliptical apertures regions to encompass saturated stars")

        # Loop over all stars
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If the star does not have saturation, continue
            if star.has_saturation: star.find_aperture(self.frame, self.config.apertures)

    # -----------------------------------------------------------------

    def remove_apertures(self):

        """
        This function ...
        :param frame:
        :param factor:
        :return:
        """

        # Inform the user
        self.log.info("Replacing aperture regions by the estimated background")

        # Loop over all stars
        for star in self.stars:

            # If the object does not have an aperture, skip it
            if not star.has_aperture: continue

            # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
            if star.on_galaxy and self.config.aperture_removal.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = self.config.aperture_removal.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            if star.on_galaxy and self.config.aperture_removal.polynomial_on_galaxy: interpolation_method = "polynomial"
            else: interpolation_method = self.config.aperture_removal.interpolation_method

            # Expansion factor
            expansion_factor = self.config.aperture_removal.expansion_factor

            # Create a source object
            # Get the parameters of the elliptical aperture
            x_center, y_center = star.aperture.positions[0]
            center = Position(x=x_center, y=y_center)

            major = star.aperture.a * expansion_factor
            minor = star.aperture.b * expansion_factor

            radius = Extent(x=major, y=minor)

            # theta is in radians
            angle = Angle(star.aperture.theta, u.rad)

            # Create a source
            source = Source.from_ellipse(self.frame, center, radius, angle, self.config.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background(interpolation_method, sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

            # Update the mask
            self.mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

    # -----------------------------------------------------------------

    def set_special(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine the full path to the special region file
        path = self.full_input_path(self.config.special_region)

        # Inform the user
        self.log.info("Setting special region from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all objects
        for skyobject in self.stars:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(self.frame.wcs)

            # Set special if position is covered by the mask
            if special_mask.masks(position): skyobject.special = True

    # -----------------------------------------------------------------

    def set_ignore(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Determine the full path to the ignore region file
        path = self.full_input_path(self.config.ignore_region)

        # Inform the user
        self.log.info("Setting region to ignore for subtraction from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        ignore_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Loop over all objects
        for skyobject in self.stars:

            # Get the position of this object in pixel coordinates
            position = skyobject.pixel_position(self.frame.wcs)

            # Ignore if position is covered by the mask
            if ignore_mask.masks(position): skyobject.ignore = True

    # -----------------------------------------------------------------

    def set_manual(self):

        """
        This function ...
        """

        # Determine the full path to the manual region file
        path = self.full_input_path(self.config.manual_region)

        # Inform the user
        self.log.info("Setting region for manual star extraction from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)

        # Loop over the shapes in the region
        for shape in region:

            # Get the center and radius of the shape (can be a circle or an ellipse)
            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            # Create a source
            source = Source.from_ellipse(self.frame, Position(x_center, y_center), Extent(x_radius, y_radius), Angle(0.0, u.deg), self.config.manual.background_outer_factor)

            # Add the source to the list of manual sources
            self.manual_sources.append(source)

    # -----------------------------------------------------------------

    def remove_manual(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Removing manually specified stars from the frame")

        # Loop over each item in the list of manual sources
        for source in self.manual_sources:

            # Estimate the background for the source
            source.estimate_background(self.config.manual.interpolation_method, self.config.manual.sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

    # -----------------------------------------------------------------

    @property
    def region(self):

        """
        This function ...
        :return:
        """

        # TODO: improve this function

        type = "sky"

        # Initialize lists
        position_list = []
        radius_list = []
        color_list = []

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.stars:

            position_list.append(star.position)

            if star.has_model:

                fwhm = star.fwhm
                color = "green"

            else:

                fwhm = default_fwhm
                color = "red"

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * self.config.region.sigma_level

            # Add the radius (in arcseconds) and the color the appropriate list
            radius_list.append(radius * self.frame.pixelscale)
            color_list.append(color)

        # Create a region
        region = Region.circles(position_list, radius_list, color_list)

        if type == "sky": return region
        elif type == "image": return region.as_imagecoord(self.frame.wcs.to_header())
        else: raise ValueError("Type should be either 'sky' or 'image'")

    # -----------------------------------------------------------------

    def write_region(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Determine the full path to the region file
        path = self.full_output_path(self.config.writing.region_path)
        annotation = self.config.writing.region_annotation

        # Inform the user
        self.log.info("Writing stars region to " + path)

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.stars:

            # Get the center in pixel coordinates
            #x_center, y_center = star.position.to_pixel(self.frame.wcs, origin=0)
            center = star.pixel_position(self.frame.wcs)

            if star.has_source:

                if star.has_model:

                    fwhm = star.fwhm
                    color = "blue"

                else:

                    fwhm = default_fwhm
                    color = "green"

            else:

                fwhm = default_fwhm
                color = "red"

            if annotation == "flux":

                if star.has_source and star.source.has_background:

                    text = "text = {" + str(int(round(star.flux))) + "}"

                else: text = ""

            elif annotation == "has_source":

                text = "text = {" + str(star.has_source) + "}"

            elif annotation == "has_background":

                if star.has_source: text = "text = {" + str(star.source.has_background) + "}"
                else: text = ""

            elif annotation is None: text = ""
            else: raise ValueError("Invalid option for annotation")

            # If the FWHM is defined, draw a circle for the star and draw a cross for its peak position (if defined)
            if fwhm is not None:

                # Calculate the radius in pixels
                radius = fwhm * statistics.fwhm_to_sigma * self.config.region.sigma_level

                # Draw a cross for the peak position
                if star.has_source and star.source.has_peak:

                    point_suffix = " # point = x " + text
                    print("image;point({},{})".format(star.source.peak.x, star.source.peak.y) + point_suffix, file=f)

                # Show a circle for the star
                color_suffix = " # color = " + color
                print("image;circle({},{},{})".format(center.x, center.y, radius) + color_suffix, file=f)

            # If the FWHM is undefined, simply draw a point for the star's position (e.g. when this function is called
            # after the fetch_stars method)
            else: print("image;point({},{})".format(center.x, center.y), file=f)

            # Aperture created from saturation mask
            if star.has_aperture:

                ap_x_center, ap_y_center = star.aperture.positions[0]
                major = star.aperture.a
                minor = star.aperture.b
                angle = star.aperture.theta * math.pi / 180

                aperture_suffix = " # color = white"

                print("image;ellipse({},{},{},{},{})".format(ap_x_center, ap_y_center, major, minor, angle) + aperture_suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_aperture_region(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the aperture region file
        path = self.full_output_path(self.config.writing.aperture_region_path)

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Loop over the stars
        for star in self.stars:

            # Skip stars that do not have an aperture
            if not star.has_aperture: continue

            ap_x_center, ap_y_center = star.aperture.positions[0]
            major = star.aperture.a
            minor = star.aperture.b
            angle = star.aperture.theta / math.pi * 180

            aperture_suffix = " # color = blue"

            # Calculate the difference between the aperture center and the star position (in number of pixels)
            difference = star.pixel_position(self.frame.wcs) - Position(ap_x_center, ap_y_center)
            aperture_suffix += " text = {" + str(difference.norm) + "}"

            print("image;ellipse({},{},{},{},{})".format(ap_x_center, ap_y_center, major, minor, angle) + aperture_suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the table file
        path = self.full_output_path(self.config.writing.table_path)

        # Inform the user
        self.log.info("Writing table to " + path)

        # Write the table to file
        self.table.write(path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.catalog_path)

        # Inform the user
        self.log.info("Writing stellar catalog to " + path)

        # Write the catalog to file
        self.catalog.write(path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        """

        # Determine the full path to the masked frame file
        path = self.full_output_path(self.config.writing.masked_frame_path)

        # Inform the user
        self.log.info("Writing masked frame to " + path)

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(path)

    # -----------------------------------------------------------------

    def write_cutouts(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the cutouts directory
        directory_path = self.full_output_path(self.config.writing.cutouts_path)

        # Inform the user
        self.log.info("Writing cutout boxes to " + directory_path)

        # Keep track of the number of stars encountered
        without_source = 0
        with_source = 0
        with_model = 0
        with_saturation = 0

        # Calculate the default FWHM based on the stars that could be fitted
        default_fwhm = self.fwhm

        # Loop over all stars
        for star in self.stars:

            # Check if saturation has been detected for this star
            if star.has_saturation:

                # Save the cutout as a FITS file
                path = os.path.join(directory_path, "star_saturation_" + str(with_saturation) + ".fits")
                star.source.save(path, self.original_frame)

                # Increment the counter of the number of stars with saturation
                with_saturation += 1

            # Check if a model has been found for this star
            elif star.has_model:

                # Save the cutout as a FITS file
                path = os.path.join(directory_path, "star_model_" + str(with_model) + ".fits")
                star.source.save(path)

                # Increment the counter of the number of stars that could be fitted
                with_model += 1

            # Check if a source was found for this star
            elif star.has_source:

                # Save the cutout as a FITS file
                path = os.path.join(directory_path, "star_source_" + str(with_source) + ".fits")
                star.source.save(path)

                # Increment the counter of the number of stars that could be detected
                with_source += 1

            # If no source was found for this star
            else:

                # Create a source for the desired sigma level and outer factor
                sigma_level = self.config.writing.cutouts_undected_sigma_level
                outer_factor = self.config.writing.cutouts_undetected_outer_factor
                source = star.source_at_sigma_level(self.frame, default_fwhm, sigma_level, outer_factor)

                # Save the cutout as a FITS file
                path = os.path.join(directory_path, "star_nosource_" + str(with_source) + ".fits")
                source.save(path)

                # Increment the counter of the number of stars without source
                without_source += 1

    # -----------------------------------------------------------------

    def write_result(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the resulting FITS file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        self.log.info("Writing resulting frame to " + path)

        # Write out the resulting frame
        self.frame.save(path)

    # -----------------------------------------------------------------

    @property
    def positions(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the galaxies
        for skyobject in self.stars:

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(skyobject.pixel_position(self.frame.wcs))

        # Return the list
        return positions

    # -----------------------------------------------------------------

    @property
    def have_source(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_source
        return count

    # -----------------------------------------------------------------

    @property
    def have_model(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_model
        return count

    # -----------------------------------------------------------------

    @property
    def have_saturation(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_saturation
        return count

    # -----------------------------------------------------------------

    @property
    def have_aperture(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_aperture
        return count

    # -----------------------------------------------------------------

    @property
    def fwhms(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fwhm of the fitted stars
        fwhms = []

        # Loop over all stars
        for star in self.stars:

            # If the star contains a model, add the fwhm of that model to the list
            if star.has_model: fwhms.append(star.fwhm)

        # Return the list
        return fwhms

    # -----------------------------------------------------------------

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fluxes of the stars
        fluxes = []

        # Loop over all stars
        for star in self.stars:

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            if star.has_source and star.source.has_background:

                # Add the flux to the list
                fluxes.append(star.flux)

        # Return the list
        return fluxes

    # -----------------------------------------------------------------

    @property
    def amplitude_differences(self):

        """
        This function ...
        :return:
        """

        # Initialize
        differences = []

        # Loop over all stars
        for star in self.stars:

            # If the star was not fitted, skip it
            if not star.has_model: continue

            # Determine the amplitude and the position of the center of the model
            amplitude_model = star.model.amplitude
            center = star.source.cutout.rel_position(fitting.center(star.model))

            # Convert into integers
            x = int(round(center.x))
            y = int(round(center.y))

            # Calculate the value of the source at the model's center position
            amplitude_source = star.source.subtracted[y, x]

            # Calculate the difference of the amplitudes
            difference = abs(amplitude_model - amplitude_source)
            rel_difference = difference / amplitude_source

            # Add the relative difference to the list
            differences.append(rel_difference)

        # Return the list of differences
        return differences

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # If the list of FWHM values is empty (the stars were not fitted yet), return None
        if len(self.fwhms) == 0: return None

        # Determine the default FWHM and return it
        if self.config.fwhm.measure == "max": return max(self.fwhms)
        elif self.config.fwhm.measure == "mean": return np.mean(self.fwhms)
        elif self.config.fwhm.measure == "median": return np.median(self.fwhms)
        else: raise ValueError("Unkown measure for determining the default FWHM")

    # -----------------------------------------------------------------

    @property
    def kernel(self):

        """
        This function ...
        :return:
        """

        # Create a Gaussian convolution kernel and return it
        sigma = self.fwhm * statistics.fwhm_to_sigma
        return Gaussian2DKernel(sigma)

    # -----------------------------------------------------------------

    @property
    def table(self):

        """
        This function ...
        :return:
        """

        # Initialize empty lists for the table columns
        catalogs = []
        ids = []
        ascensions = []
        declinations = []
        sources = []
        models = []
        fwhms = []

        # Loop over all stars
        for star in self.stars:

            catalogs.append(star.catalog)
            ids.append(star.id)
            ascensions.append(star.position.ra.value)
            declinations.append(star.position.dec.value)
            sources.append(star.has_source)
            models.append(star.has_model)
            if star.has_model: fwhms.append(star.fwhm)
            else: fwhms.append(None)

        # Create and return the table
        data = [catalogs, ids, ascensions, declinations, sources, models, fwhms]
        names = ('CATALOG', 'ID', 'RA', 'DEC', 'Source', 'Model', 'FWHM')
        meta = {'name': 'stars'}
        return tables.new(data, names, meta)

    # -----------------------------------------------------------------

    @property
    def catalog(self):

        """
        This function ...
        :return:
        """

        # Initialize empty lists for the table columns
        catalogs = []
        ids = []
        right_ascensions = []
        declinations = []
        ra_errors = []
        dec_errors = []
        magnitudes = {}
        magnitude_errors = {}
        on_galaxy_list = []
        confidence_levels = []

        ignore_list = []
        have_source_list = []
        have_model_list = []

        # Loop over the stars to see for which bands magnitudes are defined
        for star in self.stars:

            for band in star.magnitudes:

                # Skip this band if we have already encountered it
                if band in magnitudes: continue

                # Create an empty list for this band
                magnitudes[band] = []
                magnitude_errors[band] = []

        # Loop over all stars
        for star in self.stars:

            catalogs.append(star.catalog)
            ids.append(star.id)
            right_ascensions.append(star.position.ra.value)
            declinations.append(star.position.dec.value)
            ra_errors.append(star.ra_error.value)
            dec_errors.append(star.dec_error.value)

            for band in magnitudes:

                if band in star.magnitudes:

                    magnitudes[band].append(star.magnitudes[band].value)
                    if band in star.magnitude_errors: magnitude_errors[band].append(star.magnitude_errors[band].value)
                    else: magnitude_errors[band].append(None)

                else:

                    magnitudes[band].append(None)
                    magnitude_errors[band].append(None)

            on_galaxy_list.append(star.on_galaxy)
            confidence_levels.append(star.confidence_level)

            ignore_list.append(star.ignore)
            have_source_list.append(star.has_source)
            have_model_list.append(star.has_model)

        # Create and return the table
        data = [catalogs, ids, right_ascensions, declinations, ra_errors, dec_errors, on_galaxy_list, confidence_levels]
        names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'On galaxy', 'Confidence level']

        magnitude_column_names = []
        for band in magnitudes:

            # Values
            #column = MaskedColumn(magnitudes[band], mask=[mag is None for mag in magnitudes[band]])
            #data.append(column)
            data.append(magnitudes[band])
            column_name = band + " magnitude"
            names.append(column_name)

            magnitude_column_names.append(column_name)

            # Errors
            #column = MaskedColumn(magnitude_errors[band], mask=[mag is None for mag in magnitude_errors[band]])
            #data.append(column)
            data.append(magnitude_errors[band])
            column_name = band + " magnitude error"
            names.append(column_name)

            magnitude_column_names.append(column_name)

        data.append(ignore_list)
        names.append("Ignored")

        data.append(have_source_list)
        names.append("Detected")

        data.append(have_model_list)
        names.append("Fitted")

        meta = {'name': 'stars'}
        table = Table(data, names=names, meta=meta, masked=True)

        # Set units
        table["Right ascension"].unit = "deg"
        table["Declination"].unit = "deg"
        table["Right ascension error"].unit = "deg"
        table["Declination error"].unit = "deg"
        for name in magnitude_column_names:
            table[name].unit = "mag"

        # Return the catalog
        return table

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Create a HDU from this frame with the image header
        hdu = pyfits.PrimaryHDU(self.frame, self.frame.wcs.to_header())

        # Create a figure canvas
        figure = plt.figure(figsize=(20, 20))

        # Create a figure from this frame
        plot = aplpy.FITSFigure(hdu, figure=figure)

        # Plot in color scale
        plot.show_colorscale()

        # Add a color bar if requested
        if self.config.plotting.show_colorbar: plot.add_colorbar()

        # Add these shapes to the plot
        plot.show_regions(self.region)

        # Show the plot
        plt.show()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # If requested, write out a table with the star properties
        if self.config.write_table: self.write_table()

        # If requested, write out the star region
        if self.config.write_region: self.write_region()

        # If requested, write out the aperture region
        if self.config.write_aperture_region: self.write_aperture_region()

        # If requested, write out the frame where the stars are masked
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the star cutout boxes
        if self.config.write_cutouts: self.write_cutouts()

        # If requested, write out the result
        if self.config.write_result: self.write_result()

        # If requested, write out the stellar catalog
        if self.config.write_catalog: self.write_catalog()

# -----------------------------------------------------------------
