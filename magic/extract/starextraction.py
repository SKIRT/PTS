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
from astropy.coordinates import Angle
from astropy.convolution import Gaussian2DKernel

# Import the relevant AstroMagic classes and modules
from ..basics import Position, Extent, Mask, Region, Ellipse
from ..core import Source
from ..sky import Star
from ..tools import statistics, fitting, regions, catalogs

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import tables

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

        # The stellar catalog
        self.catalog = None

        # The statistics table
        self.statistics = None

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

        # 4. If specified, remove manually selected stars
        if self.config.manual_region is not None: self.set_and_remove_manual()

        # 5. Set the statistics
        if not self.config.fetching.use_statistics_file: self.set_statistics()

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
        self.catalog = catalog

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
        self.log.info("Clearing the star extractor ...")

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
        if self.catalog is None:

            # Create a list of stars based on online catalogs
            if self.config.fetching.use_catalog_file: self.import_catalog()
            else: self.fetch_catalog()

        # If an input catalog was given
        self.load_stars()

        # Import statistics if statistics file is specified
        if self.config.fetching.use_statistics_file:

            # Import the statistics
            self.import_statistics()

            # If requested, remove the stars
            if self.config.remove: self.remove_from_statistics()

        else:

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

        # Find saturated stars in the frame
        self.find_saturation()

        # If requested, remove saturation in the frame
        if self.config.remove_saturation: self.remove_saturation()

        # If requested, remove apertures that encompass the saturation sources
        if self.config.saturation.remove_apertures: self.remove_contours()

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

    def import_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_input_path(self.config.fetching.catalog_path)

        # Inform the user
        self.log.info("Importing stellar catalog from file: " + path)

        # Load the catalog
        self.catalog = tables.from_file(path)

    # -----------------------------------------------------------------

    def load_stars(self):

        """
        This function creates the star list from the star catalog.
        :return:
        """

        # Inform the user
        self.log.info("Loading the stars from the catalog ...")

        # Get masks
        special_mask = self.special_mask
        ignore_mask = self.ignore_mask

        # Copy the list of galaxies, so that we can removed already encounted galaxies (TODO: change this to use
        # an 'encountered' list as well
        encountered_galaxies = [False] * len(self.galaxy_extractor.galaxies)

        galaxy_pixel_position_list = []
        galaxy_type_list = []
        for galaxy in self.galaxy_extractor.galaxies:

            galaxy_pixel_position_list.append(galaxy.pixel_position(self.frame.wcs, self.config.transformation_method))
            if galaxy.principal: galaxy_type_list.append("principal")
            elif galaxy.companion: galaxy_type_list.append("companion")
            else: galaxy_type_list.append("other")

        # Keep track of the distances between the stars and the galaxies
        distances = []

        on_galaxy_column = []

        # Create the list of stars
        for i in range(len(self.catalog)):

            # Get the star properties
            catalog = self.catalog["Catalog"][i]
            star_id = self.catalog["Id"][i]
            ra = self.catalog["Right ascension"][i]
            dec = self.catalog["Declination"][i]
            ra_error = self.catalog["Right ascension error"][i] * u.mas
            dec_error = self.catalog["Declination error"][i] * u.mas
            confidence_level = self.catalog["Confidence level"][i]

            # Check for which bands magnitudes are defined
            magnitudes = {}
            magnitude_errors = {}
            for name in self.catalog.colnames:
                if "magnitude" in name:
                    band = name.split(" magnitude")[0]
                    magnitudes[band] = self.catalog[name][i] * u.mag
                    magnitude_errors[band] = self.catalog[name + " error"][i] * u.mag

            # Create a sky coordinate for the star position
            position = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='fk5')

            # Create a star object
            star = Star(i, catalog=catalog, id=star_id, position=position, ra_error=ra_error,
                        dec_error=dec_error, magnitudes=magnitudes, magnitude_errors=magnitude_errors)

            # Get the position of the star in pixel coordinates
            pixel_position = star.pixel_position(self.frame.wcs, self.config.transformation_method)

            # -- Checking for foreground or surroudings of galaxy --

            if "On galaxy" in self.catalog.colnames: star_on_galaxy = self.catalog["On galaxy"][i]
            else:

                # Check whether this star is on top of the galaxy, and label it so (by default, star.on_galaxy is False)
                if self.galaxy_extractor is not None: star_on_galaxy = self.galaxy_extractor.principal.contains(pixel_position)
                else: star_on_galaxy = False
                on_galaxy_column.append(star_on_galaxy)

            # -- Cross-referencing with the galaxies in the frame --

            # Loop over all galaxies
            if self.config.fetching.cross_reference_with_galaxies:

                # If a match is found with one of the galaxies, skip this star
                if matches_galaxy_position(pixel_position, galaxy_pixel_position_list, galaxy_type_list, encountered_galaxies, self.config.fetching.min_distance_from_galaxy, distances): continue

            # Set other attributes
            star.on_galaxy = star_on_galaxy
            star.confidence_level = confidence_level

            # Enable track record if requested
            if self.config.track_record: star.enable_track_record()

            # Set attributes based on masks (special and ignore)
            if special_mask is not None: star.special = special_mask.masks(pixel_position)
            if ignore_mask is not None: star.ignore = ignore_mask.masks(pixel_position)

            # If the input mask masks this star's position, skip it (don't add it to the list of stars)
            if self.input_mask is not None and self.input_mask.masks(pixel_position): continue

            # Don't add stars which are indicated as 'not stars'
            if self.config.manual_indices.not_stars is not None and i in self.config.manual_indices.not_stars: continue

            # Add the star to the list
            self.stars.append(star)

        # Add the 'on_galaxy' column to the catalog if necessary
        if "On galaxy" not in self.catalog.colnames: self.catalog["On galaxy"] = on_galaxy_column

        # Inform the user
        if self.config.fetching.cross_reference_with_galaxies: self.log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(distance) for distance in sorted(distances)[:10]))

    # -----------------------------------------------------------------

    def import_statistics(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the statistics file
        path = self.full_input_path(self.config.fetching.statistics_path)

        # Inform the user
        self.log.info("Importing stellar statistics from file: " + path)

        # Load the catalog
        self.statistics = tables.from_file(path)

        encountered = [False] * len(self.stars)

        # Loop over all entries in the statistics table
        for i in range(len(self.statistics)):

            index = self.statistics["Star index"][i]

            # Loop over all stars
            for j in range(len(self.stars)):

                # Skip this star if it has been encountered
                if encountered[j]: continue

                if self.stars[j].index == index:

                    encountered[j] = True

                    # TODO: work in progress

                    #self.stars[j].source =

                    break

            # If break is not encountered
            else: raise RuntimeError("The statistics file does not correspond to the catalog (star index does not exist)")

            #"Star index", "Detected", "Fitted", "Saturated", "FWHM", "Ignore"

    # -----------------------------------------------------------------

    def remove_from_statistics(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def fetch_catalog(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Fetching star positions from online catalogs ...")

        # Check whether the 'catalogs' setting defines a single catalog name or a list of such names
        if isinstance(self.config.fetching.catalogs, basestring): catalog_list = [self.config.fetching.catalogs]
        elif isinstance(self.config.fetching.catalogs, config.Sequence): catalog_list = self.config.fetching.catalogs
        else: raise ValueError("Invalid option for 'catalogs', should be a string or a list of strings")

        # Create the star catalog
        self.catalog = catalogs.create_star_catalog(self.frame, catalog_list)

        # Inform the user
        self.log.debug("Number of stars: " + str(len(self.catalog)))

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Inform the user
        self.log.info("Looking for sources near the star positions ...")

        # Loop over all stars in the list
        for star in self.stars:

            # If this sky object should be ignored, skip it
            if star.ignore: continue

            # Find a source
            try: star.find_source(self.frame, self.config.detection)
            except Exception as e:

                #import traceback

                self.log.error("Error when finding source")
                print(type(e))
                print(e)
                #traceback.print_exc()

                if self.config.plot_track_record_if_exception:

                    if star.has_track_record: star.track_record.plot()
                    else: self.log.warning("Track record is not enabled")

                self.log.error("Continuing with next source")

        # Inform the user
        self.log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.stars), self.have_source / len(self.stars) * 100.0))

    # -----------------------------------------------------------------

    def fit_stars(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Fitting analytical profiles to the sources ...")

        # Loop over all stars in the list
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Check if the star has a source (has been detected)
            if not star.has_source and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                center, radius, angle = star.ellipse_parameters(self.frame.wcs, self.frame.pixelscale, self.config.fitting.initial_radius)

                # Create a source object
                ellipse = Ellipse(center, radius, angle)
                source = Source.from_ellipse(self.frame, ellipse, self.config.fitting.background_outer_factor)

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
        self.log.info("Removing the stars from the frame ...")

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
            force = self.config.manual_indices.remove_stars is not None and star.index in self.config.manual_indices.remove_stars
            star.remove(self.frame, self.mask, self.config.removal, default_fwhm, force=force)

    # -----------------------------------------------------------------

    def find_saturation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Looking for saturated stars ...")

        # Inform the user on the number of stars that have a source
        self.log.debug("Number of stars with source = " + str(self.have_source))

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm

        # Set the number of stars where saturation was removed to zero initially
        success = 0

        # Loop over all stars
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If the index of this star is in the 'not_saturation' list, skip it
            if self.config.manual_indices.not_saturation is not None and star.index in self.config.manual_indices.not_saturation: continue

            # If remove_foreground is disabled and the star's position falls within the galaxy mask, we skip it
            if not self.config.saturation.remove_foreground and self.galaxy_extractor.mask.masks(star.pixel_position(self.frame.wcs)): continue

            # If a model was not found for this star, skip it unless the remove_if_not_fitted flag is enabled
            if not star.has_model and not self.config.saturation.remove_if_not_fitted: continue
            if star.has_model: assert star.has_source

            # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
            if not star.has_source and not self.config.saturation.remove_if_undetected: continue

            # Find a saturation source and remove it from the frame
            star.find_saturation(self.frame, self.original_frame, self.config.saturation, default_fwhm)
            success += star.has_saturation

        # Inform the user
        self.log.debug("Found saturation in " + str(success) + " out of " + str(self.have_source) + " stars with source ({0:.2f}%)".format(success / self.have_source * 100.0))

    # -----------------------------------------------------------------

    def remove_saturation(self):

        """
        This function ...
        :return:
        """

        # Loop over all stars
        for star in self.stars:

            # Skip stars for which saturation was not detected
            if not star.has_saturation: continue

            # Remove the saturation of this star in the frame
            star.remove_saturation(self.frame, self.mask, self.config.saturation)

    # -----------------------------------------------------------------

    def remove_contours(self):

        """
        This function ...
        :param frame:
        :param factor:
        :return:
        """

        # Inform the user
        self.log.info("Replacing regions within elliptical contours with the estimated background ...")

        # Loop over all stars
        for star in self.stars:

            # If the object does not have an contour, skip it
            if not star.has_contour: continue

            # Determine whether we want the background to be sigma-clipped when interpolating over the (saturation) source
            if star.on_galaxy and self.config.saturation.aperture_removal.no_sigma_clip_on_galaxy: sigma_clip = False
            else: sigma_clip = self.config.saturation.aperture_removal.sigma_clip

            # Determine whether we want the background to be estimated by a polynomial if we are on the galaxy
            if star.on_galaxy and self.config.saturation.aperture_removal.polynomial_on_galaxy: interpolation_method = "polynomial"
            else: interpolation_method = self.config.saturation.aperture_removal.interpolation_method

            # Expansion factor
            expansion_factor = self.config.saturation.aperture_removal.expansion_factor

            # Create a source
            ellipse = Ellipse(star.contour.center, star.contour.radius * expansion_factor, star.contour.angle)
            source = Source.from_ellipse(self.frame, ellipse, self.config.saturation.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background(interpolation_method, sigma_clip)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

            # Update the mask
            self.mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

    # -----------------------------------------------------------------

    @property
    def special_mask(self):

        """
        This function ...
        :param path:
        :return:
        """

        # If no special region is defined
        if self.config.special_region is None: return None

        # Determine the full path to the special region file
        path = self.full_input_path(self.config.special_region)

        # Inform the user
        self.log.info("Setting special region from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        special_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Return the mask
        return special_mask

    # -----------------------------------------------------------------

    @property
    def ignore_mask(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # If no ignore region is defined
        if self.config.ignore_region is None: return None

        # Determine the full path to the ignore region file
        path = self.full_input_path(self.config.ignore_region)

        # Inform the user
        self.log.info("Setting region to ignore for subtraction from " + path)

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)
        ignore_mask = Mask(region.get_mask(shape=self.frame.shape))

        # Return the mask
        return ignore_mask

    # -----------------------------------------------------------------

    def set_manual(self):

        """
        This function ...
        """

        # Determine the full path to the manual region file
        path = self.full_input_path(self.config.manual_region)

        # Inform the user
        self.log.info("Setting region for manual star extraction from " + path + " ...")

        # Load the region and create a mask from it
        region = Region.from_file(path, self.frame.wcs)

        # Loop over the shapes in the region
        for shape in region:

            # Get the center and radius of the shape (can be a circle or an ellipse)
            x_center, y_center, x_radius, y_radius, angle = regions.ellipse_parameters(shape)

            # Create a source
            ellipse = Ellipse(Position(x_center, y_center), Extent(x_radius, y_radius), Angle(angle, u.deg))
            source = Source.from_ellipse(self.frame, ellipse, self.config.manual.background_outer_factor)

            # Add the source to the list of manual sources
            self.manual_sources.append(source)

    # -----------------------------------------------------------------

    def remove_manual(self):

        """
        This function ...
        """

        # Inform the user
        self.log.info("Removing manually specified stars from the frame ...")

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

    def write_star_region(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Determine the full path to the star region file
        path = self.full_output_path(self.config.writing.star_region_path)

        # Inform the user
        self.log.info("Writing star region to " + path + " ...")

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm

        # Loop over all galaxies
        for star in self.stars:

            # Get the center in pixel coordinates
            center = star.pixel_position(self.frame.wcs, self.config.transformation_method)

            # Determine the color, based on the detection level
            if star.has_model: color = "blue"
            elif star.has_source: color = "green"
            else: color = "red"

            # Determine the FWHM
            fwhm = default_fwhm if not star.has_model else star.fwhm

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * self.config.removal.sigma_level

            # Convert the star index to a string
            text = str(star.index)

            # Show a circle for the star
            suffix = " # "
            color_suffix = "color = " + color
            text_suffix = "text = {" + text + "}"
            suffix += color_suffix + " " + text_suffix
            print("image;circle({},{},{})".format(center.x, center.y, radius) + suffix, file=f)

            # Draw a cross for the peak position
            if star.has_source and star.source.has_peak:

                suffix = " # "
                point_suffix = "point = x"
                suffix += point_suffix
                print("image;point({},{})".format(star.source.peak.x, star.source.peak.y) + suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_saturation_region(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the saturation region file
        path = self.full_output_path(self.config.writing.saturation_region_path)

        # Inform the user
        self.log.info("Writing saturation region to " + path + " ...")

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Loop over all stars
        for star in self.stars:

            # Skip stars without saturation
            if not star.has_saturation: continue

            # Convert the star index to a string
            text = str(star.index)

            # Get aperture properties
            center = star.contour.center
            major = star.contour.major
            minor = star.contour.minor
            angle = star.contour.angle.degree

            # Write to region file
            suffix = " # "
            color_suffix = "color = white"
            text_suffix = "text = {" + text + "}"
            suffix += color_suffix + " " + text_suffix
            print("image;ellipse({},{},{},{},{})".format(center.x, center.y, major, minor, angle) + suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_catalog(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.catalog_path)

        # Inform the user
        self.log.info("Writing stellar catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.catalog, path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the statistics file
        path = self.full_output_path(self.config.writing.statistics_path)

        # Inform the user
        self.log.info("Writing stellar statistics to " + path + " ...")

        # Write the catalog to file
        tables.write(self.statistics, path)

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        """

        # Determine the full path to the masked frame file
        path = self.full_output_path(self.config.writing.masked_frame_path)

        # Inform the user
        self.log.info("Writing masked frame to " + path + " ...")

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

        sigma_level = 3.0
        outer_factor = 1.5
        method = "polynomial"

        shape = Extent(21, 21)

        # Determine the full path to the cutouts directory
        directory_path = self.full_output_path(self.config.writing.cutouts_path)

        # Inform the user
        self.log.info("Writing cutout boxes to " + directory_path + " ...")

        # Calculate the default FWHM based on the stars that could be fitted
        default_fwhm = self.fwhm

        # Loop over all stars
        for star in self.stars:

            # -- Saturation sources ---

            # Check if saturation has been detected for this star
            if star.has_saturation:

                # Determine the path
                path = os.path.join(directory_path, "saturation_" + str(star.index) + ".fits")

                # Save the saturation source as a FITS file
                star.saturation.save(path)

            # -- PSF sources ---

            # Check if a model has been found for this star
            if star.has_model:

                # Determine the path
                path = os.path.join(directory_path, "star-fitted_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.save(path)

            # Check if a source was found for this star
            elif star.has_source:

                # Determine the path
                path = os.path.join(directory_path, "star-detected_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.save(path)

            # If no source was found for this star
            else:

                # Determine the path
                path = os.path.join(directory_path, "star-undetected_" + str(star.index) + ".fits")

                # Create a source for the desired sigma level and outer factor
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the cutout as a FITS file
                source.save(path)

    # -----------------------------------------------------------------

    def write_result(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the resulting FITS file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        self.log.info("Writing resulting frame to " + path + " ...")

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
    def have_contour(self):

        """
        This function ...
        :return:
        """

        count = 0
        for star in self.stars: count += star.has_contour
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

    def set_statistics(self):

        """
        This function ...
        :return:
        """

        index_column = []
        have_source_column = []
        have_model_column = []
        have_saturation_column = []

        # Peak
        x_peak_column = []
        y_peak_column = []

        # Fitting -> FWHM
        fwhm_column = []

        # Saturation -> aperture
        x_centroid_column = []
        y_centroid_column = []
        a_column = []
        b_column = []
        angle_column = []

        # Ignore
        ignore_column = []

        # Other
        #not_star_column = []
        #force_column = []
        #not_saturation_column = []

        # Loop over all stars
        for star in self.stars:

            index_column.append(star.index)
            have_source_column.append(star.has_source)
            have_model_column.append(star.has_model)
            have_saturation_column.append(star.has_saturation)

            if star.has_source and star.source.has_peak:

                x_peak_column.append(star.source.peak.x)
                y_peak_column.append(star.source.peak.y)

            else:

                x_peak_column.append(None)
                y_peak_column.append(None)

            fwhm_column.append(star.fwhm if star.has_model else None)

            if star.has_saturation:

                contour_position = star.contour.center
                x_centroid_column.append(contour_position.x)
                y_centroid_column.append(contour_position.y)
                a_column.append(star.contour.major)
                b_column.append(star.contour.minor)
                angle_column.append(star.contour.degree)

            else:

                x_centroid_column.append(None)
                y_centroid_column.append(None)
                a_column.append(None)
                b_column.append(None)
                angle_column.append(None)

            ignore_column.append(star.ignore)
            #not_star_column.append()
            #force_column.append()
            #not_saturation_column.append()

        # Create data structure and set column names
        data = [index_column, have_source_column, have_model_column, have_saturation_column,
                x_peak_column, y_peak_column, fwhm_column, x_centroid_column, y_centroid_column, a_column, b_column,
                angle_column, ignore_column]
        names = ["Star index", "Detected", "Fitted", "Saturated", "Peak x position", "Peak y position", "FWHM",
                 "Aperture x centroid", "Aperture y centroid", "Aperture a length", "Aperture b length",
                 "Aperture angle", "Ignore"]

        # Create the statistics table
        self.statistics = tables.new(data, names)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # If requested, write out the stellar catalog
        if self.config.write_catalog: self.write_catalog()

        # If requested, write out statistics
        if self.config.write_statistics: self.write_statistics()

        # If requested, write out the star region
        if self.config.write_star_region: self.write_star_region()

        # If requested, write out the saturation region
        if self.config.write_saturation_region: self.write_saturation_region()

        # If requested, write out the frame where the stars are masked
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the star cutout boxes
        if self.config.write_cutouts: self.write_cutouts()

        # If requested, write out the result
        if self.config.write_result: self.write_result()

# -----------------------------------------------------------------

def matches_galaxy_position(position, position_list, type_list, encountered, min_distances, distances=None):

    """
    This function ...
    :param position:
    :param galaxy_extractor:
    :param encountered_galaxies:
    :return:
    """

    for j in range(len(encountered)):

        # Ignore already encountered galaxies (an other star is already identified with it)
        if encountered[j]: continue

        # Calculate the pixel position of the galaxy
        galaxy_position = position_list[j]

        # Calculate the distance between the star's position and the galaxy's center
        difference = galaxy_position - position
        distance = difference.norm

        # Add the star-galaxy distance to the list of distances
        if distances is not None: distances.append(distance)

        # The principal galaxy/galaxies
        if type_list[j] == "principal":

            # Check whether the star-galaxy distance is smaller than a certain threshold
            if distance <= min_distances.principal: return True

        # Companion galaxies
        elif type_list[j] == "companion":

            if distance <= min_distances.companion:

                # Indicate that the current star has been identified with the galaxy with index j
                encountered[j] = True
                return True

        # All other galaxies in the frame
        else:

            if distance <= min_distances.other:

                # Indicate that the current star has been identified with the galaxy with index j
                encountered[j] = True
                return True

    # Return False if none of the galaxies provided a match
    return False

# -----------------------------------------------------------------
