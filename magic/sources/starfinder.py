#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.starfinder Contains the StarFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit
from astropy.convolution import Gaussian2DKernel

# Import the relevant PTS classes and modules
from ..basics.vector import Extent
from ..basics.region import Region
from ..basics.geometry import Coordinate, Circle, Ellipse
from ..basics.skygeometry import SkyCoordinate
from ..core.frame import Frame
from ..core.source import Source
from ..object.star import Star
from ..tools import statistics, fitting
from ...core.basics.configurable import Configurable
from ...core.tools import tables, filesystem
from ...core.tools.logging import log

# -----------------------------------------------------------------

class StarFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(StarFinder, self).__init__(config, "magic")

        # -- Attributes --

        # Initialize an empty list for the stars
        self.stars = []

        # The image frame
        self.frame = None

        # The mask covering objects that require special attention
        self.special_mask = None

        # The mask of of pixels that should be ignored
        self.ignore_mask = None

        # The mask of bad pixels
        self.bad_mask = None

        # The stellar catalog
        self.catalog = None

        # The statistics table
        self.statistics = None

        # Reference to the galaxy finder
        self.galaxy_finder = None

        # The segmentation map of stars
        self.segments = None

        # The regions of stars and saturation sources
        self.star_region = None
        self.saturation_region = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_finder, catalog, special=None, ignore=None, bad=None):

        """
        This function ...
        :param frame:
        :param galaxy_finder:
        :param catalog:
        :param special:
        :param ignore:
        :param bad:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_finder, catalog, special, ignore, bad)

        # 2. Find the stars
        self.find_stars()

        # 3. Create the star region
        self.create_star_region()

        # 3. If requested, find and remove saturated stars
        if self.config.find_saturation:

            self.find_saturation()
            self.create_saturation_region()

        # 4. Set the statistics
        self.set_statistics()

        # 5. Create the segmentation map
        self.create_segments()

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_finder, catalog, special_mask=None, ignore_mask=None, bad_mask=None):

        """
        This function ...
        :param frame:
        :param galaxy_finder:
        :param catalog:
        :param special_mask:
        :param ignore_mask:
        :param bad_mask:
        """

        # Call the setup function of the base class
        super(StarFinder, self).setup()

        # Make a local reference to the frame
        self.frame = frame

        self.catalog = catalog

        # Special and ignore masks
        self.special_mask = special_mask
        self.ignore_mask = ignore_mask
        self.bad_mask = bad_mask

        # Make a local reference to the galaxy finder
        self.galaxy_finder = galaxy_finder

        # Create an empty frame for the segments
        self.segments = Frame.zeros_like(self.frame)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the star finder ...")

        # Clear the list of stars
        self.stars = []

        # Clear the image frame
        self.frame = None

    # -----------------------------------------------------------------

    def find_stars(self):

        """
        This function ...
        :return:
        """

        # Load the stars from the stellar catalog
        self.load_stars()

        # For each star, find a corresponding source in the image
        self.find_sources()

        # Fit analytical models to the stars
        if not self.config.use_frame_fwhm or self.frame.fwhm is None: self.fit_stars()

        # Set the final soruces
        self.adjust_sources()

    # -----------------------------------------------------------------

    def load_stars(self):

        """
        This function creates the star list from the star catalog.
        :return:
        """

        # Inform the user
        log.info("Loading the stars from the catalog ...")

        # Copy the list of galaxies, so that we can removed already encounted galaxies (TODO: change this to use
        # an 'encountered' list as well
        encountered_galaxies = [False] * len(self.galaxy_finder.galaxies)

        galaxy_pixel_position_list = []
        galaxy_type_list = []
        for galaxy in self.galaxy_finder.galaxies:

            galaxy_pixel_position_list.append(galaxy.pixel_position(self.frame.wcs))
            if galaxy.principal: galaxy_type_list.append("principal")
            elif galaxy.companion: galaxy_type_list.append("companion")
            else: galaxy_type_list.append("other")

        # Keep track of the distances between the stars and the galaxies
        distances = []

        on_galaxy_column = [False] * len(self.catalog)

        # Create the list of stars
        for i in range(len(self.catalog)):

            # Get the star properties
            catalog = self.catalog["Catalog"][i]
            star_id = self.catalog["Id"][i]
            ra = self.catalog["Right ascension"][i]
            dec = self.catalog["Declination"][i]
            ra_error = self.catalog["Right ascension error"][i] * Unit("mas")
            dec_error = self.catalog["Declination error"][i] * Unit("mas")
            confidence_level = self.catalog["Confidence level"][i]

            # Check for which bands magnitudes are defined
            magnitudes = {}
            magnitude_errors = {}
            for name in self.catalog.colnames:
                if "magnitude" in name:
                    band = name.split(" magnitude")[0]
                    magnitudes[band] = self.catalog[name][i] * Unit("mag")
                    magnitude_errors[band] = self.catalog[name + " error"][i] * Unit("mag")

            # Create a sky coordinate for the star position
            position = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame="fk5")

            # If the stars falls outside of the frame, skip it
            if not self.frame.contains(position): continue

            # Create a star object
            star = Star(i, catalog=catalog, id=star_id, position=position, ra_error=ra_error,
                        dec_error=dec_error, magnitudes=magnitudes, magnitude_errors=magnitude_errors)

            # Get the position of the star in pixel coordinates
            pixel_position = star.pixel_position(self.frame.wcs)

            # -- Checking for foreground or surroudings of galaxy --

            if "On galaxy" in self.catalog.colnames: star_on_galaxy = self.catalog["On galaxy"][i]
            else:

                # Check whether this star is on top of the galaxy, and label it so (by default, star.on_galaxy is False)
                if self.galaxy_finder is not None: star_on_galaxy = self.galaxy_finder.principal.contains(pixel_position)
                else: star_on_galaxy = False
                on_galaxy_column[i] = star_on_galaxy

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
            if self.special_mask is not None: star.special = self.special_mask.masks(pixel_position)
            if self.ignore_mask is not None: star.ignore = self.ignore_mask.masks(pixel_position)

            # If the input mask masks this star's position, skip it (don't add it to the list of stars)
            #if "bad" in self.image.masks and self.image.masks.bad.masks(pixel_position): continue
            if self.bad_mask is not None and self.bad_mask.masks(pixel_position): continue

            # Don't add stars which are indicated as 'not stars'
            if self.config.manual_indices.not_stars is not None and i in self.config.manual_indices.not_stars: continue

            # Add the star to the list
            self.stars.append(star)

        # Add the 'on_galaxy' column to the catalog if necessary
        if "On galaxy" not in self.catalog.colnames: self.catalog["On galaxy"] = on_galaxy_column

        # Inform the user
        if self.config.fetching.cross_reference_with_galaxies: log.debug("10 smallest distances 'star - galaxy': " + ', '.join("{0:.2f}".format(distance) for distance in sorted(distances)[:10]))

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for sources near the star positions ...")

        # Loop over all stars in the list
        for star in self.stars:

            # If this sky object should be ignored, skip it
            if star.ignore: continue

            # Find a source
            try: star.find_source(self.frame, self.config.detection)
            except Exception as e:

                import traceback
                log.error("Error when finding source")
                print(type(e))
                print(e)
                traceback.print_exc()

                if self.config.plot_track_record_if_exception:

                    if star.has_track_record: star.track_record.plot()
                    else: log.warning("Track record is not enabled")

                log.error("Continuing with next source ...")

        # Inform the user
        log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_source, len(self.stars), self.have_source / len(self.stars) * 100.0))

    # -----------------------------------------------------------------

    def fit_stars(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Fitting analytical profiles to the sources ...")

        # Loop over all stars in the list
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Check if the star has a source (has been detected)
            if not star.has_source and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                ellipse = star.ellipse(self.frame.wcs, self.frame.xy_average_pixelscale, self.config.fitting.initial_radius)

                # Create a source object
                source = Source.from_ellipse(self.frame, ellipse, self.config.fitting.background_outer_factor)

            else: source = None

            # Find a model
            if star.has_source or source is not None: star.fit_model(self.config.fitting, source)

        # If requested, perform sigma-clipping to the list of FWHM's to filter out outliers
        if self.config.fitting.sigma_clip_fwhms:

            mean, median, stddev = statistics.sigma_clipped_statistics(self.fwhms_pix, self.config.fitting.fwhm_sigma_level)
            lower = median - self.config.fitting.fwhm_sigma_level * stddev
            upper = median + self.config.fitting.fwhm_sigma_level * stddev

            # Loop over all stars for which a model was found
            for star in self.stars:

                # Ignore stars without model
                if not star.has_model: continue

                # Remove the model if its FWHM is clipped out
                if star.fwhm > upper or star.fwhm < lower: star.model = None

        # Inform the user
        log.debug("Found a model for {0} out of {1} stars with source ({2:.2f}%)".format(self.have_model, self.have_source, self.have_model/self.have_source*100.0))

    # -----------------------------------------------------------------

    def remove_stars(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Removing the stars from the frame ...")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm_pix

        # Inform the user
        log.debug("Default FWHM used when star could not be fitted: {0:.2f} pixels".format(default_fwhm))

        # Loop over all stars in the list
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Remove the star in the frame
            star.remove(self.frame, self.mask, self.config.removal, default_fwhm)

    # -----------------------------------------------------------------

    def adjust_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the star sources to the same sigma level ...")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm_pix

        # Loop over all stars
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If this star does not have a source, skip it
            if not star.has_source: continue

            # Create a source for the desired sigma level and outer factor
            star.source = star.source_at_sigma_level(self.frame, default_fwhm, self.config.source_psf_sigma_level, self.config.source_outer_factor)

    # -----------------------------------------------------------------

    def find_saturation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for saturated stars ...")

        # Check whether sources are found
        with_source = self.have_source
        if with_source == 0: raise RuntimeError("Not a single source was found")

        # Inform the user on the number of stars that have a source
        log.debug("Number of stars with source = " + str(with_source))

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm_pix

        # Set the number of stars where saturation was removed to zero initially
        success = 0

        star_mask = self.star_region.to_mask(self.frame.xsize, self.frame.ysize)

        if self.config.saturation.only_brightest:

            fluxes = sorted(self.get_fluxes(without_background=True))

            # Percentage method
            if self.config.saturation.brightest_method == "percentage":

                # Get the number of fluxes lower than the percentage of highest fluxes
                percentage = self.config.saturation.brightest_level
                fraction = 0.01 * percentage
                count_before = int((1.0-fraction)*len(fluxes))

                # Determine the flux threshold
                flux_threshold = fluxes[count_before-1]

            # Sigma clipping method
            elif self.config.saturation.brightest_method == "sigma clipping":

                # Determine the sigma level
                sigma_level = self.config.saturation.brightest_level

                # Determine the flux threshold
                flux_threshold = statistics.cutoff(fluxes, "sigma_clip", sigma_level)

            # Invalid option
            else: raise ValueError("Brightest method should be 'percentage' or 'sigma clipping'")

        else: flux_threshold = None

        # Loop over all stars
        for star in self.stars:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If a flux threshold is defined
            if flux_threshold is not None:

                # No source, skip right away
                if not star.has_source: continue

                # Determine the flux of this star
                if not star.source.has_background: star.source.estimate_background()
                flux = star.get_flux(without_background=True)

                # Skip this star if its flux is lower than the threshold
                if flux < flux_threshold: continue

            # If a model was not found for this star, skip it unless the remove_if_not_fitted flag is enabled
            if not star.has_model and not self.config.saturation.remove_if_not_fitted: continue
            if star.has_model: assert star.has_source

            # Note: DustPedia stars will always get a 'source' during removal (with star.source_at_sigma_level) so star.has_source will already pass

            # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
            if not star.has_source and not self.config.saturation.remove_if_undetected: continue

            # Find a saturation source and remove it from the frame
            star.find_saturation(self.frame, self.config.saturation, default_fwhm, star_mask)
            success += star.has_saturation

        # Inform the user
        log.debug("Found saturation in " + str(success) + " out of " + str(self.have_source) + " stars with source ({0:.2f}%)".format(success / self.have_source * 100.0))

    # -----------------------------------------------------------------

    def create_star_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating star region ...")

        # Initialize the region
        self.star_region = Region()

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm_pix

        # Loop over all galaxies
        for star in self.stars:

            # Get the center in pixel coordinates
            center = star.pixel_position(self.frame.wcs)

            # Determine the color, based on the detection level
            if star.has_model: color = "blue"
            elif star.has_source: color = "green"
            else: color = "red"

            # Determine the FWHM
            fwhm = default_fwhm if not star.has_model else star.fwhm

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * self.config.source_psf_sigma_level

            # Convert the star index to a string
            text = str(star.index)

            # Create meta information
            meta = {"color": color, "text": text}

            # Create the shape and add it to the region
            shape = Circle(center, radius, meta=meta)
            self.star_region.append(shape)

            # Add a position for the peak position
            if star.has_source and star.source.has_peak:

                # Create meta information for the position
                meta = {"point": "x"}

                # Create the position and add it to the region
                position = Coordinate(star.source.peak.x, star.source.peak.y, meta=meta)
                self.star_region.append(position)

    # -----------------------------------------------------------------

    def create_saturation_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating saturation region ...")

        # Initialize the region
        self.saturation_region = Region()

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

            radius = Extent(major, minor)

            # Create meta information
            meta = {"color": "white", "text": text}

            # Create the ellipse and add it to the region
            ellipse = Ellipse(center, radius, angle, meta=meta)
            self.saturation_region.append(ellipse)

    # -----------------------------------------------------------------

    def create_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the segmentation map for the stars ...")

        # Loop over all stars
        for star in self.stars:

            # Stars with saturation
            if star.has_saturation:

                # Add the saturation segment to the segmentation map
                self.segments[star.saturation.y_slice, star.saturation.x_slice][star.saturation.mask] = star.index

            # Stars without saturation
            else:

                # Skip stars without a source
                if not star.has_source: continue

                # Add the star segment to the segmentation map
                self.segments[star.source.y_slice, star.source.x_slice][star.source.mask] = star.index

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
        log.info("Writing cutout boxes to " + directory_path + " ...")

        # Calculate the default FWHM based on the stars that could be fitted
        default_fwhm = self.fwhm_pix

        # Loop over all stars
        for star in self.stars:

            # -- Saturation sources ---

            # Check if saturation has been detected for this star
            if star.has_saturation:

                # Determine the path
                path = filesystem.join(directory_path, "saturation_" + str(star.index) + ".fits")

                # Save the saturation source as a FITS file
                star.saturation.save(path, origin=self.name)

            # -- PSF sources ---

            # Check if a model has been found for this star
            if star.has_model:

                # Determine the path
                path = filesystem.join(directory_path, "star-fitted_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.save(path, origin=self.name)

            # Check if a source was found for this star
            elif star.has_source:

                # Determine the path
                path = filesystem.join(directory_path, "star-detected_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.save(path, origin=self.name)

            # If no source was found for this star
            else:

                # Determine the path
                path = filesystem.join(directory_path, "star-undetected_" + str(star.index) + ".fits")

                # Create a source for the desired sigma level and outer factor
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the cutout as a FITS file
                source.save(path, origin=self.name)

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
            if star.has_model:
                fwhm_pix = star.fwhm * Unit("pix")
                fwhm_arcsec = fwhm_pix * self.frame.xy_average_pixelscale.to("arcsec/pix")
                fwhms.append(fwhm_arcsec)

        # Return the list
        return fwhms

    # -----------------------------------------------------------------

    @property
    def fwhms_pix(self):

        """
        This function ...
        :return:
        """

        return [(fwhm / self.frame.xy_average_pixelscale.to("arcsec/pix")).to("pix").value for fwhm in self.fwhms]

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

    def get_fluxes(self, without_background=False):

        """
        This function ...
        :param without_background:
        :return:
        """

        # Initialize a list to contain the fluxes of the stars
        fluxes = []

        # Loop over all stars
        for star in self.stars:

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            if star.has_source and star.source.has_background:

                # Add the flux to the list
                fluxes.append(star.get_flux(without_background))

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

        # If requested, always use the FWHM defined by the frame object
        if self.config.use_frame_fwhm and self.frame.fwhm is not None: return self.frame.fwhm.to("arcsec")

        # If the list of FWHM values is empty (the stars were not fitted yet), return None
        fwhms = self.fwhms
        if len(fwhms) == 0: return None

        fwhm_values = [fwhm.to("arcsec").value for fwhm in fwhms]

        # Determine the default FWHM and return it
        if self.config.fwhm.measure == "max":
            return max(fwhm_values) * Unit("arcsec") * self.config.fwhm.scale_factor
        elif self.config.fwhm.measure == "mean":
            return np.mean(fwhm_values) * Unit("arcsec") * self.config.fwhm.scale_factor
        elif self.config.fwhm.measure == "median":
            return np.median(fwhm_values) * Unit("arcsec") * self.config.fwhm.scale_factor
        else: raise ValueError("Unkown measure for determining the default FWHM")

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        return (self.fwhm / self.frame.xy_average_pixelscale.to("arcsec/pix")).value

    # -----------------------------------------------------------------

    @property
    def kernel(self):

        """
        This function ...
        :return:
        """

        # Create a Gaussian convolution kernel and return it
        sigma = self.fwhm_pix * statistics.fwhm_to_sigma
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
                angle_column.append(star.contour.angle.degree)

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

def matches_galaxy_position(position, position_list, type_list, encountered, min_distances, distances=None):

    """
    This function ...
    :param position:
    :param position_list:
    :param type_list:
    :param encountered:
    :param min_distances:
    :param distances:
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
