#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.point Contains the PointSourceFinder class.

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
from ..region.list import PixelRegionList
from ..basics.coordinate import PixelCoordinate, SkyCoordinate
from ..region.point import PixelPointRegion
from ..region.circle import PixelCircleRegion
from ..region.ellipse import PixelEllipseRegion
from ..core.frame import Frame
from ..core.source import Source
from ..object.star import Star
from ..tools import statistics, fitting
from ...core.basics.configurable import Configurable
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ..tools import plotting
from .list import StarList
from ...core.basics.map import Map
from ..basics.stretch import PixelStretch
from ...core.basics.table import SmartTable

# -----------------------------------------------------------------

class PointSourceTable(SmartTable):

    """
    This class ...
    """

    column_info = [("RA", float, "deg", "right ascension"),
                   ("DEC", float, "deg", "declination"),
                   ("Flux", float, "Jy", "flux for the point source"),
                   ("Flux error", float, "Jy", "error on the flux value"),
                   ("FWHM", float, "arcsec", "FWHM of the point source")]

    # -----------------------------------------------------------------

    def add_source(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        # Values for the row
        values = []

        # Get point source properties

        # Add a row
        self.add_row(values)

    # -----------------------------------------------------------------

    @property
    def mean_fwhm(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def median_fwhm(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def fwhm_stddev(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

class PointSourceFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(PointSourceFinder, self).__init__(config)

        # -- Attributes --

        # Initialize the list of sources
        self.sources = []

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

        # The galaxy list
        self.galaxies = None

        # The segmentation map of stars
        self.segments = None

        # The regions
        self.regions = None
        self.saturation_regions = None

        # The point source table
        self.table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the sources from the catalog
        self.load_sources()

        # 3. For each star, find a corresponding source in the image
        self.detect_sources()

        # 4. Fit analytical models to the stars
        if not self.config.use_frame_fwhm or self.frame.fwhm is None: self.fit_psf()

        # 5. Set the final sources
        self.adjust_sources()

        # 6. Create the region list
        self.create_regions()

        # 7. If requested, find saturated and diffracted point sources
        if self.config.find_saturation: self.detect_saturation()

        # 8. Create region list of saturated sources
        if self.config.find_saturation: self.create_saturation_region()

        # 9. Create the segmentation map
        self.create_segments()

        # 10. Create the table
        self.create_table()

        # 11. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the setup function of the base class
        super(PointSourceFinder, self).setup()

        # Make a local reference to the frame
        self.frame = kwargs.pop("frame")

        # Make a local reference to the catalog
        self.catalog = kwargs.pop("catalog")

        # Special and ignore masks
        self.special_mask = kwargs.pop("special_mask", None)
        self.ignore_mask = kwargs.pop("ignore_mask", None)
        self.bad_mask = kwargs.pop("bad_mask", None)

        # Get the galaxy list
        self.galaxies = kwargs.pop("galaxies")

        # Create an empty frame for the segments
        self.segments = Frame.zeros_like(self.frame)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the point source finder ...")

        # Create a new list of sources
        self.sources = []

        # Clear the image frame
        self.frame = None

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function creates the source list from the catalog
        :return:
        """

        # Inform the user
        log.info("Loading the sources from the catalog ...")

        # Copy the list of galaxies, so that we can removed already encounted galaxies
        # (TODO: change this to use an 'encountered' list as well
        encountered_galaxies = [False] * len(self.galaxies)

        galaxy_pixel_position_list = []
        galaxy_type_list = []
        for galaxy in self.galaxies:

            galaxy_pixel_position_list.append(galaxy.pixel_position(self.frame.wcs))
            if galaxy.principal: galaxy_type_list.append("principal")
            elif galaxy.companion: galaxy_type_list.append("companion")
            else: galaxy_type_list.append("other")

        # Keep track of the distances between the stars and the galaxies
        distances = []

        on_galaxy_column = [False] * len(self.catalog)

        # Loop over the sources in the catalog
        for index in range(len(self.catalog)):

            # Get position
            position = self.catalog.get_position(index)

            # If the stars falls outside of the frame, take None for the source
            if not self.frame.contains(position): source = None

            else:

                # Calculate the pixel position of the galaxy in the frame
                pixel_position = position.to_pixel(self.frame.wcs)

                # Create a source
                source = self.catalog.create_source(index)

                # Check whether 'special'
                special = self.special_mask.masks(pixel_position) if self.special_mask is not None else False
                cutout = self.frame.cutout_around(pixel_position, 15) if special else None

                # Check whether 'ignore'
                ignore = self.ignore_mask.masks(pixel_position) if self.ignore_mask is not None else False

                # Set attributes based on masks (special and ignore)
                source.special = special
                source.ignore = ignore

                special = False

                # -- Checking for foreground or surroudings of galaxy --

                if "On galaxy" in self.catalog.colnames: star_on_galaxy = self.catalog["On galaxy"][index]
                else:

                    # Check whether this star is on top of the galaxy, and label it so (by default, star.on_galaxy is False)
                    if self.galaxies is not None: star_on_galaxy = self.galaxies.principal.contains(pixel_position)
                    else: star_on_galaxy = False
                    on_galaxy_column[index] = star_on_galaxy

                    if special: plotting.plot_box(cutout, title="On galaxy" if star_on_galaxy else "Not on galaxy")

                # -- Cross-referencing with the galaxies in the frame --

                # Loop over all galaxies to cross-referenc
                if self.config.fetching.cross_reference_with_galaxies and star_on_galaxy:

                    # If a match is found with one of the galaxies, skip this star
                    if matches_galaxy_position(pixel_position, galaxy_pixel_position_list, galaxy_type_list, encountered_galaxies, self.config.fetching.min_distance_from_galaxy, distances):

                        if special: plotting.plot_box(cutout, "Matches galaxy position (distance < " + str(self.config.fetching.min_distance_from_galaxy) + ")")
                        source = None

                if source is not None:

                    # Set other attributes
                    source.on_galaxy = star_on_galaxy

                    # Enable track record if requested
                    #if self.config.track_record: star.enable_track_record()

                    # If the input mask masks this star's position, skip it (don't add it to the list of stars)
                    #if "bad" in self.image.masks and self.image.masks.bad.masks(pixel_position): continue
                    if self.bad_mask is not None and self.bad_mask.masks(pixel_position):

                        if special: plotting.plot_box(cutout, "Covered by bad mask")
                        source = None

                    # Don't add stars which are indicated as 'not stars'
                    #if self.config.manual_indices.not_stars is not None and i in self.config.manual_indices.not_stars:
                    #    if special: plotting.plot_box(cutout, "Indicated as 'not a star'")
                    #    continue

            # Add the source to the list
            self.sources.append(source)

        # Add the 'on_galaxy' column to the catalog if necessary
        #if "On galaxy" not in self.catalog.colnames: self.catalog["On galaxy"] = on_galaxy_column

        # Inform the user
        if self.config.fetching.cross_reference_with_galaxies: log.debug("10 smallest distances 'point source - extended source': " + ', '.join("{0:.2f}".format(distance) for distance in sorted(distances)[:10]))

    # -----------------------------------------------------------------

    def detect_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Detecting the sources ...")

        # Loop over all sources
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this sky object should be ignored, skip it
            if source.ignore: continue

            # Find a source
            try: source.detect(self.frame, self.config.detection)
            except Exception as e:

                import traceback
                log.error("Error when finding source")
                print(type(e))
                print(e)
                traceback.print_exc()

                #if self.config.plot_track_record_if_exception:
                    #if source.has_track_record: star.track_record.plot()
                    #else: log.warning("Track record is not enabled")

                log.error("Continuing with next source ...")

        # Inform the user
        log.debug("Found a source for {0} out of {1} objects ({2:.2f}%)".format(self.have_detection, len(self.sources), self.have_detection / len(self.sources) * 100.0))

    # -----------------------------------------------------------------

    def fit_psf(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Fitting PSF profiles to the point sources ...")

        # Loop over all sources in the list
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this star should be ignored, skip it
            if star.ignore: continue

            # Check if the star has a source (has been detected)
            if not star.has_source and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                ellipse = star.ellipse(self.frame.wcs, self.frame.average_pixelscale, self.config.fitting.initial_radius)

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

    def adjust_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the source detections to the same sigma level ...")

        # Calculate the default FWHM, for the stars for which a model was not found
        default_fwhm = self.fwhm_pix

        # Loop over all sources
        for star in self.sources:

            # If this star should be ignored, skip it
            if star.ignore: continue

            # If this star does not have a source, skip it
            if not star.has_source: continue

            # Create a source for the desired sigma level and outer factor
            star.source = star.source_at_sigma_level(self.frame, default_fwhm, self.config.source_psf_sigma_level, self.config.source_outer_factor)

    # -----------------------------------------------------------------

    def create_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating point source regions ...")

        # Initialize the region
        self.regions = PixelRegionList()

        # Calculate the default FWHM (calculated based on fitted stars)
        default_fwhm = self.fwhm_pix

        # Loop over all sources
        for source in self.sources:

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
            shape = PixelCircleRegion(center, radius, meta=meta)
            self.regions.append(shape)

            # Add a position for the peak position
            if star.has_source and star.source.has_peak:

                # Create meta information for the position
                meta = {"point": "x"}

                # Create the position and add it to the region
                position = PixelPointRegion(star.source.peak.x, star.source.peak.y, meta=meta)
                self.regions.append(position)

    # -----------------------------------------------------------------

    def detect_saturation(self):

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

        # Create star mask
        star_mask = self.regions.to_mask(self.frame.xsize, self.frame.ysize)

        # Only brightest method
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

        # Otherwise, no flux threshold
        else: flux_threshold = None

        # Loop over all sources
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this star should be ignored, skip it
            if source.ignore: continue

            # If a flux threshold is defined
            if flux_threshold is not None:

                # No source, skip right away
                if not source.has_detection: continue

                # Determine the flux of this star
                if not source.detection.has_background: source.detection.estimate_background()
                flux = source.detection.get_flux(without_background=True)

                # Skip this star if its flux is lower than the threshold
                if flux < flux_threshold: continue

            # If a model was not found for this star, skip it unless the remove_if_not_fitted flag is enabled
            if not source.has_model and not self.config.saturation.remove_if_not_fitted: continue
            if star.has_model: assert star.has_source

            # Note: DustPedia stars will always get a 'source' during removal (with star.source_at_sigma_level) so star.has_source will already pass

            # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
            if not star.has_source and not self.config.saturation.remove_if_undetected: continue

            # Find a saturation source and remove it from the frame
            star.find_saturation(self.frame, self.config.saturation, default_fwhm, star_mask)
            success += star.has_saturation

        # Inform the user
        log.debug("Found saturation in " + str(success) + " out of " + str(self.have_source) + " sources with detection ({0:.2f}%)".format(success / self.have_source * 100.0))

    # -----------------------------------------------------------------

    def create_saturation_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating saturation region ...")

        # Initialize the region
        self.saturation_regions = PixelRegionList()

        # Loop over all sources
        for source in self.sources:

            # Skip stars without saturation
            if not star.has_saturation: continue

            # Convert the star index to a string
            text = str(star.index)

            # Get aperture properties
            center = star.contour.center
            semimajor = star.contour.semimajor
            semiminor = star.contour.semiminor
            angle = star.contour.angle

            radius = PixelStretch(semimajor, semiminor)

            # Create meta information
            meta = {"color": "white", "text": text}

            # Create the ellipse and add it to the region
            ellipse = PixelEllipseRegion(center, radius, angle, meta=meta)
            self.saturation_regions.append(ellipse)

    # -----------------------------------------------------------------

    def create_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the segmentation map ...")

        # Loop over all sources
        for source in self.sources:

            # Stars with saturation
            if source.has_saturation:

                # Add the saturation segment to the segmentation map
                self.segments[star.saturation.y_slice, star.saturation.x_slice][star.saturation.mask] = star.index

            # Stars without saturation
            else:

                # Skip stars without a source
                if not star.has_source: continue

                # Add the star segment to the segmentation map
                self.segments[star.source.y_slice, star.source.x_slice][star.source.mask] = star.index

    # -----------------------------------------------------------------

    def create_table(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the table
        self.write_table()

        # Write the region lists
        self.write_regions()

        # Write the segmentation maps
        self.write_segments()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        pass

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
                path = fs.join(directory_path, "saturation_" + str(star.index) + ".fits")

                # Save the saturation source as a FITS file
                star.saturation.saveto(path, origin=self.name)

            # -- PSF sources ---

            # Check if a model has been found for this star
            if star.has_model:

                # Determine the path
                path = fs.join(directory_path, "star-fitted_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.saveto(path, origin=self.name)

            # Check if a source was found for this star
            elif star.has_source:

                # Determine the path
                path = fs.join(directory_path, "star-detected_" + str(star.index) + ".fits")

                # Create source
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the source as a FITS file
                source.saveto(path, origin=self.name)

            # If no source was found for this star
            else:

                # Determine the path
                path = fs.join(directory_path, "star-undetected_" + str(star.index) + ".fits")

                # Create a source for the desired sigma level and outer factor
                source = star.source_at_sigma_level(self.original_frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

                # Estimate the background
                sigma_clip = not star.on_galaxy
                source.estimate_background(method, sigma_clip)

                # Save the cutout as a FITS file
                source.saveto(path, origin=self.name)

    # -----------------------------------------------------------------

    @property
    def positions(self):

        """
        This function ...
        :return:
        """

        return self.stars.get_positions(self.frame.wcs)

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
                fwhm_arcsec = fwhm_pix * self.frame.average_pixelscale.to("arcsec/pix")
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

        return [(fwhm / self.frame.average_pixelscale.to("arcsec/pix")).to("pix").value for fwhm in self.fwhms]

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

        return (self.fwhm / self.frame.average_pixelscale.to("arcsec/pix")).value

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

    def get_statistics(self):

        """
        This function ...
        :return:
        """

        # Create the statistics
        statistics = Map()
        statistics.fwhm = self.fwhm
        return statistics

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
                a_column.append(star.contour.semimajor)
                b_column.append(star.contour.semiminor)
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
