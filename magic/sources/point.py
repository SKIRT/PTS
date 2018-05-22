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
from collections import OrderedDict

# Import astronomical modules
from astropy.modeling.models import Gaussian2D

# Import the relevant PTS classes and modules
from ..basics.vector import Extent
from ..region.list import PixelRegionList
from ..basics.coordinate import SkyCoordinate
from ..region.point import PixelPointRegion
from ..region.circle import PixelCircleRegion
from ..region.ellipse import PixelEllipseRegion
from ..core.frame import Frame
from ..core.detection import Detection
from ..tools import statistics, fitting
from ...core.basics.configurable import Configurable
from ...core.tools import tables, arrays
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ..tools import plotting
from ..basics.stretch import PixelStretch
from ...core.basics.table import SmartTable
from ...core.units.parsing import parse_unit as u
from ..basics.vector import Pixel
from ..core.kernel import ConvolutionKernel

# -----------------------------------------------------------------

class PointSourceTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["RA"] = (float, u("deg"), "right ascension")
    _column_info["DEC"] = (float, u("deg"), "declination")
    _column_info["Detected"] = (bool, None, "Has source detected")
    _column_info["Flux"] = (float, u("Jy"), "flux for the point source")
    _column_info["Flux error"] = (float, u("Jy"), "error on the flux value")
    _column_info["FWHM"] = (float, u("arcsec"), "FWHM of the point source")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSourceTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_source(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        if source is not None:

            # Inform the user
            log.info("Adding source " + str(source.index) + " to the table of point sources ...")

            # Get point source properties
            ra = source.position.ra
            dec = source.position.dec
            detected = source.has_detection
            flux = None
            flux_error = None
            fwhm = source.fwhm

            # Construct the row
            values = [ra, dec, detected, flux, flux_error, fwhm]

        else: values = [None, None, None, None, None, None]

        # Add a row
        self.add_row(values)

    # -----------------------------------------------------------------

    def get_position(self, index):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(ra=self["RA"][index] * self["RA"].unit, dec=self["DEC"][index] * self["DEC"].unit, unit="deg", frame="fk5")

    # -----------------------------------------------------------------

    def is_detected(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Detected"][index]

    # -----------------------------------------------------------------

    def get_flux(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Flux", index)

    # -----------------------------------------------------------------

    def get_flux_error(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_quantity("Flux error", index)

    # -----------------------------------------------------------------

    def get_fwhm(self, index):

        """
        This function ...
        :return:
        """

        return self.get_quantity("FWHM", index)

    # -----------------------------------------------------------------

    def fwhms(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self["FWHM"], unit=unit)
        else: return arrays.array_as_list(self["FWHM"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    @property
    def mean_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhms = self.fwhms(asarray=True)
        not_nan = np.logical_not(np.isnan(fwhms))
        return np.mean(fwhms[not_nan])

    # -----------------------------------------------------------------

    @property
    def median_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhms = self.fwhms(asarray=True)
        not_nan = np.logical_not(np.isnan(fwhms))
        return np.median(fwhms[not_nan])

    # -----------------------------------------------------------------

    @property
    def fwhm_stddev(self):

        """
        This function ...
        :return:
        """

        fwhms = self.fwhms(asarray=True)
        not_nan = np.logical_not(np.isnan(fwhms))
        return np.std(fwhms[not_nan])

# -----------------------------------------------------------------

class PointSourceFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PointSourceFinder, self).__init__(*args, **kwargs)

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

        # The mask of the principal galaxy
        self.principal_mask = None

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # 2. Load the sources from the catalog
        self.load_sources()

        # 3. For each star, find a corresponding source in the image
        if not self.config.weak: self.detect_sources()
        else: self.set_detections()

        # 4. Fit analytical models to the stars
        if not self.has_psf: self.fit_psf()

        # 5. Set the final sources
        if self.has_sources: self.adjust_sources()

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

        # Principal galaxy mask
        self.principal_mask = kwargs.pop("principal_mask", None)

        # Get the galaxy list
        self.galaxies = kwargs.pop("galaxies")

        # Create an empty frame for the segments
        self.segments = Frame.zeros_like(self.frame)

        # Initialize the table
        self.table = PointSourceTable()

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

                # Get corresponding pixel
                pixel = Pixel.for_coordinate(pixel_position)

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
                    #if self.galaxies is not None: star_on_galaxy = self.galaxies.principal.contains(pixel_position)
                    if self.principal_mask is not None: star_on_galaxy = self.principal_mask[pixel.y, pixel.x]
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

    def set_detections(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the detections ...")

        # Loop over all sources
        for source in self.sources:

            # Skip None
            if source is None: continue

            # If this source should be ignored, skip it
            if source.ignore: continue

            # Get the parameters of the circle
            radius = PixelStretch(self.config.detection.initial_radius, self.config.detection.initial_radius)
            ellipse = source.ellipse(self.frame.wcs, radius)

            # Create a source object
            detection = Detection.from_ellipse(self.frame, ellipse, self.config.detection.background_outer_factor)
            source.detection = detection

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
            if source.ignore: continue

            # Check if the star has been detected
            if not source.has_detection and self.config.fitting.fit_if_undetected:

                # Get the parameters of the circle
                ellipse = source.ellipse(self.frame.wcs, self.frame.average_pixelscale, self.config.fitting.initial_radius)

                # Create a detection object
                detection = Detection.from_ellipse(self.frame, ellipse, self.config.fitting.background_outer_factor)

            # No deection
            else: detection = None

            # Find a model, if detection was found
            if source.has_detection or detection is not None: source.fit_model(self.config.fitting, detection)

        # If requested, perform sigma-clipping to the list of FWHM's to filter out outliers
        if self.config.fitting.sigma_clip_fwhms and len(self.fwhms_pix_valid) > 0:

            #print(self.fwhms_pix)
            mean, median, stddev = statistics.sigma_clipped_statistics(self.fwhms_pix_valid, self.config.fitting.fwhm_sigma_level)
            lower = median - self.config.fitting.fwhm_sigma_level * stddev
            upper = median + self.config.fitting.fwhm_sigma_level * stddev

            # Loop over all sources for which a model was found
            for source in self.sources:

                # Skip Nones
                if source is None: continue

                # Ignore sources without model
                if not source.has_model: continue

                # Remove the model if its FWHM is clipped out
                if source.fwhm > upper or source.fwhm < lower: source.model = None

        # Inform the user
        if self.have_detection > 0: log.debug("Found a model for {0} out of {1} stars with a detection ({2:.2f}%)".format(self.have_model, self.have_detection, self.have_model/self.have_detection*100.0))

    # -----------------------------------------------------------------

    @property
    def nsources(self):

        """
        This function ...
        :return:
        """

        return len(self.sources)

    # -----------------------------------------------------------------

    @property
    def has_sources(self):

        """
        This function ...
        :return:
        """

        return self.nsources > 0

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

        # Check whether default FWHM is defined
        if default_fwhm is None:
            log.warning("Could not determine the FWHM (no stars could be modeled to a PSF)")
            return

        # Loop over all sources
        for source in self.sources:

            # If no source
            if source is None: continue

            # If this source should be ignored, skip it
            if source.ignore: continue

            # If this source does not have a source, skip it
            if not source.has_detection: continue

            # Create a detection for the desired sigma level and outer factor
            source.detection = source.detection_at_sigma_level(self.frame, default_fwhm, self.config.source_psf_sigma_level, self.config.source_outer_factor)

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

            # No source
            if source is None: continue

            # Get the center in pixel coordinates
            #center = source.pixel_position(self.frame.wcs)

            # Determine the center position of the detection (center of model if present, otherwise position of the star)
            if source.detection is not None:

                # If the star has been modeled succesfully, use the center position of the model
                # Otherwise, use the source's peak
                if source.psf_model is not None: center = fitting.center(source.psf_model)
                elif source.detection.has_peak: center = source.detection.peak
                else:

                    log.warning("Star source does not have peak")
                    center = source.pixel_position(self.frame.wcs)

            # Calculate the pixel coordinate of the star's position
            else: center = source.pixel_position(self.frame.wcs)

            # Determine the color, based on the detection level
            if source.has_model: color = "blue"
            elif source.has_detection: color = "green"
            else: color = "red"

            # Determine the FWHM
            fwhm = default_fwhm if not source.has_model else source.fwhm

            # Calculate the radius in pixels
            if fwhm is not None: radius = fwhm * statistics.fwhm_to_sigma * self.config.source_psf_sigma_level
            else: radius = 0.

            # Convert the source index to a string
            text = str(source.index)

            # Create meta information
            meta = {"color": color, "text": text, "index": source.index}

            # Create the shape and add it to the region
            shape = PixelCircleRegion(center, radius, meta=meta)
            self.regions.append(shape)

            # Add a position for the peak position
            if source.has_detection and source.detection.has_peak:

                # Create meta information for the position
                meta = {"point": "x"}

                # Create the position and add it to the region
                position = PixelPointRegion(source.detection.peak.x, source.detection.peak.y, meta=meta)
                self.regions.append(position)

    # -----------------------------------------------------------------

    def detect_saturation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for saturated stars ...")

        # Check whether detections are found
        with_detection = self.have_detection
        if with_detection == 0:
            # Stop searching for saturation, but don't break the whole process
            print("Not a single source was found")
            self.config.find_saturation = False
        else:

            # Inform the user on the number of stars that have a detection
            log.debug("Number of stars with detection = " + str(with_detection))

            # Calculate the default FWHM, for the stars for which a model was not found
            default_fwhm = self.fwhm_pix

            # Set the number of stars where saturation was removed to zero initially
            success = 0

            # Create point source mask
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
                if source.has_model: assert source.has_detection

                # Note: DustPedia stars will always get a 'source' during removal (with star.source_at_sigma_level) so star.has_source will already pass

                # If a source was not found for this star, skip it unless the remove_if_undetected flag is enabled
                if not source.has_detection and not self.config.saturation.remove_if_undetected: continue

                # Find a saturation source and remove it from the frame
                if default_fwhm is None:
                    log.warning("Could not determine the FWHM (no stars could be modeled to a PSF)")
                else:
                    source.find_saturation(self.frame, self.config.saturation, default_fwhm, star_mask)
                    success += source.has_saturation

            # Inform the user
            log.debug("Found saturation in " + str(success) + " out of " + str(self.have_detection) + " sources with detection ({0:.2f}%)".format(success / self.have_detection * 100.0))

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

            # No source
            if source is None: continue

            # Skip stars without saturation
            if not source.has_saturation: continue

            # Convert the star index to a string
            text = str(source.index)

            # Get aperture properties
            center = source.contour.center
            semimajor = source.contour.semimajor
            semiminor = source.contour.semiminor
            angle = source.contour.angle

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

            # No source
            if source is None: continue

            # Stars with saturation
            if source.has_saturation:

                # Add the saturation segment to the segmentation map
                self.segments[source.saturation.y_slice, source.saturation.x_slice][source.saturation.mask.data] = source.index

            # Stars without saturation
            else:

                # Skip stars without a detection
                if not source.has_detection: continue

                # Add the star segment to the segmentation map
                self.segments[source.detection.y_slice, source.detection.x_slice][source.detection.mask.data] = source.index

    # -----------------------------------------------------------------

    def create_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the table of point sources ...")

        # Loop over the sources
        for source in self.sources:

            # No source?
            #if source is None: continue

            # Add source
            self.table.add_source(source)

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

        # Inform the user
        log.info("Writing the table ...")

        # Determine path
        path = self.output_path_file("point_sources.dat")

        # Write
        self.table.saveto(path)

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
        directory_path = self.output_path_directory(self.config.writing.cutouts_path)

        # Inform the user
        log.info("Writing cutout boxes to " + directory_path + " ...")

        # Calculate the default FWHM based on the stars that could be fitted
        default_fwhm = self.fwhm_pix

        # Loop over all stars
        for star in self.sources:

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
                source = star.source_at_sigma_level(self.frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

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
                source = star.source_at_sigma_level(self.frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

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
                source = star.source_at_sigma_level(self.frame, default_fwhm, sigma_level, outer_factor, use_default_fwhm=True, shape=shape)

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

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the sources
        for source in self.sources:

            # Skip None
            if source is None: positions.append(None)

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(source.pixel_position(self.frame.wcs))

        # Return the list
        return positions

    # -----------------------------------------------------------------

    @property
    def have_detection(self):

        """
        This function ...
        :return:
        """

        count = 0
        for source in self.sources:
            if source is None: continue
            count += source.has_detection
        return count

    # -----------------------------------------------------------------

    @property
    def have_model(self):

        """
        This function ...
        :return:
        """

        count = 0
        for source in self.sources:
            if source is None: continue
            count += source.has_model
        return count

    # -----------------------------------------------------------------

    @property
    def have_saturation(self):

        """
        This function ...
        :return:
        """

        count = 0
        for source in self.sources:
            if source is None: continue
            count += source.has_saturation
        return count

    # -----------------------------------------------------------------

    @property
    def have_contour(self):

        """
        This function ...
        :return:
        """

        count = 0
        for source in self.sources:
            if source is None: continue
            count += source.has_contour
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

        # Loop over all sources
        for source in self.sources:
            
            if source is None: fwhms.append(None)

            # If the star contains a model, add the fwhm of that model to the list
            elif source.has_model:

                fwhm_pix = source.fwhm
                fwhm_arcsec = fwhm_pix * self.frame.average_pixelscale.to("arcsec")
                fwhms.append(fwhm_arcsec)

            else: fwhms.append(None)

        # Return the list
        return fwhms

    # -----------------------------------------------------------------

    @property
    def fwhms_valid(self):

        """
        This function ...
        :return: 
        """

        return [value for value in self.fwhms if value is not None]

    # -----------------------------------------------------------------

    @property
    def fwhms_pix(self):

        """
        This function ...
        :return:
        """

        return [(fwhm / self.frame.average_pixelscale.to("arcsec")).to("").value if fwhm is not None else None for fwhm in self.fwhms]

    # -----------------------------------------------------------------

    @property
    def fwhms_pix_valid(self):

        """
        This function ...
        :return: 
        """

        return [value for value in self.fwhms_pix if value is not None]

    # -----------------------------------------------------------------

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the fluxes of the sources
        fluxes = []

        # Loop over all sources
        for source in self.sources:

            if source is None: fluxes.append(None)

            # If the star contains a source and the background of this source has been subtracted, calculate the flux
            elif source.has_detection and source.detection.has_background:

                # Add the flux to the list
                fluxes.append(source.flux)

            else: fluxes.append(None)

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

        # Loop over all sources
        for source in self.sources:

            # If no source, skip
            if source is None: continue

            # If the source was not fitted, skip it
            if not source.has_model: continue

            # Determine the amplitude and the position of the center of the model
            amplitude_model = source.model.amplitude
            center = source.detection.cutout.rel_position(fitting.center(source.model))

            # Convert into integers
            x = int(round(center.x))
            y = int(round(center.y))

            # Calculate the value of the source at the model's center position
            amplitude_source = source.detection.subtracted[y, x]

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
        fwhms = self.fwhms_valid
        if len(fwhms) == 0: return None

        fwhm_values = [fwhm.to("arcsec").value for fwhm in fwhms]

        # Determine the default FWHM and return it
        if self.config.fwhm.measure == "max":
            return max(fwhm_values) * u("arcsec") * self.config.fwhm.scale_factor
        elif self.config.fwhm.measure == "mean":
            return np.mean(fwhm_values) * u("arcsec") * self.config.fwhm.scale_factor
        elif self.config.fwhm.measure == "median":
            return np.median(fwhm_values) * u("arcsec") * self.config.fwhm.scale_factor
        else: raise ValueError("Unkown measure for determining the default FWHM")

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """
        return (self.fwhm / self.frame.average_pixelscale.to("arcsec")).value if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def psf(self):

        """
        This function ...
        :return:
        """

        #if self.config.use_frame_fwhm and self.frame.fwhm is not None:

        if self.fwhm_pix is None: return None

        # Create a Gaussian convolution kernel and return it
        sigma = self.fwhm_pix * statistics.fwhm_to_sigma
        model = Gaussian2D(1. / (2 * np.pi * sigma ** 2), 0, sigma, sigma)
        return model

    # -----------------------------------------------------------------

    @property
    def has_psf(self):

        """
        This function ...
        :return: 
        """

        return self.psf is not None

    # -----------------------------------------------------------------

    @property
    def kernel(self):

        """
        This function ...
        :return: 
        """

        # Create the kernel
        kernel = ConvolutionKernel.from_model(self.psf, to_filter=self.frame.psf_filter)

        # Return the kernel
        return kernel

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
        for source in self.sources:

            index_column.append(source.index)
            have_source_column.append(source.has_detection)
            have_model_column.append(source.has_model)
            have_saturation_column.append(source.has_saturation)

            if source.has_detection and source.detection.has_peak:

                x_peak_column.append(source.detection.peak.x)
                y_peak_column.append(source.detection.peak.y)

            else:

                x_peak_column.append(None)
                y_peak_column.append(None)

            fwhm_column.append(source.fwhm if source.has_model else None)

            if source.has_saturation:

                contour_position = source.contour.center
                x_centroid_column.append(contour_position.x)
                y_centroid_column.append(contour_position.y)
                a_column.append(source.contour.semimajor)
                b_column.append(source.contour.semiminor)
                angle_column.append(source.contour.angle.degree)

            else:

                x_centroid_column.append(None)
                y_centroid_column.append(None)
                a_column.append(None)
                b_column.append(None)
                angle_column.append(None)

            ignore_column.append(source.ignore)
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
