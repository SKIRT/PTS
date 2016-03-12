#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.trainedextractor Contains the TrainedExtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import ndimage

# Import astronomical modules
from photutils import detect_sources
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..core.source import Source
from ..basics.mask import Mask
from ..basics.region import Region
from ..basics.geometry import Ellipse
from ..tools import statistics, masks, plotting, general
from ..analysis import sources
from ..train import Classifier
from ..object.star import Star
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log

# -----------------------------------------------------------------

class TrainedFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(TrainedFinder, self).__init__(config, "magic")

        # -- Attributes --

        # The image (with bad and sources masks)
        self.image = None

        # Masks ...
        self.special_mask = None
        self.ignore_mask = None

        # Set the segmentation map to None initially
        self.segments = None

        # List of sources
        self.sources = []

        # The classifier
        self.classifier = Classifier()

        # List of stars
        self.stars = []

        # Local references to the galaxy and star finder
        self.galaxy_finder = None
        self.star_finder = None

        self.region = Region()

    # -----------------------------------------------------------------

    def run(self, image, galaxy_finder=None, star_finder=None, special=None, ignore=None):

        """
        This function ...
        :param image:
        :param galaxy_finder:
        :param star_finder:
        :param special:
        :param ignore:
        :return:
        """

        # 1. Call the setup function
        self.setup(image, galaxy_finder, star_finder, special, ignore)

        # 2. Find sources
        self.find_sources()

        self.create_region()

        # 3. Remove sources
        #if self.config.remove: self.remove_sources()

        # 4. Classify sources
        #if self.config.classify: self.classify_sources()

    # -----------------------------------------------------------------

    def setup(self, image, galaxy_finder=None, star_finder=None, special_mask=None, ignore_mask=None):

        """
        This function ...
        :param image:
        :param galaxy_finder:
        :param star_finder:
        :param special_mask:
        :param ignore_mask:
        :return:
        """

        # Call the setup function of the base class
        super(TrainedFinder, self).setup()

        # Make a local reference to the image
        self.image = image
        self.special_mask = special_mask
        self.ignore_mask = ignore_mask

        # Make local references to the galaxy and star extractors
        self.galaxy_finder = galaxy_finder
        self.star_finder = star_finder

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for sources in the frame ...")

        # Find sources by locating peaks
        if self.config.detection.method == "peaks": self.find_sources_peaks()

        # Find sources by segmenting the image
        elif self.config.detection.method == "segmentation": self.find_sources_segmentation()

        # Find sources by using SExtractor
        elif self.config.detection.method == "sextractor": self.find_sources_sextractor()

        # Unknown source detection method
        else: raise ValueError("Unknown source detection method")

    # -----------------------------------------------------------------

    def remove_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing the sources from the frame ...")

        # Loop over all sources
        for source in self.sources:

            special = self.special_mask.masks(source.center) if self.special_mask is not None else False

            # Estimate the background
            interpolation_method = "local_mean"
            sigma_clip = True
            source.estimate_background(interpolation_method, sigma_clip)

            if special: source.plot(title="Estimated background for source")

            # Replace the frame with the estimated background
            source.background.replace(self.image.frames.primary, where=source.mask)

            if special: plotting.plot_box(self.image.frames.primary[source.cutout.y_slice, source.cutout.x_slice], title="Replaced frame inside this box")

    # -----------------------------------------------------------------

    def create_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing elliptical contours to encompass sources ...")

        # Return the list of apertures
        contours = sources.find_contours(self.segments, self.segments, self.config.detection.apertures.sigma_level)

        # Add shapes to region
        for contour in contours: self.region.append(contour)

    # -----------------------------------------------------------------

    def classify_sources(self):

        """
        This function ...
        :return:
        """

        fwhms = self.star_finder.fwhms
        if len(fwhms) > 0:
            min_fwhm = min(fwhms)
            max_fwhm = max(fwhms)
        else:
            if self.star_finder.config.use_frame_fwhm and self.image.fwhm is not None:
                fwhm = self.image.fwhm.to("arcsec").value / self.image.frames.primary.xy_average_pixelscale.to("arcsec/pix").value
            else: fwhm = self.star_finder.fwhm
            min_fwhm = fwhm * 0.5
            max_fwhm = fwhm * 1.5

        # Loop over all sources
        for source in self.sources:

            special = self.special_mask.masks(source.center) if self.special_mask is not None else False

            # Find peaks
            peaks = source.locate_peaks(3.0)

            if source.has_peak or self.config.classification.fitting.fit_if_undetected:

                # Create sky coordinate from the peak position
                position = SkyCoord.from_pixel(source.peak.x, source.peak.y, self.image.wcs, mode="wcs")

                # Create star object
                index = None
                star = Star(index, catalog=None, id=None, position=position, ra_error=None, dec_error=None)

                # Set the source
                star.source = source

                # Set whether this star should be treated as special
                star.special = special

                # Try to fit star
                star.fit_model(self.config.classification.fitting)

                # If a model is found, add the star to the list
                if star.has_model:

                    if star.fwhm <= max_fwhm and star.fwhm >= min_fwhm:
                        self.stars.append(star)
                    #else: print(star.fwhm)

            # Test whether this source corresponds to a star
            #if self.classifier.is_star(source):

                # Create a Star instance
                #star = Star()

                # Add the star to the list of stars
                #self.stars.append(star)

            # Test whether the source corresponds to a galaxy
            #elif self.classifier.is_galaxy(source):

                # Create a Galaxy instance
                #galaxy = Galaxy()

                # Add the galaxy to the list of galaxies
                #self.galaxies.append(galaxy)

            # Not a star or a galaxy
            #else: log.debug("The source does not correspond to a star or galaxy")

    # -----------------------------------------------------------------

    def remove_contours(self):

        """
        This function ...
        :return:
        """

        # Loop over all other apertures
        for contour in self.contours:

            # Configuration settings
            sigma_clip = self.config.aperture_removal.sigma_clip
            interpolation_method = self.config.aperture_removal.interpolation_method
            expansion_factor = self.config.aperture_removal.expansion_factor

            # Create a source
            ellipse = Ellipse(contour.center, contour.radius * expansion_factor, contour.angle)
            source = Source.from_ellipse(self.image.frames.primary, ellipse, self.config.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background("local_mean", True)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.image.frames.primary, where=source.mask)

    # -----------------------------------------------------------------

    def find_sources_peaks(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def find_sources_segmentation(self):

        """
        This function ...
        :return:
        """

        mask = Mask(self.galaxy_finder.segments) + Mask(self.star_finder.segments)

        frame = self.image.frames.primary.copy()

        frame[mask] = 0.0

        # Create the sigma-clipped mask
        if "bad" in self.image.masks: mask += self.image.masks.bad
        clipped_mask = statistics.sigma_clip_mask(frame, 3.0, mask)

        # Calculate the median sky value and the standard deviation
        median = np.median(np.ma.masked_array(frame, mask=clipped_mask).compressed())
        stddev = np.ma.masked_array(frame, mask=clipped_mask).std()

        # Calculate the detection threshold
        threshold = median + (3.0 * stddev)

        #kernel = self.star_finder.kernel # doesn't work when there was no star extraction on the image, self.star_finder does not have attribute image thus cannot give image.fwhm
        if self.star_finder.config.use_frame_fwhm and self.image.fwhm is not None:

            fwhm = self.image.fwhm.to("arcsec").value / self.image.frames.primary.xy_average_pixelscale.to("arcsec/pix").value
            sigma = fwhm * statistics.fwhm_to_sigma
            kernel = Gaussian2DKernel(sigma)

        else: kernel = self.star_finder.kernel

        try:
            # Create a segmentation map from the frame
            self.segments = Frame(detect_sources(frame, threshold, npixels=5, filter_kernel=kernel).data)
        except RuntimeError:

            log.debug("Runtime error during detect_sources ...")
            log.debug("kernel = " + str(kernel))

            conv_mode = 'constant'
            conv_val = 0.0
            image = ndimage.convolve(self.image.frames.primary, kernel.array, mode=conv_mode, cval=conv_val)

            log.debug("median = " + str(median))
            log.debug("stddev = " + str(stddev))

            #print("image=", image)
            log.debug("image.ndim = " + str(image.ndim))
            log.debug("type image = " + type(image))
            log.debug("image.shape = "+ str(image.shape))
            log.debug("threshold = " + str(threshold))
            image = image > threshold
            log.debug("image.ndim = " + str(image.ndim))
            log.debug("type image = " + str(type(image)))
            log.debug("image.shape = " + str(image.shape))

        # Eliminate the principal galaxy and companion galaxies from the segments
        if self.galaxy_finder is not None:

            # Determine the mask that covers the principal and companion galaxies
            eliminate_mask = self.galaxy_finder.principal_mask + self.galaxy_finder.companion_mask

            # NEW: PLUS: Eliminate the segments covered by the 'ignore mask'
            if self.ignore_mask is not None: eliminate_mask += self.ignore_mask

            # Check where the galaxy mask overlaps with the segmentation map
            overlap = masks.intersection(self.segments, eliminate_mask)
            if np.any(overlap):

                # Check which indices are present in the overlap map
                possible = np.array(range(1, np.max(overlap) + 1))
                present = np.in1d(possible, overlap)
                indices = possible[present]

                # Remove the galaxies from the segmentation map
                for index in indices: self.segments[self.segments == index] = 0

        # Find apertures
        #contours = self.find_contours()

        #self._create_sources_list(contours)

    # -----------------------------------------------------------------

    def _create_sources_list(self, contours):

        """
        This function ...
        :param contours:
        :return:
        """

        # Construct sources
        for contour in contours:

            # Special ...
            special = self.special_mask.masks(contour.center) if self.special_mask is not None else False

            background_factor = 1.5

            # If the aperture has to be rescaled
            #aperture.a *= 1.0
            #aperture.b *= 1.0

            # No: use the FWHM ! Hmm.. or not: saturation ?

            # Create a source from the aperture
            source = Source.from_ellipse(self.image.frames.primary, contour, background_factor)

            if special: source.plot(title="Source created from contour around segment")

            y_min = source.cutout.y_min
            y_max = source.cutout.y_max
            x_min = source.cutout.x_min
            x_max = source.cutout.x_max

            # Create source mask from the segmentation map
            #mask = Mask(self.segments[y_min:y_max, x_min:x_max])

            #label = self.segments[y_min:y_max, x_min:x_max][]

            segments_cutout = self.segments[y_min:y_max, x_min:x_max]

            label = segments_cutout[int(round(0.5*segments_cutout.shape[0])), int(round(0.5*segments_cutout.shape[1]))]

            # If the center pixel is identified as being part of the background (the center does not
            # correspond to a segment), find the different indices that are within the segments_cutout
            if label == 0:

                # Check which indices are present in the overlap map
                #possible = np.array(range(1, np.max(segments_cutout) + 1))
                #present = np.in1d(possible, segments_cutout)
                #indices = possible[present]

                # Loop with a spiral from the center and check which is the first non-zero index that is encountered
                for spiral_x,spiral_y in general.spiral(segments_cutout.shape[1], segments_cutout.shape[0]):

                    current_label = segments_cutout[spiral_y, spiral_x]
                    if current_label != 0.0:
                        label = current_label
                        break

            # If the label is still zero, don't create a source for this contour (skip the code below)
            if label != 0:

                mask = Mask(segments_cutout == label)

                if special: plotting.plot_box(mask, title="Mask created from segment")

                mask = mask.fill_holes()

                if special: plotting.plot_box(mask, title="Filled holes in mask created from segment")

                # Set the source mask
                source.mask = mask

                if special: source.plot(title="Source after mask has been replaced with segment mask")

                # Add the source
                self.sources.append(source)

    # -----------------------------------------------------------------

    def find_sources_sextractor(self):

        """
        This function ...
        :return:
        """

        # Create a SExtractor instance
        sextractor = SExtractor()

        # Run SExtractor on the image frame
        sextractor.run(self.image.frames.primary)

    # -----------------------------------------------------------------

    def write_star_region(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the star region file
        path = self.full_output_path(self.config.writing.star_region_path)

        # Inform the user
        log.info("Writing star region to " + path + " ...")

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Loop over all galaxies
        for star in self.stars:

            # Get the center in pixel coordinates
            center = star.pixel_position(self.image.wcs, "wcs")

            # Determine the color, based on the detection level
            color = "blue"

            #if star.has_source: region.append(star.source.contour, color)
            #else: region.append(star.ellipse())

            # Determine the FWHM
            fwhm = star.fwhm

            # Calculate the radius in pixels
            #radius = fwhm * statistics.fwhm_to_sigma * self.config.removal.sigma_level
            radius = fwhm * statistics.fwhm_to_sigma * 3.0

            # Show a circle for the star
            suffix = " # "
            color_suffix = "color = " + color
            suffix += color_suffix
            print("image;circle({},{},{})".format(center.x+1, center.y+1, radius) + suffix, file=f)

            # Draw a cross for the peak position
            suffix = " # "
            point_suffix = "point = x"
            suffix += point_suffix
            print("image;point({},{})".format(star.source.peak.x+1, star.source.peak.y+1) + suffix, file=f)

        # Close the file
        f.close()

# -----------------------------------------------------------------
