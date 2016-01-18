#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
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
from astropy.coordinates import Angle
from photutils import detect_sources
from astropy.coordinates import SkyCoord

# Import the relevant AstroMagic classes and modules
from ..core import Frame, Source
from ..basics import Position, Extent, Mask, Ellipse
from ..tools import statistics, masks, catalogs
from ..analysis import SExtractor, sources
from ..train import Classifier
from ..sky import Galaxy, Star

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class TrainedExtractor(Configurable):

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
        super(TrainedExtractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image frame and mask
        self.frame = None
        self.input_mask = None

        # Set the segmentation map to None initially
        self.segments = None

        # List of sources
        self.sources = []

        # The classifier
        self.classifier = Classifier()

        # List of stars
        self.stars = []

        # Stellar catalog
        self.catalog = []

        # Mask
        self.mask = None

        # Set the catalog to None initially
        self.catalog = None

    # -----------------------------------------------------------------

    def run(self, frame, input_mask, galaxyextractor=None, starextractor=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, input_mask, galaxyextractor, starextractor)

        # 2. Find sources
        self.find_sources()

        # 3. Remove sources
        if self.config.remove: self.remove_sources()

        # 4. Classify sources
        if self.config.classify: self.classify_sources()

        # 5. Build catalogs
        if self.config.build_catalogs: self.build_catalog()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame, input_mask, galaxyextractor=None, starextractor=None):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(TrainedExtractor, self).setup()

        # Make a local reference to the frame
        self.frame = frame
        self.input_mask = input_mask

        # Create a mask with shape equal to the shape of the frame
        self.mask = Mask.from_shape(self.frame.shape)

        # Make local references to the galaxy and star extractors
        self.galaxy_extractor = galaxyextractor
        self.star_extractor = starextractor

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Looking for sources in the frame ...")

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
        self.log.info("Removing the other sources from the frame ...")

        # Loop over all sources
        for source in self.sources:

            # Estimate the background
            interpolation_method = "local_mean"
            sigma_clip = True
            source.estimate_background(interpolation_method, sigma_clip)

            # Replace the frame with the estimated background
            source.background.replace(self.frame, where=source.mask)

            # Update the mask
            self.mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

    # -----------------------------------------------------------------

    def find_contours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Constructing elliptical contours to encompass sources ...")

        # Return the list of apertures
        return sources.find_contours(self.frame, self.segments, self.config.detection.apertures.sigma_level)

    # -----------------------------------------------------------------

    def classify_sources(self):

        """
        This function ...
        :return:
        """

        fwhms = self.star_extractor.fwhms
        min_fwhm = min(fwhms)
        max_fwhm = max(fwhms)

        # Loop over all sources
        for source in self.sources:

            # Find peaks
            peaks = source.locate_peaks(3.0)

            if source.has_peak or self.config.classification.fitting.fit_if_undetected:

                # Create sky coordinate from the peak position
                position = SkyCoord.from_pixel(source.peak.x, source.peak.y, self.frame.wcs, mode="wcs")

                # Create star object
                index = None
                star = Star(index, catalog="DustPedia", id=None, position=position, ra_error=None, dec_error=None)

                star.source = source

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
            #else: self.log.debug("The source does not correspond to a star or galaxy")

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
            source = Source.from_ellipse(self.frame, ellipse, self.config.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background("local_mean", True)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

            # Update the mask
            self.mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

    # -----------------------------------------------------------------

    def build_catalog(self):

        """
        This function ...
        :return:
        """

        # Create the stellar catalog
        self.catalog = catalogs.from_stars(self.stars)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # If requested, write out the frame where the sources are masked
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, ...
        if self.config.write_star_region: self.write_star_region()

        # If requested, write out the segmentation map
        if self.config.write_segments: self.write_segments()

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

        # Create the sigma-clipped mask
        clipped_mask = statistics.sigma_clip_mask(self.frame, 3.0, self.input_mask)

        # Calculate the median sky value and the standard deviation
        median = np.median(np.ma.masked_array(self.frame, mask=clipped_mask).compressed())
        stddev = np.ma.masked_array(self.frame, mask=clipped_mask).std()

        # Calculate the detection threshold
        threshold = median + (3.0 * stddev)

        kernel = self.star_extractor.kernel

        try:
            # Create a segmentation map from the frame
            self.segments = detect_sources(self.frame, threshold, npixels=5, filter_kernel=kernel).data
        except RuntimeError:

            print("kernel=", kernel)
            #print("self.frame=", self.frame)
            print("self.frame.ndim=", self.frame.ndim)

            conv_mode = 'constant'
            conv_val = 0.0
            image = ndimage.convolve(self.frame, kernel.array, mode=conv_mode, cval=conv_val)

            print("median=", median)
            print("stddev=", stddev)

            #print("image=", image)
            print("image.ndim=", image.ndim)
            print("type image=", type(image))
            print("image.shape=", image.shape)
            print("threshold=", threshold)
            image = image > threshold
            print("image.ndim=", image.ndim)
            print("type image=", type(image))
            print("image.shape=", image.shape)

        # Write the segmentation map to file
        #Frame(self.segments).save(self.config.writing.other_segments_path[:-5]+"_or.fits")

        # Eliminate the principal galaxy and companion galaxies from the segments
        if self.galaxy_extractor is not None:

            # Determine the mask that covers the principal and companion galaxies
            galaxy_mask = self.galaxy_extractor.principal_mask + self.galaxy_extractor.companion_mask

            # Check where the galaxy mask overlaps with the segmentation map
            overlap = masks.intersection(self.segments, galaxy_mask)
            if not np.any(overlap): return

            # Check which indices are present in the overlap map
            possible = np.array(range(1, np.max(overlap) + 1))
            present = np.in1d(possible, overlap)
            indices = possible[present]

            # Remove the galaxies from the segmentation map
            for index in indices: self.segments[self.segments == index] = 0

        # Find apertures
        contours = self.find_contours()

        # Construct sources
        for contour in contours:

            background_factor = 1.5

            # If the aperture has to be rescaled
            #aperture.a *= 1.0
            #aperture.b *= 1.0

            # No: use the FWHM ! Hmm.. or not: saturation ?

            # Create a source from the aperture
            source = Source.from_ellipse(self.frame, contour, background_factor)

            y_min = source.cutout.y_min
            y_max = source.cutout.y_max
            x_min = source.cutout.x_min
            x_max = source.cutout.x_max

            # Create source mask from the segmentation map
            mask = Mask(self.segments[y_min:y_max, x_min:x_max])
            mask = mask.fill_holes()

            # Set the source mask
            source.mask = mask

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
        sextractor.run(self.frame)

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the masked frame file path
        path = self.full_output_path(self.config.writing.masked_frame_path)

        # Inform the user
        self.log.info("Writing masked frame to " + path + " ...")

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(path)

    # -----------------------------------------------------------------

    def write_star_region(self):

        """
        This function ...
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

        # Loop over all galaxies
        for star in self.stars:

            # Get the center in pixel coordinates
            center = star.pixel_position(self.frame.wcs, "wcs")

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
            print("image;circle({},{},{})".format(center.x, center.y, radius) + suffix, file=f)

            # Draw a cross for the peak position
            suffix = " # "
            point_suffix = "point = x"
            suffix += point_suffix
            print("image;point({},{})".format(star.source.peak.x, star.source.peak.y) + suffix, file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the segmentation file
        path = self.full_output_path(self.config.writing.segments_path)

        # Inform the user
        self.log.info("Writing the segmentation map to " + path + " ...")

        # Save the segmentation map
        Frame(self.segments).save(path)

# -----------------------------------------------------------------
