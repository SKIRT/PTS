#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.other Contains the OtherSourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
#from scipy import ndimage
from skimage.feature import match_template

# Import astronomical modules
from photutils import detect_sources
from astropy.coordinates import SkyCoord

# Import the relevant PTS classes and modules
from ..core.frame import Frame
from ..core.source import Source
from ..basics.mask import Mask
from ..tools import statistics, masks, plotting, general
from ..analysis import sources
from ..object.star import Star
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..core.detection import Detection
from ..region.list import PixelRegionList

# -----------------------------------------------------------------

class OtherSourceFinder(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(OtherSourceFinder, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The image frame
        self.frame = None

        # Masks ...
        self.special_mask = None
        self.ignore_mask = None
        self.bad_mask = None

        # Set the segmentation map to None initially
        self.segments = None

        # List of sources
        self.sources = []

        # The classifier
        #self.classifier = Classifier()

        # List of galaxies and stars (from classification)
        #self.galaxies = GalaxyList()
        #self.stars = StarList()

        # The galaxy and star lists
        self.galaxies = None
        self.stars = None

        # The PSF
        self.kernel = None

        # The galaxy and star segments
        self.galaxy_segments = None
        self.star_segments = None

        # The principal mask and companion mask
        self.principal_mask = None
        self.companion_mask = None

        # The region
        self.region = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Find sources
        self.find_sources()

        # 3. Create the region
        self.create_region()

        # 4. Do dilation if requested
        if self.config.dilate: self.dilate_sources()

        # 3. Remove sources
        #if self.config.remove: self.remove_sources()

        # 4. Classify sources
        #if self.config.classify: self.classify_sources()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(OtherSourceFinder, self).setup()

        # Make a local reference to the image frame
        self.frame = kwargs.pop("frame")

        # Masks
        self.special_mask = kwargs.pop("special_mask", None)
        self.ignore_mask = kwargs.pop("ignore_mask", None)
        self.bad_mask = kwargs.pop("bad_mask", None)

        # Make local references to the galaxy and star extractors
        self.galaxies = kwargs.pop("galaxies", None)
        self.stars = kwargs.pop("stars", None)

        # The PSF kernel
        self.kernel = kwargs.pop("kernel", None)

        # The galaxy and star segments
        self.galaxy_segments = kwargs.pop("galaxy_segments", None)
        self.star_segments = kwargs.pop("star_segments", None)

        # Set the principal mask
        self.principal_mask = kwargs.pop("principal_mask", None)
        self.companion_mask = kwargs.pop("companion_mask", None)

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

        # Find sources using experimental method
        elif self.config.detection.method == "experimental": self.find_sources_experimental()

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
            source.background.replace(self.frame, where=source.mask)

            if special: plotting.plot_box(self.frame[source.cutout.y_slice, source.cutout.x_slice], title="Replaced frame inside this box")

    # -----------------------------------------------------------------

    def create_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Constructing elliptical contours to encompass sources ...")

        # Initialize the region
        self.region = PixelRegionList()

        # Return the list of apertures
        contours = sources.find_contours(self.segments._data, self.segments._data, self.config.detection.apertures.sigma_level)

        # Add shapes to region
        for contour in contours: self.region.append(contour)

    # -----------------------------------------------------------------

    def dilate_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("dilating the source segments ...")

        import math

        # Loop over the shapes in the region
        for shape in self.region:

            # Get the integer label for this shape
            label = int(shape.label)

            # Create a source for this label
            source = Source.from_shape(self.frame, shape, self.config.source_outer_factor)

            # Replace the source mask
            segments_cutout = self.segments[source.y_slice, source.x_slice]
            source.mask = Mask(segments_cutout == label).fill_holes()

            ## CODE FOR DILATION (FROM SOURCES MODULE)

            source = source.zoom_out(self.config.dilation_factor, self.segments, keep_original_mask=True)

            mask_area = np.sum(source.mask)
            area_dilation_factor = self.config.dilation_factor ** 2.
            new_area = mask_area * area_dilation_factor

            ## Circular mask approximation

            # ellipse = find_contour(source.mask.astype(float), source.mask)
            # radius = ellipse.radius.norm

            mask_radius = math.sqrt(mask_area / math.pi)
            new_radius = math.sqrt(new_area / math.pi)

            kernel_radius = new_radius - mask_radius

            # Replace mask
            source.mask = source.mask.disk_dilation(radius=kernel_radius)

            ## SET BACK INTO SEGMENTATION MAP

            self.segments[source.y_slice, source.x_slice][source.mask] = label

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
            if self.star_finder.config.use_frame_fwhm and self.frame.fwhm is not None:
                fwhm = self.frame.fwhm.to("arcsec").value / self.frame.average_pixelscale.to("arcsec").value
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
                position = SkyCoord.from_pixel(source.peak.x, source.peak.y, self.frame.wcs, mode="wcs")

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
            source = Detection.from_ellipse(self.frame, ellipse, self.config.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background("local_mean", True)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

    # -----------------------------------------------------------------

    def find_sources_peaks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources based on local peak detection ...")

    # -----------------------------------------------------------------

    def find_sources_segmentation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources based on image segmentation ...")

        # Create mask
        if self.star_segments is not None: mask = Mask(self.galaxy_segments._data) + Mask(self.star_segments._data)
        else: mask = Mask(self.galaxy_segments._data)

        # Mask the data?
        data = self.frame.copy()
        data[mask] = 0.0

        #mask = Mask(self.galaxy_finder.segments)
        #star_mask = Mask(self.star_finder.segments)

        #data = interpolation.in_paint(self.image.frames.primary, star_mask) # Interpolate over stars
        #data[mask] = 0.0 # set galaxies to zero

        # Create the sigma-clipped mask
        if self.bad_mask is not None: mask += self.bad_mask

        clipped_mask = statistics.sigma_clip_mask(data, self.config.detection.segmentation.clipping_sigma_level, mask)

        # Calculate the median sky value and the standard deviation
        median = np.median(np.ma.masked_array(data, mask=clipped_mask).compressed())
        stddev = np.ma.masked_array(data, mask=clipped_mask).std()

        # Calculate the detection threshold
        threshold = median + (self.config.detection.segmentation.sigma_level * stddev)

        #try:
        # Create a segmentation map from the frame
        self.segments = Frame(detect_sources(data, threshold, npixels=5, filter_kernel=self.kernel).data)
        #except RuntimeError as e:

            #log.error("Runtime error during detect_sources ...")

            #print(e)
            #traceback.print_exc()

            #log.debug("kernel = " + str(kernel))

            #conv_mode = 'constant'
            #conv_val = 0.0
            #image = ndimage.convolve(data, kernel.array, mode=conv_mode, cval=conv_val)

            #log.debug("median = " + str(median))
            #log.debug("stddev = " + str(stddev))

            #print("image=", image)
            #log.debug("image.ndim = " + str(image.ndim))
            #log.debug("type image = " + type(image))
            #log.debug("image.shape = "+ str(image.shape))
            #log.debug("threshold = " + str(threshold))
            #image = image > threshold
            #log.debug("image.ndim = " + str(image.ndim))
            #log.debug("type image = " + str(type(image)))
            #log.debug("image.shape = " + str(image.shape))

        # Create eliminate mask
        eliminate_mask = Mask.empty_like(self.frame)

        # Eliminate the principal galaxy and companion galaxies from the segments
        if self.galaxies is not None:

            #principal_mask = self.galaxies.get_principal_mask(self.frame)
            #companion_mask = self.galaxies.get_companion_mask(self.frame)

            # Add mask of principal galaxy
            if self.principal_mask is None: log.warning("Principal mask is not defined")
            else: eliminate_mask += self.principal_mask

            # Add mask of companion galaxy
            if self.companion_mask is None: log.warning("Companion mask is not defined")
            else: eliminate_mask += self.companion_mask

            # Determine the mask that covers the principal and companion galaxies
            #eliminate_mask += principal_mask + companion_mask

        # NEW: PLUS: Eliminate the segments covered by the 'ignore mask'
        if self.ignore_mask is not None: eliminate_mask += self.ignore_mask

        # NEW: Eliminate the segments covered by the 'bad mask'
        if self.bad_mask is not None: eliminate_mask += self.bad_mask

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
            source = Detection.from_ellipse(self.frame, contour, background_factor)

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

        # Inform the user
        log.info("Finding sources using SExtractor ...")

        # Create a SExtractor instance
        sextractor = SExtractor()

        # Run SExtractor on the image frame
        sextractor.run(self.frame)

    # -----------------------------------------------------------------

    def find_sources_experimental(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources based on template matching ...")

        kernel = self.kernel._array

        result = match_template(self.frame, kernel, pad_input=True)

        #plotting.plot_box(result)

        #y, x = np.unravel_index(np.argmax(result), result.shape)
        #source = Source.around_coordinate(self.frame, Coordinate(x,y), radius=5, factor=1.3)
        #source.plot()

        #from ...core.basics.distribution import Distribution
        #distribution = Distribution.from_values(result.flatten())
        #distribution.plot()

        from photutils import find_peaks

        mask = Mask(self.galaxy_segments) + Mask(self.star_segments)

        peaks = find_peaks(result, 0.8, box_size=5, mask=mask)

        index = 1

        self.segments = Frame.zeros_like(self.frame)

        # Loop over the peaks
        for peak in peaks:

            # Calculate the absolute x and y coordinate of the peak
            x = peak['x_peak']
            y = peak['y_peak']
            coordinate = Coordinate(x,y)

            source = Detection.around_coordinate(self.frame, coordinate, radius=5, factor=1.3)

            self.segments[source.y_slice, source.x_slice][source.mask] = index

            index += 1

            #self.sources.append(source)

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
            print("image;circle({},{},{})".format(center.x+1, center.y+1, radius) + suffix, file=f)

            # Draw a cross for the peak position
            suffix = " # "
            point_suffix = "point = x"
            suffix += point_suffix
            print("image;point({},{})".format(star.source.peak.x+1, star.source.peak.y+1) + suffix, file=f)

        # Close the file
        f.close()

# -----------------------------------------------------------------
