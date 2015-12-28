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

# Import the relevant AstroMagic classes and modules
from .core import Frame, Source
from .basics import Position, Extent, Mask
from .tools import statistics, masks
from .analysis import SExtractor

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable

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

        # Initialize an empty list for the apertures constructed from the sources
        self.apertures = []

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

        # 3. Find apertures
        if self.config.find_apertures: self.find_apertures()

        # 4. Remove apertures
        if self.config.remove_apertures: self.remove_apertures()

        # 3. Classify the sources
        if self.config.classify: self.classify()

        # 4. Writing
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
        self.mask = Mask(np.zeros_like(self.frame))

        # Make local references to the galaxy and star extractors
        self.galaxy_extractor = galaxyextractor
        self.star_extractor = starextractor

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

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
        self.log.info("Removing the other sources from the frame")

        # Interpolate over the segments
        mask = self.segments > 0
        interpolated = self.frame.interpolated(mask, self.config.removal.interpolation_method)

        # Adapt the frame
        self.frame[mask] = interpolated[mask]

        # Update the mask
        self.mask[mask] = True

    # -----------------------------------------------------------------

    def find_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Constructing elliptical aperture regions to encompass other contaminating sources")

        # Find apertures for the other sources
        #from scipy import ndimage
        from photutils import segment_properties, properties_table
        #from photutils.segmentation import SegmentProperties
        from photutils import EllipticalAperture

        # Get the segment properties
        # Since there is only one segment in the source.mask (the center segment), the props
        # list contains only one entry (one galaxy)
        properties_list = segment_properties(np.asarray(self.frame), self.segments)

        # Below we perform some steps exactly as in the photutils segment_properties function, but we want to
        # avoid calling that function since it also calls _prepare_data which calls _check_units in turn,
        # and because our frames have a 'unit' attribute it crashes for some reason (saying that if one of the
        # [frames, errors, background] has a unit, the other must have too, although we don't provide these errors
        # and background arrays ... More investigation in this should follow.
        #label_ids = np.unique(self.segments[self.segments > 0])
        #filtered_data = None
        #label_slices = ndimage.find_objects(self.segments)
        #properties_list = []
        #for i, label_slice in enumerate(label_slices):
        #    label = i + 1    # consecutive even if some label numbers are missing
        #    # label_slice is None for missing label numbers
        #    if label_slice is None or label not in label_ids:
        #        continue
        #    segm_props = SegmentProperties(
        #        self.frame, self.segments, label, label_slice=label_slice, error=None,
        #        effective_gain=None, mask=None, background=None,
        #        wcs=self.frame.wcs, filtered_data=filtered_data, data_prepared=True)
        #    properties_list.append(segm_props)

        #table = properties_table(properties)
        for properties in properties_list:

            # Obtain the position, orientation and extent
            position = (properties.xcentroid.value, properties.ycentroid.value)
            a = properties.semimajor_axis_sigma.value * self.config.apertures.sigma_level
            b = properties.semiminor_axis_sigma.value * self.config.apertures.sigma_level
            theta = properties.orientation.value

            ellipticity = (a-b)/b

            # Create the aperture
            if ellipticity < self.config.apertures.max_ellipticity: self.apertures.append(EllipticalAperture(position, a, b, theta=theta))

        # Plotting the apertures
        #from astropy.visualization import SqrtStretch
        #from astropy.visualization.mpl_normalize import ImageNormalize
        #import matplotlib.pylab as plt

        #norm = ImageNormalize(stretch=SqrtStretch())

        #plt.figure()

        #plt.imshow(self.segments, origin='lower', cmap='jet')
        #for aperture in self.other_apertures: aperture.plot(color='blue', lw=1.5, alpha=0.5)

        #plt.show()

    # -----------------------------------------------------------------

    def remove_apertures(self):

        """
        This function ...
        :return:
        """

        # Loop over all other apertures
        for aperture in self.apertures:

            # Configuration settings
            sigma_clip = self.config.aperture_removal.sigma_clip
            interpolation_method = self.config.aperture_removal.interpolation_method
            expansion_factor = self.config.aperture_removal.expansion_factor

            # Create a source object
            # Get the parameters of the elliptical aperture
            x_center, y_center = aperture.positions[0]
            center = Position(x=x_center, y=y_center)

            major = aperture.a * expansion_factor
            minor = aperture.b * expansion_factor

            radius = Extent(x=major, y=minor)

            # theta is in radians
            angle = Angle(aperture.theta, u.rad)

            # Create a source
            source = Source(self.frame, center, radius, angle, self.config.aperture_removal.background_outer_factor)

            # Estimate the background for the source
            source.estimate_background("local_mean", True)

            # Replace the frame in the appropriate area with the estimated background
            source.background.replace(self.frame, where=source.mask)

            # Update the mask
            self.mask[source.cutout.y_slice, source.cutout.x_slice] += source.mask

    # -----------------------------------------------------------------

    def classify(self):

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

        # If requested, write out the frame where the sources are masked
        if self.config.write_masked_frame: self.write_masked_frame()

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
            self.segments = detect_sources(self.frame, threshold, npixels=5, filter_kernel=kernel)
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

        # Inform the user
        self.log.info("Writing masked frame to " + self.config.writing.masked_frame_path)

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(self.config.writing.masked_frame_path)

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing the segmentation map to " + self.config.writing.segments_path)

        # Save the segmentation map
        Frame(self.segments).save(self.config.writing.segments_path)

# -----------------------------------------------------------------
