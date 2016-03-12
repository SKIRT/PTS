#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.extractor Contains the SourceExtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.mask import Mask
from ..core.source import Source
from ..tools import interpolation
from ...core.tools.logging import log
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class SourceExtractor(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(SourceExtractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image frame
        self.frame = None

        # The output path
        self.output_path = None

        # The mask of nans
        self.nan_mask = None

        # Regions
        self.galaxy_region = None
        self.star_region = None
        self.saturation_region = None
        self.other_region = None

        # Segmentation maps
        self.galaxy_segments = None
        self.star_segments = None
        self.other_segments = None

        # The total mask of removed sources
        self.mask = None

        self.sources = []

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SourceExtractor instance
        extractor = cls()

        # Return the extractor
        return extractor

    # -----------------------------------------------------------------

    def run(self, frame, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments):

        """
        This function ...
        :param frame:
        :param galaxy_region:
        :param star_region:
        :param saturation_region:
        :param other_region:
        :param galaxy_segments:
        :param star_segments:
        :param other_segments:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments)

        # 2. Create the mask
        #self.create_mask()

        # 2. Load the sources
        self.load_sources()

        # 3. Remove the sources
        self.remove_sources()

        # 4. Set nans back into the frame
        self.set_nans()

    # -----------------------------------------------------------------

    def setup(self, frame, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments):

        """
        This function ...
        :param frame:
        :param galaxy_region:
        :param star_region:
        :param saturation_region:
        :param other_region:
        :param star_segments:
        :param other_segments:
        :return:
        """

        # Set the image frame
        self.frame = frame

        # Regions
        self.galaxy_region = galaxy_region
        self.star_region = star_region
        self.saturation_region = saturation_region
        self.other_region = other_region

        # Segmentation maps
        self.galaxy_segments = galaxy_segments
        self.star_segments = star_segments
        self.other_segments = other_segments

        # Initialize the mask
        self.mask = Mask.empty_like(self.frame)

        # Create a mask of the pixels that are NaNs
        self.nan_mask = Mask.is_nan(self.frame)
        self.frame[self.nan_mask] = 0.0

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the mask of sources to extract ...")

        # Debugging info
        log.debug("Adding galaxies to the mask ...")

        # Add the galaxies to the mask
        self.mask += self.galaxy_segments == 3
        if self.config.remove_companions: self.mask += self.galaxy_segments == 2

        # Debugging info
        log.debug("Adding stars (including saturation) to the mask ...")

        # Initialize a list for the star indices that are present in the star region
        star_indices = set()

        # Loop over all stars in the region
        for shape in self.star_region:

            # Ignore shapes without text, these should be just the positions of the peaks
            if "text" not in shape.meta: continue

            # Ignore shapes with color red (stars without source)
            if shape.meta["color"] == "red": continue

            # Get the star index
            index = int(shape.meta["text"])
            star_indices.add(index)

            # Get the star position
            #position = shape.center
            # Check whether the star is a foreground star
            #if self.principal_mask.masks(position):
                # FOREGROUND STARS: create sources, because the background is estimated by fitting a polynomial to
                # pixels in annulus around shape
                #source_mask = star_segments == index # not necessary, the source is created from the ellipse so
                # the star mask is created again
                # Create a source from the shape
                #source = Source.from_ellipse(self.image.frames.primary, shape, 1.3)
                #source.estimate_background("polynomial", sigma_clip=True)
                # Add the source to the dictionary, with a key that is the star index (so that the source for this star
                # can be replaced by the saturation source, if any
                #self.foreground_stars_sources[index] = source
            # Not a foreground star
            #else: self.other_stars_mask += star_segments == index

            # Create a mask from the shape and add it to the total mask
            self.mask += shape.to_mask(self.frame.xsize, self.frame.ysize)

        # Add the saturation sources
        # Loop over the shapes in the saturation region
        for shape in self.saturation_region:

            # Shapes without text are drawn by the user, always add them to the mask
            if "text" not in shape.meta: self.mask += shape.to_mask(self.frame.xsize, self.frame.ysize)

            else:

                # Get the star index
                index = int(shape.meta["text"])

                # If this star index is not in the star_indices list (the star is removed from the star region by the user),
                # ignore it (don't add the saturation mask for it)
                if index not in star_indices: continue

                # Get the star position
                #position = shape.center
                # Check whether the star is a foreground star
                #if self.principal_mask.masks(position):
                    # Create a source from the shape
                    #source = Source.from_ellipse(self.image.frames.primary, shape, 1.3)
                    # Set the source mask
                    #source.mask = saturation_segments == index
                    # Estimate the background
                    #source.estimate_background("polynomial", True)
                    # Replace the star source by the saturation source
                    #self.foreground_stars_sources[index] = source
                # Not a foreground star
                #else: self.other_stars_mask += saturation_segments == index

                # Add the segment to the mask
                self.mask += self.star_segments == index

        # Debugging info
        log.debug("Adding other sources to the mask ...")

        # Add all segments to the mask
        self.mask += Mask(self.other_segments)

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function ...
        :return:
        """

        # Load the galaxy sources
        self.load_galaxy_sources()

        # Load the star sources
        self.load_star_sources()

        # Load the other sources
        self.load_other_sources()

    # -----------------------------------------------------------------

    def load_galaxy_sources(self):

        """
        This function ...
        :return:
        """

        # Loop over the shapes in the galaxy region
        for shape in self.galaxy_region:

            # Shapes without text are in this case just coordinates
            if "text" not in shape.meta: continue

            # Get the coordinate of the center for this galaxy
            center = shape.center

            # Check the label of the corresponding segment
            label = self.galaxy_segments[int(center.y), int(center.x)]

            if label == 3 or (label == 2 and self.config.remove_companions):

                # Create a source and add it to the list
                source = Source.from_shape(self.frame, shape, 1.3)
                self.sources.append(source)

    # -----------------------------------------------------------------

    def load_star_sources(self):

        """
        This function ...
        :return:
        """

        # Loop over all stars in the region
        for shape in self.star_region:

            # Ignore shapes without text, these should be just the positions of the peaks
            if "text" not in shape.meta: continue

            # Ignore shapes with color red (stars without source)
            if shape.meta["color"] == "red": continue

            # Get the star index
            index = int(shape.meta["text"])

            # Look whether a saturation source is present
            saturation_source = None

            # Add the saturation sources
            # Loop over the shapes in the saturation region
            for j in range(len(self.saturation_region)):

                saturation_shape = self.saturation_region[j]

                if "text" not in saturation_shape.meta: continue

                saturation_index = int(saturation_shape.meta["text"])

                if index != saturation_index: continue
                else:
                    # Remove the saturation shape from the region
                    saturation_shape = self.saturation_region.pop(j)

                    # Create saturation source
                    saturation_source = Source.from_shape(self.frame, saturation_shape, 1.3)

                    # Replace the saturation mask
                    segments_cutout = self.star_segments[saturation_source.y_slice, saturation_source.x_slice]
                    saturation_source.mask = Mask(segments_cutout == index)

                    # Break the loop
                    break

            if saturation_source is not None:

                self.sources.append(saturation_source)

            else:

                source = Source.from_shape(self.frame, shape, 1.3)

                self.sources.append(source)

        # Remove remaining shapes: drawn by the user, always remove them
        #for shape in self.saturation_region:

            # Shapes without text are drawn by the user, always add them to the mask
            #self.mask += shape.to_mask(self.frame.xsize, self.frame.ysize)

            # Create source
            #source = Source.from_shape(self.frame, shape, 1.3)

            # Add the source to the list
            #self.sources.append(source)

    # -----------------------------------------------------------------

    def load_other_sources(self):

        """
        This function ...
        :return:
        """

        # Loop over the shapes in the other sources region
        for shape in self.other_region:

            label = int(shape.meta["text"])

            # Create a source
            source = Source.from_shape(self.frame, shape, 1.3)

            # Replace the source mask
            segments_cutout = self.other_segments[source.y_slice, source.x_slice]
            source.mask = Mask(segments_cutout == label)

            # Add the source to the list
            self.sources.append(source)

    # -----------------------------------------------------------------

    def remove_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the frame over the masked pixels ...")

        # Interpolate
        #self.frame[:] = interpolation.inpaint_biharmonic(self.frame, self.mask)

        nsources = len(self.sources)
        count = 0

        for source in self.sources:

            # Debugging
            log.debug("Estimating background and replacing the frame pixels of source " + str(count) + " of " + str(nsources) + " ...")

            # Estimate the background
            source.estimate_background("biharmonic", True)

            # Replace the pixels by the background
            source.background.replace(self.frame, where=source.mask)

            count += 1

    # -----------------------------------------------------------------

    def set_nans(self):

        """
        This function ...
        :return:
        """

        # Set the NaN pixels to zero in the frame
        self.frame[self.nan_mask] = float("nan")

# -----------------------------------------------------------------
