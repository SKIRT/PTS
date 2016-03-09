#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.extraction.simpleextraction Contains the SimpleExtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant AstroMagic classes and modules
from ..core import Frame, Source, Image
from ..basics import Mask, Region
from ..tools import masks

# Import the relevant PTS classes and modules
from ...core.tools import filesystem, tables
from ...core.tools.logging import log

# -----------------------------------------------------------------

class SimpleExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :param config:
        :return:
        """

        # -- Attributes --

        # The image
        self.image = None

        # The output path
        self.output_path = None

        # The regions
        self.galaxy_region = None
        self.star_region = None
        self.saturation_region = None

    # -----------------------------------------------------------------

    def run(self, image_path, output_path):

        """
        This function ...
        :param image_path:
        :param output_path:
        :return:
        """

        # 1. Call the setup function
        self.setup(image_path, output_path)

        # 2. Load the regions
        self.load_regions()

        # 3. Remove the galaxies
        self.remove_galaxies()

        # 4. Remove the stars
        self.remove_stars()

        # 5. Remove the saturation
        self.remove_saturation()

    # -----------------------------------------------------------------

    def setup(self, image_path, output_path):

        """
        This function ...
        :param image_path:
        :param output_path:
        :return:
        """

        # Set the image
        self.image = Image.from_file(image_path)

        # Set the output path
        self.output_path = output_path

    # -----------------------------------------------------------------

    def load_regions(self):

        """
        This function ...
        :return:
        """

        self.load_galaxy_region()

        self.load_star_region()

        self.load_saturation_region()

    # -----------------------------------------------------------------

    def load_galaxy_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy region ...")

        galaxy_region_path = filesystem.join(self.output_path, "galaxies.reg")
        self.galaxy_region = Region.from_file(galaxy_region_path)

        #log.info("Loading the galaxy mask ...")

        #galaxy_mask_path = filesystem.join(self.output_path, "galaxies.")

    # -----------------------------------------------------------------

    def load_star_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the star region ...")

        # Star region
        star_region_path = filesystem.join(self.output_path, "stars.reg")
        self.star_region = Region.from_file(star_region_path)

        # Star segmentation map
        star_segments_path = filesystem.join(self.output_path, "star_segments.fits")
        self.star_segments = Frame.from_file(star_segments_path)

    # -----------------------------------------------------------------

    def load_saturation_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the saturation region")

        # Saturation region
        saturation_region_path = filesystem.join(self.output_path, "saturation.reg")
        self.saturation_region = Region.from_file(saturation_region_path)

        # Saturation segmentation map
        saturation_segments_path = filesystem.join(self.output_path, "saturation_segments.fits")
        self.saturation_segments = Frame.from_file(saturation_region_path)

    # -----------------------------------------------------------------

    def remove_galaxies(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a mask from the galaxy region ...")

        # Create a mask from the region
        galaxy_mask = self.galaxy_region.to_mask(frame.xsize, frame.ysize)

    # -----------------------------------------------------------------

    def remove_stars(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def remove_saturation(self):

        """
        This function ...
        :return:
        """

        # Check where the galaxy mask overlaps with the segmentation map
        overlapping, not_overlapping, overlapping_segments, not_overlapping_segments = masks.split_overlap(mask, galaxy_mask, return_segments=True)

        # Set mask for sources outside of the galaxy contour
        mask = not_overlapping

        # Inform the user
        log.info("Interpolating over the galaxy ...")

        # Find contours
        sigma_level = 4.0
        from pts.magic.analysis import sources
        contours = sources.find_contours(overlapping_segments, overlapping_segments, sigma_level)

        # Construct sources
        for contour in contours:

            # Create a source from the aperture
            source = Source.from_ellipse(frame, contour, 1.3)

            y_min = source.cutout.y_min
            y_max = source.cutout.y_max
            x_min = source.cutout.x_min
            x_max = source.cutout.x_max

            # Set the source mask
            overlapping_segments_box = overlapping_segments[y_min:y_max, x_min:x_max]
            center_label = overlapping_segments_box[int(round(0.5*overlapping_segments_box.shape[0])), int(round(0.5*overlapping_segments_box.shape[1]))]
            source_mask = Mask(overlapping_segments_box == center_label)
            #source.mask = overlapping[y_min:y_max, x_min:x_max] simple way, but: may contain the masks of nearby sources (nearby shapes in the region)
            source.mask = source_mask

            source.estimate_background("polynomial", True)

            # Replace the frame with the estimated background
            source.background.replace(frame, where=source.mask)

    # -----------------------------------------------------------------

    

# -----------------------------------------------------------------
