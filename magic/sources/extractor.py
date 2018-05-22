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

# Import standard modules
import math
import numpy as np

# Import the relevant PTS classes and modules
from ..basics.mask import Mask
from ..region.ellipse import PixelEllipseRegion
from ..basics.coordinate import PixelCoordinate
from ..core.detection import Detection
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from ..tools import masks
from ...core.basics.animation import Animation
from ..region.list import PixelRegionList
from ..core.image import Image
from ..core.frame import Frame
from ...core.tools import filesystem as fs
from ..region.list import load_as_pixel_region_list
from ..region.point import PixelPointRegion
from pts.core.tools.utils import lazyproperty
from ..core.mask import Mask as newMask

# -----------------------------------------------------------------

class SourceExtractor(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SourceExtractor, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The image frame
        self.frame = None

        # The original minimum and maximum value
        self.minimum_value = None
        self.maximum_value = None

        # The mask of nans
        self.nan_mask = None

        # Regions
        self.galaxy_region = None
        self.star_region = None
        self.saturation_region = None
        self.other_region = None

        # The animation
        self.animation = None

        # Special mask
        self.special_mask = None

        # Segmentation maps
        self.galaxy_segments = None
        self.star_segments = None
        self.other_segments = None

        # The total mask of removed sources
        self.mask = None

        # The list of sources
        self.sources = []
        self.labels = []

        # STATISTICS
        self.ngalaxy_sources = 0
        self.nstar_sources = 0
        self.nother_sources = 0
        self.nforeground = 0
        self.nfailed = 0
        self.nsuccess = 0
        self.nwith_saturation = 0

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Load the sources
        self.load_sources()

        # 3. Create the mask of all sources to be removed
        self.create_mask()

        # 4. For each source, check the pixels in the background that belong to an other source
        self.set_cross_contamination()

        # 5. Remove the sources
        self.remove_sources()

        # 6. Fix extreme values that showed up during the interpolation steps
        self.fix_extreme_values()

        # 7. Set nans back into the frame
        self.set_nans()

        # 8. Writing
        if self.config.output is not None and self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SourceExtractor, self).setup(**kwargs)

        # Set the image frame
        if "frame" in kwargs: self.frame = kwargs.pop("frame")
        else: self.load_frame()

        # Load the region lists
        self.load_regions(**kwargs)

        # Load the segmentation maps
        self.load_segments(**kwargs)

        # Initialize the mask
        self.mask = Mask.empty_like(self.frame)

        # Remember the minimum and maximum value
        self.minimum_value = np.nanmin(self.frame)
        self.maximum_value = np.nanmax(self.frame)

        # Create a mask of the pixels that are NaNs
        self.nan_mask = Mask.is_nan(self.frame)
        self.frame[self.nan_mask] = 0.0

        # Make a reference to the animation
        self.animation = kwargs.pop("animation", None)

        # Create mask from special region
        if "special_region" in kwargs:
            special_region = kwargs.pop("special_region")
            self.special_mask = Mask.from_region(special_region, self.frame.xsize, self.frame.ysize) if special_region is not None else None

        # If making animation is enabled
        if self.config.animation:
            self.animation = Animation()
            self.animation.fps = 1

    # -----------------------------------------------------------------

    def load_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the image frame ...")

        # load
        self.frame = Frame.from_file(self.config.image)

    # -----------------------------------------------------------------

    @property
    def frame_name(self):

        """
        This function ...
        :return:
        """

        return self.frame.name

    # -----------------------------------------------------------------

    @property
    def filter(self):

        """
        This function ...
        :return:
        """

        return self.frame.filter

    # -----------------------------------------------------------------

    @property
    def filter_name(self):

        """
        This function ...
        :return:
        """

        return self.frame.filter_name

    # -----------------------------------------------------------------

    def load_regions(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Loading the regions ...")

        # Load the galaxy region
        if "galaxy_region" in kwargs: self.galaxy_region = kwargs.pop("galaxy_region")
        else:
            if "name" in kwargs: galaxy_region_path = self.input_path_file("galaxies_" + kwargs["name"] + ".reg")
            else:

                galaxy_region_path = self.input_path_file("galaxies.reg")

                # Check
                if not fs.is_file(galaxy_region_path):

                    if self.filter_name is not None:

                        galaxy_region_path = self.input_path_file("galaxies_" + self.filter_name + ".reg")
                        if not fs.is_file(galaxy_region_path):
                            log.warning("No galaxy regions file could be found")
                            galaxy_region_path = None

                    elif self.frame_name is not None:

                        galaxy_region_path = self.input_path_file("galaxies_" + self.frame_name + ".reg")
                        if not fs.is_file(galaxy_region_path):
                            log.warning("No galaxy regions file could be found")
                            galaxy_region_path = None

                    else:
                        log.warning("No galaxy regions file could be found")
                        galaxy_region_path = None

            # Load the galaxy regions
            if galaxy_region_path is not None: self.galaxy_region = load_as_pixel_region_list(galaxy_region_path, self.frame.wcs)

        # Load the star region
        if "star_region" in kwargs: self.star_region = kwargs.pop("star_region")
        else:
            if "name" in kwargs: star_region_path = self.input_path_file("stars_" + kwargs["name"] + ".reg")
            else:

                star_region_path = self.input_path_file("stars.reg")

                # Check
                if not fs.is_file(star_region_path):

                    if self.filter_name is not None:

                        star_region_path = self.input_path_file("stars_" + self.filter_name + ".reg")
                        if not fs.is_file(star_region_path):
                            log.warning("No star regions file could be found")
                            star_region_path = None

                    elif self.frame_name is not None:

                        star_region_path = self.input_path_file("stars_" + self.frame_name + ".reg")
                        if not fs.is_file(star_region_path):
                            log.warning("No star regions file could be found")
                            star_region_path = None

                    else:
                        log.warning("No star regions file could be found")
                        star_region_path = None

            # Load the star regions
            if star_region_path is not None: self.star_region = load_as_pixel_region_list(star_region_path, self.frame.wcs)

        # Load the saturation region
        if "saturation_region" in kwargs: self.saturation_region = kwargs.pop("saturation_region")
        else:
            if "name" in kwargs: saturation_region_path = self.input_path_file("saturation_" + kwargs["name"] + ".reg")
            else:

                saturation_region_path = self.input_path_file("saturation.reg")

                # Check
                if not fs.is_file(saturation_region_path):

                    if self.filter_name is not None:

                        saturation_region_path = self.input_path_file("saturation_" + self.filter_name + ".reg")
                        if not fs.is_file(saturation_region_path):
                            log.warning("No saturation regions file could be found")
                            saturation_region_path = None

                    elif self.frame_name is not None:

                        saturation_region_path = self.input_path_file("saturation_" + self.frame_name + ".reg")
                        if not fs.is_file(saturation_region_path):
                            log.warning("No saturation regions file could be found")
                            saturation_region_path = None
                    else:
                        log.warning("No saturation regions file could be found")
                        saturation_region_path = None

            # Load the saturation regions
            if saturation_region_path is not None: self.saturation_region = load_as_pixel_region_list(saturation_region_path, self.frame.wcs)

        # Load the region of other sources
        if "other_region" in kwargs: self.other_region = kwargs.pop("other_region")
        else:
            if "name" in kwargs: other_region_path = self.input_path_file("other_sources_" + kwargs["name"] + ".reg")
            else:

                other_region_path = self.input_path_file("other_sources.reg")

                # Check
                if not fs.is_file(other_region_path):

                    if self.filter_name is not None:

                        other_region_path = self.input_path_file("other_sources_" + self.filter_name + ".reg")
                        if not fs.is_file(other_region_path):
                            log.warning("No other regions file could be found")
                            other_region_path = None

                    elif self.frame_name is not None:

                        other_region_path = self.input_path_file("other_sources_" + self.frame_name + ".reg")
                        if not fs.is_file(other_region_path):
                            log.warning("No other regions file could be found")
                            other_region_path = None
                    else:
                        log.warning("No other regions file could be found")
                        other_region_path = None

            # Load the other regions
            if other_region_path is not None: self.other_region = load_as_pixel_region_list(other_region_path, self.frame.wcs)

        # Debugging
        if self.galaxy_region is not None: log.debug("Galaxy regions: PRESENT")
        else: log.debug("Galaxy regions: NOT PRESENT")
        if self.star_region is not None: log.debug("Star regions: PRESENT")
        else: log.debug("Star regions: NOT PRESENT")
        if self.saturation_region is not None: log.debug("Saturation regions: PRESENT")
        else: log.debug("Saturation regions: NOT PRESENT")
        if self.other_region is not None: log.debug("Other regions: PRESENT")
        else: log.debug("Other regions: NOT PRESENT")

    # -----------------------------------------------------------------

    def load_segments(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Loading the segmentation maps ...")

        # Load the image with segmentation maps
        segments = None
        if "segments" in kwargs: segments = kwargs.pop("segments")
        else:
            if "name" in kwargs: segments_path = self.input_path_file("segments_" + kwargs["name"] + ".fits")
            else:

                segments_path = self.input_path_file("segments.fits")

                if not fs.is_file(segments_path):

                    if self.filter_name is not None:

                        segments_path = self.input_path_file("segments_" + self.filter_name + ".fits")
                        if not fs.is_file(segments_path):
                            log.warning("No segmentation maps found, will be using regions to define the to be extracted patches")
                            segments_path = None

                    elif self.frame_name is not None:

                        segments_path = self.input_path_file("segments_" + self.frame_name + ".fits")
                        if not fs.is_file(segments_path):
                            log.warning("No segmentation maps found, will be using regions to define the to be extracted patches")
                            segments_path = None

                    else:
                        log.warning("No segmentation maps found, will be using regions to define the to be extracted patches")
                        segments_path = None

            # Load the segments
            if segments_path is not None: segments = Image.from_file(segments_path, no_filter=True)

        # If segments is not None
        if segments is not None:

            #print(segments.frames.keys())

            # Get the segmentation maps
            self.galaxy_segments = segments.frames["extended"] if "extended" in segments.frames else None
            self.star_segments = segments.frames["point"] if "point" in segments.frames else None
            self.other_segments = segments.frames["other"] if "other" in segments.frames else None

        # Debugging
        if self.galaxy_segments is not None: log.debug("Galaxy segments: PRESENT")
        else: log.debug("Galaxy segments: NOT PRESENT")
        if self.star_segments is not None: log.debug("Star segments: PRESENT")
        else: log.debug("Star segments: NOT PRESENT")
        if self.other_segments is not None: log.debug("Other segments: PRESENT")
        else: log.debug("Other segments: NOT PRESENT")

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sources ...")

        # Load the galaxy sources
        self.load_galaxy_sources()

        # Load the star sources
        if self.star_region is not None: self.load_star_sources()

        # Load the other sources
        if self.other_region is not None: self.load_other_sources()

    # -----------------------------------------------------------------

    def load_galaxy_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy sources ...")

        # Loop over the shapes in the galaxy region
        for shape in self.galaxy_region:

            #print(shape.label)
            #print(shape.meta)

            # Shapes without text are in this case just coordinates
            #if "text" not in shape.meta: continue
            if shape.label is None: continue

            # Debugging
            log.debug("Adding galaxy '" + shape.label + "' ...")

            # Segments are passed
            if self.galaxy_segments is not None:

                # Get the coordinate of the center for this galaxy
                center = shape.center

                #print("here")

                # Check the label of the corresponding segment
                label = self.galaxy_segments[int(center.y), int(center.x)]

                if label == 3 or (label == 2 and self.config.remove_companions):

                    # Create a source
                    source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

                    # Check whether it is a 'special' source
                    source.special = self.special_mask.masks(center) if self.special_mask is not None else False

                    self.ngalaxy_sources += 1

                    # Add the source to the list
                    self.sources.append(source)
                    self.labels.append(label)

            elif "principal" not in shape.label:

                # Create a source
                source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

                # Check whether it is a special source
                source.special = self.special_mask.masks(shape.center) if self.special_mask is not None else False

                self.ngalaxy_sources += 1

                # Add the source to the list
                self.sources.append(source)
                self.labels.append(shape.label)

    # -----------------------------------------------------------------

    def load_star_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the star sources ...")

        # Loop over all stars in the region
        for shape in self.star_region:

            #print(shape.label)

            # Ignore shapes without text, these should be just the positions of the peaks
            #if "text" not in shape.meta: continue
            if shape.label is None: continue

            # Ignore shapes with color red (stars without source)
            if shape.appearance["color"] == "red": continue

            # Get the star index
            index = int(shape.label)

            # Debugging
            log.debug("Adding star " + str(index) + " ...")

            # Get the saturation source
            saturation_source = self.find_saturation_source(index)

            # Check whether the star is a 'special' region
            special = self.special_mask.masks(shape.center) if self.special_mask is not None else False

            # Saturation source was found
            if saturation_source is not None:

                self.nwith_saturation += 1

                ## DILATION
                if self.config.dilate_saturation: self.dilate_saturation_source(saturation_source)

                # Set the source to be the saturation source
                source = saturation_source

            # Create a new source from the shape
            else: source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

            # Set special flag
            source.special = special

            # Increment
            self.nstar_sources += 1

            # Add it to the list
            self.sources.append(source)

            #
            self.labels.append(index)

    # -----------------------------------------------------------------

    def find_saturation_source(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        # Deubgging
        log.debug("Finding a saturation source for star " + str(index) + " ...")

        # Look whether a saturation source is present
        saturation_source = None

        # Check whether the star is a foreground star
        #if self.principal_mask.masks(shape.center): foreground = True

        # If there is a saturation region
        if self.saturation_region is not None:

            # Add the saturation sources
            # Loop over the shapes in the saturation region
            for j in range(len(self.saturation_region)):

                saturation_shape = self.saturation_region[j]

                #if "text" not in saturation_shape.meta: continue
                if saturation_shape.label is None: continue

                saturation_index = int(saturation_shape.label)

                if index != saturation_index: continue
                else:

                    # Remove the saturation shape from the region
                    saturation_shape = self.saturation_region.pop(j)

                    # Create saturation source
                    saturation_source = Detection.from_shape(self.frame, saturation_shape, self.config.source_outer_factor)

                    # Replace the saturation mask
                    segments_cutout = self.star_segments[saturation_source.y_slice, saturation_source.x_slice]
                    saturation_mask = Mask(segments_cutout == index)
                    saturation_source.mask = saturation_mask.fill_holes()

                    # Break the loop
                    break

        # Return the saturation source
        return saturation_source

    # -----------------------------------------------------------------

    def dilate_saturation_source(self, saturation_source):

        """
        This function ...
        :param saturation_source:
        :return: 
        """

        # factor = saturation_dilation_factor
        dilation_factor = self.config.saturation_dilation_factor

        saturation_source = saturation_source.zoom_out(dilation_factor, self.frame, keep_original_mask=True)

        mask_area = np.sum(saturation_source.mask)
        area_dilation_factor = dilation_factor ** 2.
        new_area = mask_area * area_dilation_factor

        ## Circular mask approximation

        # ellipse = find_contour(source.mask.astype(float), source.mask)
        # radius = ellipse.radius.norm

        mask_radius = math.sqrt(mask_area / math.pi)
        new_radius = math.sqrt(new_area / math.pi)

        kernel_radius = new_radius - mask_radius

        # Replace mask
        saturation_source.mask = saturation_source.mask.disk_dilation(radius=kernel_radius)

    # -----------------------------------------------------------------

    def load_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the other sources ...")

        # Loop over the shapes in the other sources region
        for shape in self.other_region:

            # This is a source found by SourceFinder
            if shape.label is not None:

                # Debugging
                log.debug("Adding other source '" + shape.label + "' ...")

                # Get integer label
                label = int(shape.label)

                # Create a source
                source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

                # Replace the source mask
                segments_cutout = self.other_segments[source.y_slice, source.x_slice]
                source.mask = Mask(segments_cutout == label).fill_holes()

                ## DILATION

                if self.config.dilate_other:

                    # DILATE SOURCE
                    # factor = other_dilation_factor

                    dilation_factor = self.config.other_dilation_factor

                    ## CODE FOR DILATION (FROM SOURCES MODULE)

                    source = source.zoom_out(dilation_factor, self.frame, keep_original_mask=True)

                    mask_area = np.sum(source.mask)
                    area_dilation_factor = dilation_factor ** 2.
                    new_area = mask_area * area_dilation_factor

                    ## Circular mask approximation

                    # ellipse = find_contour(source.mask.astype(float), source.mask)
                    # radius = ellipse.radius.norm

                    mask_radius = math.sqrt(mask_area / math.pi)
                    new_radius = math.sqrt(new_area / math.pi)

                    kernel_radius = new_radius - mask_radius

                    # Replace mask
                    source.mask = source.mask.disk_dilation(radius=kernel_radius)

                ## END DILATION CODE

            # This is a shape drawn by the user and added to the other sources region
            # # Create a source
            else: source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

            # Check whether source is 'special'
            source.special = self.special_mask.masks(shape.center) if self.special_mask is not None else False

            # Increment
            self.nother_sources += 1

            # Add the source to the list
            self.sources.append(source)
            self.labels.append(label)

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the mask of all sources to be removed ...")

        # Loop over all sources
        #for source in self.sources:
        index = 0
        while index < len(self.sources):

            # Get the current source
            source = self.sources[index]

            # If these pixels are already masked by an overlapping source (e.g. saturation), remove this source,
            # otherwise the area will be messed up
            current_mask_cutout = self.mask[source.y_slice, source.x_slice]
            if current_mask_cutout.covers(source.mask):
                self.sources.pop(index)
                continue

            # Adapt the mask
            self.mask[source.y_slice, source.x_slice] += source.mask

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def set_cross_contamination(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("For each source, checking which pixels in the neighborhood are contaminated by other sources ...")

        # Loop over all sources
        for source in self.sources:

            # Create the contamination mask for this source
            other_sources_mask = Mask.empty_like(source.cutout)
            other_sources_mask[source.background_mask] = self.mask[source.y_slice, source.x_slice][source.background_mask]
            source.contamination = other_sources_mask

    # -----------------------------------------------------------------

    @property
    def nsources(self):

        """
        This function ...
        :return: 
        """

        return len(self.sources)

    # -----------------------------------------------------------------

    def remove_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the frame over the masked pixels ...")

        nsources = len(self.sources)
        count = 0

        # Set principal ellipse for the source extraction animation
        if self.animation is not None: self.animation.principal_shape = self.principal_shape

        # Loop over all sources and remove them from the frame
        for label, source in zip(self.labels, self.sources):

            # Debugging
            log.debug("Estimating background and replacing the frame pixels of source " + str(count+1) + " of " + str(nsources) + " ...")

            # Check whether the source is in front of the principal galaxy
            #foreground = self.principal_mask.masks(source.center)
            if self.principal_mask is not None: foreground = masks.overlap(self.principal_mask[source.y_slice, source.x_slice], source.mask)
            else: foreground = False

            if foreground: self.nforeground += 1

            # SKip foreground if requested
            if self.config.only_foreground and not foreground: continue

            # Disable sigma-clipping for estimating background when the source is foreground to the principal galaxy (to avoid clipping the galaxy's gradient)
            sigma_clip = self.config.sigma_clip if not foreground else False

            # Debugging
            log.debug("Sigma-clipping enabled for estimating background gradient for this source" if sigma_clip else "Sigma-clipping disabled for estimating background gradient for this source")

            # If these pixels are already replaced by an overlapping source (e.g. saturation), skip this source,
            # otherwise the area will be messed up
            #current_mask_cutout = self.mask[source.y_slice, source.x_slice]
            #if current_mask_cutout.covers(source.mask):
            #    count += 1
            #    continue
            ## ==> this is now also done in create_mask

            # Estimate the background
            try:
                source.estimate_background(self.config.interpolation_method, sigma_clip=sigma_clip)
            except ValueError: # ValueError: zero-size array to reduction operation minimum which has no identity
                # in: limits = (np.min(known_points), np.max(known_points)) [inpaint_biharmonic]
                self.nfailed += 1
                count += 1
                continue

            # Adapt the mask
            #self.mask[source.y_slice, source.x_slice] += source.mask # this is now done beforehand, in the create_mask function

            # Add frame to the animation
            if self.animation is not None and (self.principal_mask is None or self.principal_mask.masks(source.center)) and self.animation.nframes <= 20:
                self.animation.add_source(source)

            # Replace the pixels by the background
            source.background.replace(self.frame, where=source.mask)

            # Increment
            self.nsuccess += 1

            #if not sigma_clip:
            #    # source.plot()

            #    plotting.plot_removal(source.cutout, source.mask, source.background,
            #                          self.frame[source.y_slice, source.x_slice])

            count += 1

    # -----------------------------------------------------------------

    def fix_extreme_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fixing extreme values that were introduced during the interpolation steps ...")

        self.frame[self.frame < self.minimum_value] = self.minimum_value
        self.frame[self.frame > self.maximum_value] = self.maximum_value

    # -----------------------------------------------------------------

    def set_nans(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting original NaN-pixels back to NaN ...")

        # Set the NaN pixels to zero in the frame
        self.frame[self.nan_mask] = float("nan")

    # -----------------------------------------------------------------

    def write(self):

        """
        THis function ...
        :return:
        """

        # Inform the suer
        log.info("Writing ...")

        # Write the animation
        if self.animation is not None: self.write_animation()

        # Write the resulting frame
        self.write_frame()

        # Write the mask
        self.write_mask()

    # -----------------------------------------------------------------

    def write_animation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the animation ...")

        # Save the animation
        path = self.output_path_file("animation.gif")
        self.animation.saveto(path)

    # -----------------------------------------------------------------

    def write_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the result ...")

        # Determine the path to the resulting FITS file
        path = self.output_path_file("extracted.fits")

        # Save the resulting image as a FITS file
        self.frame.saveto(path)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mask ...")

        # Determine the path to the mask
        path = self.output_path_file("mask.fits")

        # Save the total mask as a FITS file
        Frame(self.mask.astype(float)).saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def principal_shape(self):

        """
        This function ...
        :return:
        """

        if self.galaxy_region is None: return None

        largest_shape = None

        # Loop over all the shapes in the galaxy region
        for shape in self.galaxy_region:

            # Skip single coordinates
            if isinstance(shape, PixelCoordinate): continue

            if shape.label is not None and "principal" in shape.label: return shape
            if "text" in shape.meta and "principal" in shape.meta["text"]: return shape

            if not isinstance(shape, PixelEllipseRegion) and not isinstance(shape, PixelPointRegion): return shape

            semimajor_axis_length = shape.semimajor
            if largest_shape is None or semimajor_axis_length > largest_shape.semimajor: largest_shape = shape

        # Return the largest shape
        return largest_shape

    # -----------------------------------------------------------------

    @lazyproperty
    def principal_mask(self):

        """
        This function ...
        :return:
        """

        if self.galaxy_segments is not None: return newMask.where(self.galaxy_segments, 1)
        elif self.principal_shape is not None: return self.principal_shape.to_mask(self.frame.xsize, self.frame.ysize)
        else: return None

# -----------------------------------------------------------------
