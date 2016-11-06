#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.batchfinder Contains the BatchSourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.convolution import Gaussian2DKernel

# Import the relevant PTS classes and modules
from .galaxyfinder import GalaxyFinder
from .starfinder import StarFinder
from .trainedfinder import TrainedFinder
from ..basics.mask import Mask
from ..catalog.builder import CatalogBuilder
from ..catalog.synchronizer import CatalogSynchronizer
from ..tools import wavelengths
from ...core.tools import tables
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ..core.dataset import DataSet
from ..catalog.importer import CatalogImporter
from ...core.tools import filesystem as fs
from ..basics.skyregion import SkyRegion
from ..core.image import Image
from ..core.frame import Frame
from ..tools import statistics

# -----------------------------------------------------------------

class BatchSourceFinder(Configurable):

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
        super(BatchSourceFinder, self).__init__(config)

        # -- Attributes --

        # The frames
        self.frames = dict()

        # The masks
        self.special_masks = dict()
        self.ignore_masks = dict()

        # Downsampled images
        self.downsampled = None
        self.original_wcs = None

        # The galactic and stellar catalog
        self.galactic_catalog = None
        self.stellar_catalog = None

        # The regions covering areas that should be ignored throughout the entire extraction procedure
        self.special_region = None
        self.ignore_region = None

        # The name of the principal galaxy
        self.galaxy_name = None

        # The regions
        self.galaxy_regions = dict()
        self.star_regions = dict()
        self.saturation_regions = dict()
        self.other_regions = dict()

        # The segmentation maps
        self.segments = dict()

        # The finders
        self.star_finder = None
        self.galaxy_finder = None
        self.trained_finder = None

        # Galaxy and star lists
        self.galaxies = dict()
        self.stars = dict()

        # The PSFs
        self.psfs = dict()

        # The statistics
        self.statistics = dict()

    # -----------------------------------------------------------------

    def add_frame(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.frames:

            wcs = self.frames[name].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegion()

        # Add the bounding boxes as sky rectangles
        for name in self.frames: boxes_region.append(self.frames[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get catalogs
        self.get_catalogs()

        # 3. Find the galaxies
        self.find_galaxies()
        
        # 3. Find the stars
        if self.config.find_stars: self.find_stars()

        # 4. Look for other sources
        if self.config.find_other_sources: self.find_other_sources()

        # 5. Build and update catalog
        #self.build_and_synchronize_catalog()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BatchSourceFinder, self).setup(**kwargs)

        # Load the frames
        self.load_frames()

        # Load special region
        self.special_region = SkyRegion.from_file(self.config.special_region) if self.config.special_region is not None else None

        # Load ignore region
        self.ignore_region = SkyRegion.from_file(self.config.ignore_region) if self.config.ignore_region is not None else None

        # Create the masks
        self.create_masks()

        # Create the finders
        self.galaxy_finder = GalaxyFinder(self.config.galaxies)
        self.star_finder = StarFinder(self.config.stars)
        self.trained_finder = TrainedFinder(self.config.other_sources)

        # DOWNSAMPLE ??

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the frame(s) ...")

        # Create new dataset
        if self.config.dataset.endswith(".fits"):

            # Load the frame
            frame = Frame.from_file(self.config.dataset)

            # Determine the name for this image
            name = str(frame.filter)

            # Add the frame
            self.add_frame(name, frame)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the frames
            self.frames = dataset.get_frames()

        # Invalid value for 'dataset'
        else: raise ValueError("Parameter 'dataset' must be filename of a dataset file (.dat) or a FITS file (.fits)")

    # -----------------------------------------------------------------

    def create_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the masks ...")

        if self.special_region is not None:

            for name in self.frames:

                # Create the mask
                special_mask = Mask.from_region(self.special_region, self.frames[name].xsize, self.frames[name].ysize)

                self.special_masks[name] = special_mask

        if self.ignore_region is not None:

            for name in self.frames:

                # Create the mask
                ignore_mask = Mask.from_region(self.ignore_region, self.frames[name].xsize, self.frames[name].ysize)

                self.ignore_masks[name] = ignore_mask

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting galactic and stellar catalogs ...")

        # Create a CatalogImporter instance
        catalog_importer = CatalogImporter()

        # Set configuration
        if self.config.galactic_catalog_file is not None:
            catalog_importer.config.galaxies.use_catalog_file = True
            catalog_importer.config.galaxies.catalog_path = self.config.galactic_catalog_file

        if self.config.stellar_catalog_file is not None:
            catalog_importer.config.stars.use_catalog_file = True
            catalog_importer.config.stars.catalog_path = self.config.stellar_catalog_file

        # Get the coordinate box and minimum pixelscale
        coordinate_box = self.bounding_box
        min_pixelscale = self.min_pixelscale

        # Run the catalog importer
        catalog_importer.run(coordinate_box=coordinate_box, pixelscale=min_pixelscale)  # work with coordinate box instead ? image.coordinate_box ?

        # Set the catalogs
        self.galactic_catalog = catalog_importer.galactic_catalog
        self.stellar_catalog = catalog_importer.stellar_catalog

        galactic_catalog_path = self.output_path_file("galaxies.cat")
        stellar_catalog_path = self.output_path_file("stars.cat")

        tables.write(self.galactic_catalog, galactic_catalog_path)
        tables.write(self.stellar_catalog, stellar_catalog_path)

    # -----------------------------------------------------------------

    def find_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the galaxies ...")

        # Loop over the images
        for name in self.frames:

            # Get the frame
            frame = self.frames[name]

            special_mask = self.special_masks[name] if name in self.special_masks else None
            ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
            bad_mask = None

            # Run the galaxy finder
            self.galaxy_finder.run(frame=frame, catalog=self.galactic_catalog, special_mask=special_mask, ignore_mask=ignore_mask, bad_mask=bad_mask)

            # Set the name of the principal galaxy
            #self.galaxy_name = self.galaxy_finder.principal.name

            # Get the galaxy region
            #galaxy_sky_region = self.galaxy_finder.finder.galaxy_sky_region
            #if galaxy_sky_region is not None:
            #    galaxy_region = galaxy_sky_region.to_pixel(image.wcs)

            if self.galaxy_finder.region is not None:

                galaxy_sky_region = self.galaxy_finder.region.to_sky(frame.wcs)

                #if self.downsampled:
                #    sky_region = self.galaxy_sky_region
                #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
                #ele: return self.galaxy_finder.region

                self.galaxy_regions[name] = galaxy_sky_region

            if self.galaxy_finder.segments is not None:

                # Create an image with the segmentation maps
                self.segments[name] = Image("segments")

                #if self.galaxy_finder.segments is None: return None
                #if self.downsampled:

                #    segments = self.galaxy_finder.segments
                #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
                #    upsampled.unpad(self.pad_x, self.pad_y)
                #    return upsampled

                #else: return self.galaxy_finder.segments

                galaxy_segments = self.galaxy_finder.segments

                # Add the segmentation map of the galaxies
                self.segments[name].add_frame(galaxy_segments, "galaxies")

            # Set galaxies
            self.galaxies[name] = self.galaxy_finder.galaxies

            # Inform the user
            log.success("Finished finding the galaxies for '" + name + "' ...")

            # Clear the galaxy finder
            self.galaxy_finder.clear()

    # -----------------------------------------------------------------
    
    def find_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the stars ...")

        # Loop over the images
        for name in self.frames:

            # Get the frame
            frame = self.frames[name]

            # Run the star finder if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
            if frame.wavelength is None or frame.wavelength < wavelengths.ranges.ir.mir.max:

                # Inform the user
                log.info("Finding the stars ...")

                special_mask = self.special_masks[name] if name in self.special_masks else None
                ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
                bad_mask = None

                # Run the star finder
                self.star_finder.run(frame=frame, galaxies=self.galaxies[name], catalog=self.stellar_catalog, special_mask=special_mask, ignore_mask=ignore_mask, bad_mask=bad_mask)

                if self.star_finder.star_region is not None:

                    star_sky_region = self.star_finder.star_region.to_sky(frame.wcs)

                    #if self.downsampled:
                    #    sky_region = self.star_sky_region
                    #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
                    #else: return self.star_finder.star_region

                    self.star_regions[name] = star_sky_region

                if self.star_finder.saturation_region is not None:

                    saturation_sky_region = self.star_finder.saturation_region.to_sky(frame.wcs)

                    #if self.downsampled:
                    #    sky_region = self.saturation_sky_region
                    #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
                    #else: return self.star_finder.saturation_region

                    self.saturation_regions[name] = saturation_sky_region

                if self.star_finder.segments is not None:

                    #if self.star_finder.segments is None: return None
                    #if self.downsampled:
                    #    segments = self.star_finder.segments
                    #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
                    #    upsampled.unpad(self.pad_x, self.pad_y)
                    #    return upsampled
                    #else: return self.star_finder.segments

                    star_segments = self.star_finder.segments

                    # Add the segmentation map of the saturated stars
                    self.segments[name].add_frame(star_segments, "stars")

                # Set the stars
                self.stars[name] = self.star_finder.stars

                # kernel = self.star_finder.kernel # doesn't work when there was no star extraction on the image, self.star_finder does not have attribute image thus cannot give image.fwhm
                # Set the kernel (PSF)
                if self.star_finder.config.use_frame_fwhm and frame.fwhm is not None:

                    fwhm = frame.fwhm.to("arcsec").value / frame.average_pixelscale.to("arcsec/pix").value
                    sigma = fwhm * statistics.fwhm_to_sigma
                    kernel = Gaussian2DKernel(sigma)

                else: kernel = self.star_finder.kernel

                # Set the PSF
                self.psfs[name] = kernel

                # Get the statistics
                self.statistics[name] = self.star_finder.get_statistics()

                # Show the FWHM
                log.info("The FWHM that could be fitted to the point sources in the " + name + " image is " + str(self.statistics[name].fwhm))

                # Inform the user
                log.success("Finished finding the stars for '" + name + "' ...")

                # Clear the star finder
                self.star_finder.clear()

            # No star subtraction for this image
            else: log.info("Finding stars will not be performed on this frame")

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources in the frame not in the catalog ...")

        # Loop over the frames
        for name in self.frames:

            # Get the frame
            frame = self.frames[name]

            # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
            if frame.wavelength is not None and frame.wavelength > wavelengths.ranges.ir.mir.max: self.trained_finder.config.classify = False
            else: self.trained_finder.config.classify = True

            special_mask = self.special_masks[name] if name in self.special_masks else None
            ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
            bad_mask = None

            # Run the trained finder just to find sources
            self.trained_finder.run(frame=frame, galaxies=self.galaxies[name], stars=self.stars[name], special_mask=special_mask,
                                    ignore_mask=ignore_mask, bad_mask=bad_mask, galaxy_segments=self.segments[name].frames.galaxies,
                                    star_segments=self.segments[name].frames.stars, kernel=self.psfs[name])

            if self.trained_finder.region is not None:

                other_sky_region = self.trained_finder.region.to_sky(frame.wcs)

                #if self.downsampled:
                #    sky_region = self.other_sky_region
                #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
                #else: return self.trained_finder.region

                # Add the region
                self.other_regions[name] = other_sky_region

            if self.trained_finder.segments is not None:

                #if self.trained_finder.segments is None: return None
                #if self.downsampled:
                #    segments = self.trained_finder.segments
                #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
                #    upsampled.unpad(self.pad_x, self.pad_y)
                #    return upsampled
                #else: return self.trained_finder.segments

                other_segments = self.trained_finder.segments

                # Add the segmentation map of the other sources
                self.segments[name].add_frame(other_segments, "other_sources")

            # Inform the user
            log.success("Finished finding other sources for '" + name + "' ...")

    # -----------------------------------------------------------------

    def build_and_synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the catalog
        if self.config.build_catalogs: self.build_catalog()

        # Synchronize the catalog
        if self.config.build_catalogs and self.config.synchronize_catalogs: self.synchronize_catalog()

    # -----------------------------------------------------------------

    def build_catalog(self):

        """
        This function ...
        :return:
        """

        # Build the stellar catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.frame.wavelength is None or self.frame.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Building the stellar catalog ...")

            # Run the catalog builder
            self.catalog_builder.run(self.frame, self.galaxy_finder, self.star_finder, self.trained_finder)

            # Inform the user
            log.success("Stellar catalog built")

    # -----------------------------------------------------------------

    def synchronize_catalog(self):

        """
        This function ...
        :return:
        """

        # Synchronize the catalog if the wavelength of this image is smaller than 25 micron (or the wavelength is unknown)
        if self.frame.wavelength is None or self.frame.wavelength < wavelengths.ranges.ir.mir.max:

            # Inform the user
            log.info("Synchronizing with the DustPedia catalog ...")

            # Run the catalog synchronizer
            self.catalog_synchronizer.run(self.frame.frames.primary, self.galaxy_name, self.catalog_builder.galactic_catalog, self.catalog_builder.stellar_catalog)

            # Inform the user
            log.success("Catalog synchronization done")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write regions
        self.write_regions()

        # Write segmentation maps
        self.write_segments()

        # 1. Write
        #self.write_galactic_catalogs()

        # 2. Write
        #self.write_stellar_catalogs()

        # 3. Write ...
        self.write_statistics()

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the regions ...")

        # Loop over the regions
        for name in self.galaxy_regions:

            #galaxy_region = galaxy_sky_region.to_pixel(image.wcs)

            # Determine the path
            path = self.output_path_file("galaxies_" + name + ".reg") if len(self.frames) > 1 else self.output_path_file("galaxies.reg")

            # Save
            self.galaxy_regions[name].to_pixel(self.frames[name].wcs).save(path)

        # Loop over the star regions
        for name in self.star_regions:

            #star_region = star_sky_region.to_pixel(image.wcs)

            # Determine the path
            path = self.output_path_file("stars_" + name + ".reg") if len(self.frames) > 1 else self.output_path_file("stars.reg")

            # Save
            self.star_regions[name].to_pixel(self.frames[name].wcs).save(path)

        # Loop over the saturation regions
        for name in self.saturation_regions:

            #saturation_region = saturation_sky_region.to_pixel(image.wcs)

            # Determine the path
            path = self.output_path_file("saturation_" + name + ".reg") if len(self.frames) > 1 else self.output_path_file("saturation.reg")

            # Save
            self.saturation_regions[name].to_pixel(self.frames[name].wcs).save(path)

        # Loop over the other regions
        for name in self.other_regions:

            #other_region = other_sky_region.to_pixel(image.wcs)

            # Determine the path
            path = self.output_path_file("other_sources_" + name + ".reg") if len(self.frames) > 1 else self.output_path_file("other_sources.reg")

            # Save
            self.other_regions[name].to_pixel(self.frames[name].wcs).save(path)

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the segmentation maps ...")

        for name in self.segments:

            # Save the FITS file with the segmentation maps
            path = self.output_path_file("segments_" + name + ".fits") if len(self.frames) > 1 else self.output_path_file("segments.fits")

            # Save
            self.segments[name].save(path)

    # -----------------------------------------------------------------

    def write_galactic_catalogs(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.galactic_catalog_path)
        if path is None:
            log.error("Galactic catalog path is not defined, skipping writing galactic catalog ...")
            return

        # Inform the user
        log.info("Writing galactic catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.catalog_builder.galactic_catalog, path)

    # -----------------------------------------------------------------

    def write_stellar_catalogs(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the catalog file
        path = self.full_output_path(self.config.writing.stellar_catalog_path)
        if path is None:
            log.error("Stellar catalog path is not defined, skipping writing stellar catalog ...")
            return

        # Inform the user
        log.info("Writing stellar catalog to " + path + " ...")

        # Write the catalog to file
        tables.write(self.catalog_builder.stellar_catalog, path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing statistics ...")

        # Loop over the image names
        for name in self.statistics:

            # Determine the path to the statistics file
            path = self.output_path_file("statistics" + name + ".dat") if len(self.frames) > 1 else self.output_path_file("statistics.dat")

            # Open the file, write the info
            with open(path, 'w') as statistics_file:
                statistics_file.write("FWHM: " + str(self.statistics[name].fwhm) + "\n")

# -----------------------------------------------------------------
