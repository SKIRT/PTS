#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.batchfinder Contains the SourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from multiprocessing import Pool

# Import astronomical modules
from astropy.convolution import Gaussian2DKernel

# Import the relevant PTS classes and modules
from .galaxyfinder import GalaxyFinder
from .starfinder import StarFinder
from .trainedfinder import TrainedFinder
from ..basics.mask import Mask
from ..tools import wavelengths
from ...core.tools import tables
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ..core.dataset import DataSet
from ..catalog.importer import CatalogImporter
from ...core.tools import filesystem as fs
from ..region.list import SkyRegionList
from ..core.image import Image
from ..core.frame import Frame
from ..tools import statistics

# -----------------------------------------------------------------

class SourceFinder(Configurable):

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
        super(SourceFinder, self).__init__(config)

        # -- Attributes --

        # The process pool
        self.pool = None

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

        # Ignore images
        self.ignore = []
        self.ignore_stars = []
        self.ignore_other_sources = []

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

        # Settings for the star finder for different bands
        self.star_finder_settings = dict()

        # Galaxy and star lists
        self.galaxies = dict()
        self.stars = dict()

        # The PSFs
        self.psfs = dict()

        # The statistics
        self.statistics = dict()

    # -----------------------------------------------------------------

    def add_frame(self, name, frame, star_finder_settings=None):

        """
        This function ...
        :param name:
        :param frame:
        :param star_finder_settings:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

        # Set the settings
        if star_finder_settings is not None: self.star_finder_settings[name] = star_finder_settings

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
        boxes_region = SkyRegionList()

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
        super(SourceFinder, self).setup(**kwargs)

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

        # Load the images (from config or input kwargs)
        if "frames" in kwargs: self.frames = kwargs.pop("frames")
        elif "dataset" in kwargs:
            dataset = kwargs.pop("dataset")
            self.frames = dataset.get_frames()
        else: self.load_frames()

        # Get the settings
        if "star_finder_settings" in kwargs: self.star_finder_settings = kwargs.pop("star_finder_settings")

        # Ignore certain images
        self.ignore = kwargs.pop("ignore", [])

        # Ignore stars in certain images
        self.ignore_stars = kwargs.pop("ignore_stars", [])

        # Ignore other sources in certain images
        self.ignore_other_sources = kwargs.pop("ignore_other_sources", [])

        # Load special region
        self.special_region = SkyRegionList.from_file(self.config.special_region) if self.config.special_region is not None else None

        # Load ignore region
        self.ignore_region = SkyRegionList.from_file(self.config.ignore_region) if self.config.ignore_region is not None else None

        # Create the masks
        self.create_masks()

        # Initialize images for the segmentation maps
        for name in self.frames: self.segments[name] = Image("segments")

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

        # Dictionary to keep result handles
        results = dict()

        # Loop over the images
        for name in self.frames:

            # Ignore if requested
            if name in self.ignore: continue

            # Get the frame
            frame = self.frames[name]

            # Get masks
            special_mask = self.special_masks[name] if name in self.special_masks else None
            ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
            bad_mask = None

            # Get configuration
            config = self.config.galaxies.copy()

            # Do the detection
            result = self.pool.apply_async(detect_galaxies, args=(frame, self.galactic_catalog, config, special_mask, ignore_mask, bad_mask,))
            results[name] = result

        # Process results
        for name in results:

            # Get result
            galaxies, region_list, segments = results[name].get()

            # Set galaxies
            self.galaxies[name] = galaxies

            # Set region list
            self.galaxy_regions[name] = region_list

            # Set segmentation map
            # Add the segmentation map of the galaxies
            self.segments[name].add_frame(segments, "galaxies")

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------
    
    def find_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the stars ...")

        # Dictionary to keep result handles
        results = dict()

        # Loop over the images
        for name in self.frames:

            # Ignore if requested
            if name in self.ignore: continue
            if name in self.ignore_stars: continue

            # Get the frame
            frame = self.frames[name]

            # Don't run the star finder if the wavelength of this image is greater than 25 micron
            if frame.wavelength is not None or frame.wavelength > wavelengths.ranges.ir.mir.max:

                # No star subtraction for this image
                log.info("Finding stars will not be performed on this frame")
                continue

            # Inform the user
            log.info("Finding the stars ...")

            # Get masks
            special_mask = self.special_masks[name] if name in self.special_masks else None
            ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
            bad_mask = None

            # Create configuration
            config = self.config.stars.copy()
            if name in self.star_finder_settings: config.set_items(self.star_finder_settings[name])

            # Do the detection
            result = self.pool.apply_async(detect_stars, args=(frame, self.galaxies[name], self.stellar_catalog, config, special_mask, ignore_mask, bad_mask,))
            results[name] = result

        # Process results
        for name in results:

            # Get result
            # stars, star_region_list, saturation_region_list, star_segments, kernel, statistics
            stars, star_region_list, saturation_region_list, star_segments, kernel, statistics = results[name].get()

            # Set stars
            self.stars[name] = stars

            # Set star region list
            self.star_regions[name] = star_region_list

            # Set saturation region list
            self.saturation_regions[name] = saturation_region_list

            # Set segmentation map
            # Add the segmentation map of the galaxies
            self.segments[name].add_frame(star_segments, "stars")

            # Set the PSF
            self.psfs[name] = kernel

            # Get the statistics
            self.statistics[name] = statistics

            # Show the FWHM
            log.info("The FWHM that could be fitted to the point sources in the " + name + " image is " + str(self.statistics[name].fwhm))

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources in the frame not in the catalog ...")

        # Dictionary to keep result handles
        results = dict()

        # Loop over the frames
        for name in self.frames:

            # Ignore if requested
            if name in self.ignore: continue
            if name in self.ignore_other_sources: continue

            # Get the frame
            frame = self.frames[name]

            # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
            if frame.wavelength is not None and frame.wavelength > wavelengths.ranges.ir.mir.max: self.trained_finder.config.classify = False
            else: self.trained_finder.config.classify = True

            # Get masks
            special_mask = self.special_masks[name] if name in self.special_masks else None
            ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
            bad_mask = None

            # Create the configuration
            config = self.config.other_sources.copy()

            # Get other input
            galaxies = self.galaxies[name]
            stars = self.stars[name]
            galaxy_segments = self.segments[name].frames.galaxies
            star_segments = self.segments[name].frames.stars
            kernel = self.psfs[name]

            # Do the detection
            # frame, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask, bad_mask
            result = self.pool.apply_async(detect_other, args=(frame, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask, bad_mask,))
            results[name] = result

        # Process results
        for name in results:

            # Get the result
            region_list, segments = results[name].get()

            # Add the region
            self.other_regions[name] = region_list

            # Add the segmentation map of the other sources
            self.segments[name].add_frame(segments, "other_sources")

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

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

def detect_galaxies(frame, catalog, config, special_mask, ignore_mask, bad_mask):

    """
    This function ...
    :param frame:
    :param catalog:
    :param config:
    :param special_mask:
    :param ignore_mask:
    :param bad_mask:
    :return:
    """

    # Create the galaxy finder
    galaxy_finder = GalaxyFinder(config)

    # Run the galaxy finder
    galaxy_finder.run(frame=frame, catalog=catalog, special_mask=special_mask, ignore_mask=ignore_mask, bad_mask=bad_mask)

    # Set the name of the principal galaxy
    # self.galaxy_name = self.galaxy_finder.principal.name

    # Get the galaxy region
    # galaxy_sky_region = self.galaxy_finder.finder.galaxy_sky_region
    # if galaxy_sky_region is not None:
    #    galaxy_region = galaxy_sky_region.to_pixel(image.wcs)

    if galaxy_finder.region is not None:

        galaxy_sky_region = galaxy_finder.region.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.galaxy_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # ele: return self.galaxy_finder.region

        # Get region list
        region_list = galaxy_sky_region

    if galaxy_finder.segments is not None:

        # if self.galaxy_finder.segments is None: return None
        # if self.downsampled:

        #    segments = self.galaxy_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled

        # else: return self.galaxy_finder.segments

        # Get the segments
        segments = galaxy_finder.segments

    # Get the galaxies
    galaxies = galaxy_finder.galaxies

    # Inform the user
    log.success("Finished finding the galaxies for '" + frame.name + "' ...")

    # Return the output
    return galaxies, region_list, segments

# -----------------------------------------------------------------

def detect_stars(frame, galaxies, catalog, config, special_mask, ignore_mask, bad_mask):

    """
    This function ...
    :param frame:
    :param galaxies:
    :param catalog:
    :param config:
    :param special_mask:
    :param ignore_mask:
    :param bad_mask:
    :return:
    """

    # Create the star finder
    finder = StarFinder(config)

    # Run the star finder
    finder.run(frame=frame, galaxies=galaxies, catalog=catalog, special_mask=special_mask, ignore_mask=ignore_mask, bad_mask=bad_mask)

    if finder.star_region is not None:

        star_sky_region = finder.star_region.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.star_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.star_finder.star_region

        star_region_list = star_sky_region

    if finder.saturation_region is not None:

        saturation_sky_region = finder.saturation_region.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.saturation_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.star_finder.saturation_region

        saturation_region_list = saturation_sky_region

    if finder.segments is not None:

        # if self.star_finder.segments is None: return None
        # if self.downsampled:
        #    segments = self.star_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled
        # else: return self.star_finder.segments

        star_segments = finder.segments

    # Set the stars
    stars = finder.stars

    # kernel = self.star_finder.kernel # doesn't work when there was no star extraction on the image, self.star_finder does not have attribute image thus cannot give image.fwhm
    # Set the kernel (PSF)
    if finder.config.use_frame_fwhm and frame.fwhm is not None:

        fwhm = frame.fwhm.to("arcsec").value / frame.average_pixelscale.to("arcsec/pix").value
        sigma = fwhm * statistics.fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma)

    else: kernel = finder.kernel

    # Inform the user
    log.success("Finished finding the stars for '" + frame.name + "' ...")

    # Return the output
    return stars, star_region_list, saturation_region_list, star_segments, kernel, statistics

# -----------------------------------------------------------------

def detect_other(frame, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask, bad_mask):

    """
    This function ...
    :param frame:
    :param config:
    :param galaxies:
    :param stars:
    :param galaxy_segments:
    :param star_segments:
    :param kernel:
    :param special_mask:
    :param ignore_mask:
    :param bad_mask:
    :return:
    """

    # Create the trained finder
    finder = TrainedFinder(config)

    # Run the trained finder just to find sources
    finder.run(frame=frame, galaxies=galaxies, stars=stars, special_mask=special_mask,
                            ignore_mask=ignore_mask, bad_mask=bad_mask,
                            galaxy_segments=galaxy_segments,
                            star_segments=star_segments, kernel=kernel)

    if finder.region is not None:

        other_sky_region = finder.region.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.other_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.trained_finder.region

    if finder.segments is not None:

        # if self.trained_finder.segments is None: return None
        # if self.downsampled:
        #    segments = self.trained_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled
        # else: return self.trained_finder.segments

        other_segments = finder.segments

    # Inform the user
    log.success("Finished finding other sources for '" + frame.name + "' ...")

    # Return the output
    return other_sky_region, other_segments

# -----------------------------------------------------------------

# From: https://github.com/fred3m/astropyp/blob/master/astropyp/phot/calibrate.py

import numpy as np
from astropy.table import Table, join

def clip_color_outliers(color, ref_color, dbscan_kwargs={}, show_plots=True):

    """
    Use DBSCAN clustering to only select the points in color-color space
    that are part of the main cluster

    Parameters
    ----------
    color: array-like
        Array of instrumental colors
    ref_color: array-like
        Array of colors from a reference catalog
    dbscan_kwargs: dict
        Dictionary of keyword arguments to pass to DBSCAN. Typical
        parameters are 'eps', the maximum separation for points in a
        group and 'min_samples', the minimum number of points for a
        cluster to be grouped together. No kwargs are required.

    Returns
    -------
    idx: array
        Indices of points that are NOT outliers
    groups: array
        Group numbers. Outliers have group number ``-1``.
    labels: array
        Label for each point in color, ref_color with the group
        number that it is a member of
    """

    try:
        from sklearn.cluster import DBSCAN
    except ImportError:
        raise ImportError("You must have sklearn installed to clip outliers")

    coords = np.array([color, ref_color])
    coords = coords.T
    db = DBSCAN(**dbscan_kwargs).fit(coords)
    groups = np.unique(db.labels_)
    if show_plots:
        import matplotlib.pyplot as plt
        # for gid in groups:
        #    idx = gid==db.labels_
        #    plt.plot(color[idx],ref_color[idx], '.')
        idx = db.labels_ >= 0
        plt.plot(color[~idx], ref_color[~idx], 'r.')
        plt.plot(color[idx], ref_color[idx], '.')
        plt.title('Outliers in red')
        plt.xlabel('instrumental color')
        plt.ylabel('reference color')
        plt.show()

    idx = db.labels_ > -1
    return idx, groups, db.labels_

# -----------------------------------------------------------------

def calibrate_color(instr_color, airmass, a, b, k1, k2):
    """
    Transform colors to a different photometric system.
    See Landolt 2007 for more.
    """
    return a + b * (1 / (1 + k2 * airmass)) * (instr_color - k1 * airmass)

# -----------------------------------------------------------------

def calculate_color_coeffs(instr_color, ref_color, airmass, init_params):
    """
    Using the procedure in Landolt 2007 we adjust the colors to a set
    of reference colors.
    """
    from scipy.optimize import curve_fit

    def get_color(measurements, a, b, k1, k2):
        # Unpack the observations and call the calibrate_color function
        # (for consistency)
        instr_color, airmass = measurements
        return calibrate_color(instr_color, airmass, a, b, k1, k2)

    results = curve_fit(get_color, [instr_color, airmass], ref_color,
                        init_params)
    return results

# -----------------------------------------------------------------

def calibrate_magnitude(instr_mag, airmass, ref_color,
                        zero, extinct, color, instr=None):
    """
    Calibrate instrumental magnitude to a standard photometric system.
    This assumes that ref_color has already been adjusted for any
    non-linearities between the instrumental magnitudes and the
    photometric system. Otherwise use `~calibrate_magnitude_full`
    """
    if instr is None:
        result = instr_mag - extinct * airmass + zero + color * ref_color
    else:
        result = instr * (instr_mag - extinct * airmass) + zero + color * ref_color
    return result

# -----------------------------------------------------------------

def calibrate_magnitude_full(instr_mag, airmass, instr_color,
                             a, b, k1, k2, zero, extinct, color, instr=None):
    """
    Using the transformation in Landolt 2007, calibrate an instrumental
    magnitude to a photometric standard.
    """
    adjusted_color = calibrate_color(instr_color, airmass, a, b, k1, k2)
    result = calibrate_magnitude(instr_mag, airmass, adjusted_color,
                                 zero, extinct, color, instr)
    return result

# -----------------------------------------------------------------

def calculate_izY_coeffs(instr_mag, instr_color, airmass,
                         ref_mag, ref_color, dbscan_kwargs={}, show_plots=True,
                         color_init_params=None, mag_init_params=None, cluster=True):
    """
    Use the procedure from Landolt 2007 to calibrate to a given
    standard catalog (including adjusting the instrumental colors to the
    standard catalog colors).
    """
    from scipy.optimize import curve_fit

    def get_mag(measurements, zero, extinct, color, instr=None):
        instr_mag, airmass, ref_color = measurements
        result = calibrate_magnitude(instr_mag, airmass, ref_color,
                                     zero, extinct, color, instr)
        return result

    # Calculate the coefficients to adjust to the standard colors
    if color_init_params is None:
        color_init_params = [2., 1., .1, .1]

    if cluster:
        idx, groups, labels = clip_color_outliers(
            instr_color, ref_color, dbscan_kwargs, show_plots)
        color_result = calculate_color_coeffs(
            instr_color[idx], ref_color[idx], airmass[idx], color_init_params)
    else:
        color_result = calculate_color_coeffs(
            instr_color, ref_color, airmass, color_init_params)
    a, b, k1, k2 = color_result[0]
    adjusted_color = calibrate_color(instr_color, airmass, a, b, k1, k2)

    if show_plots:
        import matplotlib.pyplot as plt
        for am in np.unique(airmass):
            aidx = airmass[idx] == am
            x = np.linspace(np.min(instr_color[idx][aidx]),
                            np.max(instr_color[idx][aidx]), 10)
            y = calibrate_color(x, am, a, b, k1, k2)
            plt.plot(instr_color[idx][aidx], ref_color[idx][aidx],
                     '.', alpha=.1)
            plt.plot(x, y, 'r')
            plt.title('Airmass={0:.2f}'.format(am))
            plt.xlabel('Instrumental Color')
            plt.ylabel('Reference Color')
            plt.show()
        plt.plot(adjusted_color, ref_color, '.', alpha=.1)
        plt.title('Calibrated Color')
        plt.xlabel('Adjusted Color')
        plt.ylabel('Standard Color')
        plt.axis('equal')
        plt.show()

    # Fit the coefficients
    measurements = [instr_mag, airmass, adjusted_color]
    if mag_init_params is None:
        mag_init_params = [25., .1, .1]
    mag_result = curve_fit(get_mag, measurements, ref_mag, mag_init_params)
    # Package the results
    if len(mag_result[0]) == 3:
        zero, extinct, color = mag_result[0]
        instr = None
    else:
        zero, extinct, color, instr = mag_result[0]
    results = color_result[0].tolist() + mag_result[0].tolist()

    if show_plots:
        mag = calibrate_magnitude_full(instr_mag, airmass, instr_color,
                                       a, b, k1, k2, zero, extinct, color, instr)
        diff = mag - ref_mag
        rms = np.sqrt(np.mean(diff) ** 2 + np.std(diff) ** 2)
        plt.plot(mag, diff, '.', alpha=.1)
        plt.title('Calibrated magnitudes, rms={0:.4f}'.format(rms))
        plt.ylim([-.15, .15])
        plt.xlabel('mag')
        plt.ylabel('diff from standard')
        plt.show()
    return results, color_result[1], mag_result[1]

# -----------------------------------------------------------------

def calculate_coeffs_by_frame(instr_mag, instr_color, airmass,
                              ref_mag, ref_color, catalog_frame, frames,
                              dbscan_kwargs={}, color_init_params=None,
                              mag_init_params=[25., .1, .1], show_plots=True):
    """
    Calculate coefficients to transform instrumental magnitudes
    to a standard photometric catalog individually for each frame.
    """
    # Clip outliers from the entire catalog
    idx, groups, labels = clip_color_outliers(
        instr_color, ref_color, dbscan_kwargs, show_plots)
    mag = np.zeros((np.sum(idx),))
    mag[:] = np.nan

    # Create a table to hold the coefficients
    frame_count = len(frames)
    a = np.zeros((frame_count,), dtype=float)
    b = np.zeros((frame_count,), dtype=float)
    k1 = np.zeros((frame_count,), dtype=float)
    k2 = np.zeros((frame_count,), dtype=float)
    zero = np.zeros((frame_count,), dtype=float)
    color = np.zeros((frame_count,), dtype=float)
    extinct = np.zeros((frame_count,), dtype=float)
    instr = np.zeros((frame_count,), dtype=float)
    frame_coeff = np.zeros((frame_count,), dtype='S4')

    # For each frame calculate the coefficients
    for n, frame in enumerate(frames):
        fidx = catalog_frame[idx] == frame
        result = calculate_izY_coeffs(
            instr_mag[idx][fidx], instr_color[idx][fidx], airmass[idx][fidx],
            ref_mag[idx][fidx], ref_color[idx][fidx],
            dbscan_kwargs, show_plots=False, cluster=False,
            color_init_params=color_init_params,
            mag_init_params=mag_init_params)
        if len(mag_init_params) == 3:
            a[n], b[n], k1[n], k2[n], zero[n], extinct[n], color[n] = result[0]
        else:
            a[n], b[n], k1[n], k2[n], zero[n], extinct[n], color[n], instr[n] = result[0]
        frame_coeff[n] = frame
    # Build the table
    if len(mag_init_params) == 3:
        result = Table([a, b, k1, k2, zero, extinct, color, frame_coeff],
                       names=('a', 'b', 'k1', 'k2', 'zero', 'extinct', 'color', 'frame'))
    else:
        result = Table([a, b, k1, k2, zero, extinct, color, instr, frame_coeff],
                       names=('a', 'b', 'k1', 'k2', 'zero', 'extinct', 'color', 'instr', 'frame'))
    return result

# -----------------------------------------------------------------

def calibrate_photometry_by_frame(instr_mag, instr_color, airmass,
                                  catalog_frame, coeffs):
    """
    Transform instrumental magnitudes to a standard photometric
    catalog using different coefficients for each frame
    """
    catalog = Table([catalog_frame, np.arange(len(instr_mag), dtype=int)],
                    names=('frame', 'index'))
    joint_tbl = join(catalog, coeffs)
    joint_tbl.sort('index')
    if 'instr' in coeffs.columns.keys():
        mag = calibrate_magnitude_full(
            instr_mag, airmass, instr_color,
            joint_tbl['a'], joint_tbl['b'], joint_tbl['k1'],
            joint_tbl['k2'], joint_tbl['zero'],
            joint_tbl['extinct'], joint_tbl['color'],
            joint_tbl['instr'])
    else:
        mag = calibrate_magnitude_full(
            instr_mag, airmass, instr_color,
            joint_tbl['a'], joint_tbl['b'], joint_tbl['k1'],
            joint_tbl['k2'], joint_tbl['zero'],
            joint_tbl['extinct'], joint_tbl['color'])
    return mag

# -----------------------------------------------------------------
