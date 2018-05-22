#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.finder Contains the SourceFinder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.convolution.kernels import Gaussian2DKernel

# Import the relevant PTS classes and modules
from .extended import ExtendedSourceFinder
from .point import PointSourceFinder
from .other import OtherSourceFinder
from ..basics.mask import Mask
from ..tools import wavelengths
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..core.dataset import DataSet
from ..region.list import SkyRegionList
from ..core.image import Image
from ..catalog.extended import ExtendedSourceCatalog
from ..catalog.point import PointSourceCatalog
from ..catalog.fetcher import CatalogFetcher
from .list import GalaxyList, StarList
from ..object.galaxy import Galaxy
from ..object.star import Star
from ...core.data.sed import ObservedSED
from ...core.filter.broad import BroadBandFilter
from ...core.tools.parallelization import ParallelTarget
from ..core.frame import Frame
from ..tools import statistics
from ...core.tools.stringify import tostr
from ..convolution.kernels import get_fwhm, has_variable_fwhm, get_average_variable_fwhm
from ...core.tools import filesystem as fs
from .tables import FWHMTable, StarTable, GalaxyTable
from ..region.rectangle import SkyRectangleRegion
from ..basics.coordinate import SkyCoordinate
from ..basics.stretch import SkyStretch

# -----------------------------------------------------------------

# FROM HELGA
# [A nearby galaxy] It is consequently contaminated by the light of thousands
# of foreground stars, especially in the UV and optical part
# of the spectrum. At longer wavelengths, the infrared emission of
# background galaxies becomes the main source of contaminating
# sources. At Herschel wavelengths, however, most of the emission
# from non-[the galaxy] point sources is negligible even at scales
# of the SPIRE 500 µm beam.
#
# The extended emission of the Milky Way Galactic Cirrus is prominently visible
# here. This dust emission can fortunately be associated with
# HI emission. Using the velocity information of HI maps, the
# Galactic cirrus can partly be disentangled from the emission of
# [the galaxy]. Paper I [of HELGA] goes into more detail about this technique.
#
# We made use of SExtractor v2.8.6 (Bertin & Arnouts 1996)
# to list the location of all point sources above a certain threshold
# (5 times the background noise level). The program simultaneously
# produces background maps that can be tweaked to represent
# the diffuse emission from M 31. In this way, we could
# replace the non-M 31 point sources with the local M 31 background
# value obtained from these maps.
#
# For each source an optimal radius was derived by comparing
# the pixel flux with the local background at increasing distance
# from the peak location. Once the pixel-to-background flux
# ratio dropped below 2, the radius was cut off at that distance.
# Based on this radius, a total flux was extracted in order to make
# colour evaluations. We constructed point source masks for the
# GALEX, SDSS, WISE, and Spitzer subsets based on different
# colour criteria.

# The GALEX and SDSS point sources were evaluated based
# on their UV colour. This technique was applied by Gil de Paz
# et al. (2007) for over 1000 galaxies and proved successful. In
# practice, we mask all sources with

# |FUV-NUV| > 0.75

# if they are detected at the 1σ level in their particular wavelength
# band. SExtractor identified 58 330 point sources in the UV fields,
# of which over 51 000 were masked in the FUV and NUV. Many
# point sources from the UV catalogue were not detected at optical
# bands, hence only 25 000 sources were masked in the SDSS
# bands. Around 7000 sources were identified as extragalactic.
# They were therefore assumed to belong to M 31 and were not
# masked.

# As an example, Fig. A.1 shows the u-band image of M 31
# before and after the mask was applied. The contamination of the
# image has been significantly reduced using the above technique.
# The point sources in the WISE and the Spitzer IRAC
# and MIPS frames were masked analogously, based on their

# IRAC colours (see below). At these wavelengths, however, the
# non-M 31 point sources are a mix of foreground stars and background
# galaxies. Furthermore, some bright sources may be associated
# with HII regions in M 31 and must not be masked. We designed
# a scheme based on the technique by Muñoz-Mateos et al.
# (2009b), which was successfully applied to the SINGS galaxies.
# Foreground stars have almost no PAH emission, while the diffuse
# ISM in galaxies shows a roughly constant F5.8/F8 ratio (Draine
# & Li 2007). Background galaxies are redshifted spirals or ellipticals
# and can consequently have a wide range in F5.8/F8. It is
# thus possible to construct a rough filter relying on the difference
# in MIR flux ratios. First, it was checked which point source extracted
# from the IRAC 3.6 µm had a non-detection at 8 µm. This
# riterion proved to be sufficient to select the foreground stars
# in the field. A second, colour-based, criterion disentangled the
# background galaxies from the HII regions:

# 0.29 < F5.8 / F8 < 0.85
# F3.6 / F5.8 < 1.58

# Figure A.2 shows the colour−colour diagram for these sources.
# The HII regions follow a more or less horizontal track at the
# lower-left part of the plot. The colour criteria for filtering out
# these HII regions were obtained empirically to ensure effective
# identification. Once identified, these star forming regions were
# consequently not masked. The resulting mask was applied to all
# IRAC and MIPS bands. Sources that were not detected at longer
# wavelengths were obviously not masked. Figure A.1 shows the
# IRAC 3.6 µm image of M 31 before and after the mask was applied.

# -----------------------------------------------------------------

class SourceFinder(Configurable):

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
        super(SourceFinder, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frames
        self.frames = dict()

        # The error maps
        self.error_maps = dict()

        # The masks
        self.special_masks = dict()
        self.ignore_masks = dict()

        # Downsampled images
        self.downsampled = None
        self.original_wcs = None

        # The catalog fetcher
        self.fetcher = CatalogFetcher()

        # The catalog of extended sources and the catalog of point sources
        self.extended_source_catalog = None
        self.point_source_catalog = None

        # Ignore images
        self.ignore = []
        self.ignore_stars = []
        self.ignore_other_sources = []

        # The regions covering areas that should be ignored throughout the entire extraction procedure
        self.special_region = None
        self.ignore_region = None

        # The name of the principal galaxy
        self.galaxy_name = None

        ### from the extended/point source finders

        # The tables
        self.extended_tables = dict()
        self.point_tables = dict()

        # The regions
        self.extended_regions = dict()
        self.point_regions = dict()
        self.saturation_regions = dict()
        self.other_regions = dict()

        # The segmentation maps
        self.segments = dict()

        ###

        ### The finished products:

        # The tables
        self.galaxy_table = None
        self.star_table = None

        # The regions
        #self.galaxy_regions = None
        #self.star_regions = None
        #self.saturation_regions = None

        # The segmentation maps
        self.galaxy_segments = None
        self.star_segments = None

        ###

        # Settings for the star finder for different bands
        self.star_finder_settings = dict()

        # Extended sources and point sources
        self.extended_sources = dict()
        self.point_sources = dict()

        # Principal masks and companion masks
        self.principal_masks = dict()
        self.companion_masks = dict()

        # Galaxy and star lists
        self.galaxies = GalaxyList()
        self.stars = StarList()

        # The FWHMs and PSFs
        self.fwhms = dict()
        self.psfs = dict()

        # The photometry table
        self.photometry = None

        # Output paths for image names
        self.output_paths = None

    # -----------------------------------------------------------------

    def add_frame(self, name, frame, star_finder_settings=None, error_map=None, output_path=None):

        """
        This function ...
        :param name:
        :param frame:
        :param star_finder_settings:
        :param error_map:
        :param output_path:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

        # If error map is given
        if error_map is not None: self.error_maps[name] = error_map

        # Set the settings
        if star_finder_settings is not None: self.star_finder_settings[name] = star_finder_settings

        # If output path is given
        if output_path is not None:

            if self.output_paths is None: self.output_paths = dict()
            self.output_paths[name] = output_path

    # -----------------------------------------------------------------

    def add_error_map(self, name, error_map):

        """
        This function ...
        :param name:
        :param error_map:
        :return:
        """

        # Check if name in frames
        if name not in self.frames: raise ValueError("Frame with the name '" + name + "' has not been added")

        # Set the error map
        self.error_maps[name] = error_map

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

    @property
    def overlap_box(self):

        """
        This function ...
        :return: 
        """

        # Edges
        min_ra = None
        max_ra = None
        min_dec = None
        max_dec = None

        # Loop over the frames
        for name in self.frames:

            # Get coordinate range
            #center, ra_span, dec_span = self.frames[name].wcs.coordinate_range
            #radius = SkyStretch(0.5 * ra_span, 0.5 * dec_span)
            #box = SkyRectangleRegion(center, radius)

            # Get wcs
            wcs = self.frames[name]

            if min_ra is None or wcs.min_ra > min_ra: min_ra = wcs.min_ra
            if max_ra is None or wcs.max_ra < max_ra: max_ra = wcs.max_ra
            if min_dec is None or wcs.min_dec > min_dec: min_dec = wcs.min_dec
            if max_dec is None or wcs.max_dec < max_dec: max_dec = wcs.max_dec

        # Determine center and radius
        # Get center and radius of the new bounding box
        center = SkyCoordinate(0.5 * (min_ra + max_ra), 0.5 * (min_dec + max_dec))
        radius = SkyStretch(0.5 * (max_ra - min_ra), 0.5 * (max_dec - min_dec))

        # Return the bounding box
        return SkyRectangleRegion(center, radius)

    # -----------------------------------------------------------------

    @property
    def catalog_coordinate_box(self):

        """
        This function ...
        :return: 
        """

        # Determine the bounding box
        if self.config.catalog_overlapping: coordinate_box = self.overlap_box
        else: coordinate_box = self.bounding_box

        # Return
        return coordinate_box

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return:
        """

        return [frame.filter for frame in self.frames.values()]

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 3. Find the galaxies
        if self.config.find_galaxies: self.find_galaxies()
        
        # 3. Find the stars
        if self.config.find_stars: self.find_stars()

        # 4. Look for other sources
        if self.config.find_other_sources: self.find_other_sources()

        # 5. Perform the photometry
        #self.do_photometry()

        # Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SourceFinder, self).setup(**kwargs)

        # Load the images (from config or input kwargs)
        if "frames" in kwargs:
            self.frames = kwargs.pop("frames")
            if "error_maps" in kwargs: self.error_maps = kwargs.pop("error_maps")
        elif "dataset" in kwargs:
            dataset = kwargs.pop("dataset")
            self.frames = dataset.get_frames()
            self.error_maps = dataset.get_errormaps()
        else: self.load_frames()

        # Get the output path
        if "output_paths" in kwargs: self.output_paths = kwargs.pop("output_paths")

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

        # Initialize the galaxy segments
        self.galaxy_segments = Image("galaxies")

        # Initialize the star segments
        self.star_segments = Image("stars")

        # Initialize the galaxy table
        self.galaxy_table = GalaxyTable(filters=self.filters)

        # Initialize the star table
        self.star_table = StarTable(filters=self.filters)

        # Load the catalogs
        self.load_catalogs(**kwargs)

    # -----------------------------------------------------------------

    def load_catalogs(self, **kwargs):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the catalogs ...")

        # Load extended sources catalog
        self.load_extended_sources_catalog(kwargs)

        # Load point sources catalog
        self.load_point_sources_catalog(kwargs)

    # -----------------------------------------------------------------

    def load_extended_sources_catalog(self, kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Loading catalog of extended sources ...")

        # From kwargs
        if "extended_source_catalog" in kwargs and kwargs["extended_source_catalog"] is not None: self.extended_source_catalog = kwargs.pop("extended_source_catalog")

        # From file
        elif self.config.extended_sources_catalog is not None: self.extended_source_catalog = ExtendedSourceCatalog.from_file(self.config.extended_sources_catalog)

        # Not specified
        else: pass

    # -----------------------------------------------------------------

    def load_point_sources_catalog(self, kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading catalog of point sources ...")

        # From kwargs
        if "point_source_catalog" in kwargs and kwargs["point_source_catalog"] is not None: self.point_source_catalog = kwargs.pop("point_source_catalog")

        # From file
        elif self.config.point_sources_catalog is not None: self.point_source_catalog = PointSourceCatalog.from_file(self.config.point_sources_catalog)

        # Not specified
        else: pass

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

            # Load the image
            image = Image.from_file(self.config.dataset)

            # Determine the name for this image
            name = str(image.filter)

            # Add the primary frame
            self.add_frame(name, image.primary)

            # Add the error map
            if image.has_errors: self.add_error_map(name, image.errors)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the frames
            self.frames = dataset.get_frames()

            # Get the error maps
            self.error_maps = dataset.get_errormaps()

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

        # Special mask
        if self.special_region is not None:

            # Loop over the frames
            for name in self.frames:

                # Create pixel region for this frame
                pixel_region = self.special_region.to_pixel(self.frames[name].wcs)

                # Create the mask
                special_mask = Mask.from_region(pixel_region, self.frames[name].xsize, self.frames[name].ysize)

                # Set the mask
                self.special_masks[name] = special_mask

        # Ignore mask
        if self.ignore_region is not None:

            # Loop over the frames
            for name in self.frames:

                # Create pixel region for this frame
                pixel_region = self.ignore_region.to_pixel(self.frames[name].wcs)

                # Create the mask
                ignore_mask = Mask.from_region(self.ignore_region, self.frames[name].xsize, self.frames[name].ysize)

                # Set the mask
                self.ignore_masks[name] = ignore_mask

    # -----------------------------------------------------------------

    def find_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the galaxies ...")

        # Fetch catalog of extended sources if necessary
        if self.extended_source_catalog is None: self.fetch_extended_sources_catalog()

        # Find extended sources
        self.find_extended_sources()

        # Make list of galaxies
        self.collect_galaxies()

        # Create galaxy table
        self.create_galaxy_table()

    # -----------------------------------------------------------------

    def fetch_extended_sources_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching catalog of extended sources ...")

        # Fetch the catalog
        self.extended_source_catalog = self.fetcher.get_extended_source_catalog(self.catalog_coordinate_box)

    # -----------------------------------------------------------------

    def find_extended_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding extended sources ...")

        # Dictionary to keep result handles
        results = dict()

        # Parallel execution
        with ParallelTarget(detect_extended_sources, self.config.nprocesses) as target:
        #with ParallelTarget(detect_extended_sources_wrapper, self.config.nprocesses) as target:

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
                config = self.config.extended.copy()

                # Save temporarily
                #temp_path = fs.create_directory_in(introspection.pts_temp_dir, time.unique_name("sourcefinder_" + name))
                #frame_path = fs.join(temp_path, "frame.fits")
                #special_mask_path = fs.join(temp_path, "special_mask.fits") if special_mask is not None else None
                #ignore_mask_path = fs.join(temp_path, "ignore_mask.fits") if ignore_mask is not None else None
                #bad_mask_path = fs.join(temp_path, "bad_mask.fits") if bad_mask is not None else None

                # Save
                #frame.saveto(frame_path)
                #if special_mask is not None: special_mask.saveto(special_mask_path)
                #if ignore_mask is not None: ignore_mask.saveto(ignore_mask_path)
                #if bad_mask is not None: bad_mask.saveto(bad_mask_path)

                #print(self.config)

                # Weak search
                config["weak"] = self.config.weak

                # Set the output path
                config["output"] = self.config.output

                # Do the detection
                # Call the target function
                result = target(frame, self.extended_source_catalog, config, special_mask, ignore_mask, bad_mask)
                #result = target(frame_path, self.extended_source_catalog, config, special_mask_path, ignore_mask_path, bad_mask_path)

                # Set the result handle
                results[name] = result

        # Process results
        for name in results:

            # Request output
            results[name].request()

            # Get result
            table, regions, segments, principal_mask, companion_mask = results[name].output

            # Set galaxies
            self.extended_tables[name] = table

            # Set region list
            self.extended_regions[name] = regions

            # Set segmentation map
            # Add the segmentation map of the galaxies
            self.segments[name].add_frame(segments, "extended")

            # Set mask of principal galaxy
            self.principal_masks[name] = principal_mask

            # Set mask of companion galaxies
            self.companion_masks[name] = companion_mask

    # -----------------------------------------------------------------

    def collect_galaxies(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Collecting galaxies ...")

        # Loop over the extended sources
        for index in range(len(self.extended_source_catalog)):

            # Get the galaxy position
            position = self.extended_source_catalog.get_position(index)

            # Create SED
            sed = ObservedSED(photometry_unit="Jy")

            # Loop over the frames
            for name in self.frames:

                # Get the flux for this frame
                flux = self.extended_tables[name].get_flux(index)

                # Add the flux to the SED
                if flux is not None: sed.add_point(self.frames[name].filter, flux)

            if len(sed) == 0: sed = None

            # Get other properties
            name = self.extended_source_catalog.get_name(index)
            redshift = self.extended_source_catalog.get_redshift(index)
            galaxy_type = self.extended_source_catalog.get_type(index)
            names = self.extended_source_catalog.get_names(index)
            distance = self.extended_source_catalog.get_distance(index)
            inclination = self.extended_source_catalog.get_inclination(index)
            d25 = self.extended_source_catalog.get_d25(index)
            major = self.extended_source_catalog.get_major(index)
            minor = self.extended_source_catalog.get_minor(index)
            posangle = self.extended_source_catalog.get_position_angle(index)
            principal = self.extended_source_catalog.is_principal(index)
            companions = self.extended_source_catalog.get_companions(index)
            parent = self.extended_source_catalog.get_parent(index)

            # Create the galaxy
            galaxy = Galaxy(index=index, position=position, sed=sed, name=name, redshift=redshift, galaxy_type=galaxy_type,
                            names=names, distance=distance, inclination=inclination, d25=d25, major=major, minor=minor,
                            position_angle=posangle, principal=principal, companions=companions, parent=parent)

            # Add the galaxy to the list
            self.galaxies.append(galaxy)

    # -----------------------------------------------------------------

    def create_galaxy_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the galaxy table ...")

        # Add the galaxies
        for galaxy in self.galaxies: self.galaxy_table.add_galaxy(galaxy)

    # -----------------------------------------------------------------

    def find_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        log.info("Finding the stars ...")

        # Fetch catalog if necessary
        if self.point_source_catalog is None: self.fetch_point_sources_catalog()

        # Find point sources
        self.find_point_sources()

        # Collect stars
        self.collect_stars()

        # Create the star table
        self.create_star_table()

    # -----------------------------------------------------------------

    def fetch_point_sources_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching catalog of point sources ...")

        # Get minimum pixelscale
        min_pixelscale = self.min_pixelscale

        # Fetch
        self.point_source_catalog = self.fetcher.get_point_source_catalog(self.catalog_coordinate_box, min_pixelscale, self.config.point.fetching.catalogs)

    # -----------------------------------------------------------------

    def find_point_sources(self):

        """
        This function ...
        :return:
        """

        # Dictionary to keep result handles
        results = dict()

        # Parallel execution
        with ParallelTarget(detect_point_sources, self.config.nprocesses) as target:

            # Loop over the images
            for name in self.frames:

                # Ignore if requested
                if name in self.ignore: continue
                if name in self.ignore_stars: continue

                # Get the frame
                frame = self.frames[name]

                # Don't run the star finder if the wavelength of this image is greater than 25 micron
                if frame.wavelength is None or frame.wavelength > wavelengths.ranges.ir.mir.max:

                    # No star subtraction for this image
                    log.info("Finding point sources will not be performed for the '" + name + "' image")
                    continue

                # Inform the user
                log.info("Finding the point sources for the '" + name + "' image ...")

                # Get masks
                special_mask = self.special_masks[name] if name in self.special_masks else None
                ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
                bad_mask = None

                # Create configuration
                config = self.config.point.copy()
                if name in self.star_finder_settings: config.set_items(self.star_finder_settings[name])

                #print(self.config)

                # Weak search
                config["weak"] = self.config.weak

                # Set the output path
                config["output"] = self.config.output

                # Call the target function
                result = target(frame, self.galaxies, self.point_source_catalog, config, special_mask, ignore_mask, bad_mask, self.principal_masks[name])

                # Add the result
                results[name] = result

        # Process results
        for name in results:

            # Request
            results[name].request()

            # Get result
            # stars, star_region_list, saturation_region_list, star_segments, kernel, statistics
            #stars, star_region_list, saturation_region_list, star_segments, kernel, statistics = results[name].get()
            table, regions, saturation_regions, segments, fwhm = results[name].output

            # Set table
            self.point_tables[name] = table

            # Set regions
            self.point_regions[name] = regions
            self.saturation_regions[name] = saturation_regions

            # Set segmentation map
            # Add the segmentation map of the galaxies
            self.segments[name].add_frame(segments, "point")

            # Set the PSF
            if fwhm is None: fwhm = self.frames[name].fwhm
            self.fwhms[name] = fwhm

            # Create a Gaussian convolution kernel and return it
            if fwhm is not None:

                # Debugging
                log.debug("The FWHM of the '" + name + "' image is " + tostr(fwhm))

                fwhm_pix = (fwhm / self.frames[name].average_pixelscale.to("arcsec")).value if fwhm is not None else None

                sigma = fwhm_pix * statistics.fwhm_to_sigma
                kernel = Gaussian2DKernel(sigma)

                # Set the PSF
                self.psfs[name] = kernel

            # Get the statistics
            #self.statistics[name] = point_source_statistics

            # Show the FWHM
            #log.info("The FWHM that could be fitted to the point sources in the " + name + " image is " + str(self.statistics[name].fwhm))

        # Set the kernels of the frames for which point sources were not searched
        for name in self.frames:

            if name in self.psfs: continue
            fwhm = self.frames[name].fwhm

            # Get the FWHM if not defined
            if fwhm is None:

                if has_variable_fwhm(self.frames[name].filter):
                    log.warning("The FWHM of the '" + name + "' image is still undefined after point source detection and FWHM for this filter is variable. Using average value (Clark et al., 2017) to proceed.")
                    fwhm = get_average_variable_fwhm(self.frames[name].filter)
                else: fwhm = get_fwhm(self.frames[name].filter)

            # Debugging
            log.debug("The FWHM of the '" + name + "' image is " + tostr(fwhm))

            # Set the FWHM
            self.fwhms[name] = fwhm

            # Create kernel
            fwhm_pix = (fwhm / self.frames[name].average_pixelscale.to("arcsec")).value
            sigma = fwhm_pix * statistics.fwhm_to_sigma
            kernel = Gaussian2DKernel(sigma)
            self.psfs[name] = kernel

    # -----------------------------------------------------------------

    def collect_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Collecting stars ...")

        # Loop over the point sources
        for index in range(len(self.point_source_catalog)):

            # Get the position of the point source
            position = self.point_source_catalog.get_position(index)

            # Create SED
            sed = ObservedSED(photometry_unit="Jy")

            # Loop over the frames
            for name in self.frames:

                # Point source detection was not performed for this image
                if name not in self.point_tables: continue

                # Get the flux for this frame
                flux = self.point_tables[name].get_flux(index)

                # Add the flux to the SED
                if flux is not None: sed.add_point(self.frames[name].filter, flux)

            # Check whether it can be identified as a star
            # If the SED cannot correspond to a star, skip this source
            if not is_stellar_sed(sed): continue

            # Get other properties
            catalog = self.point_source_catalog.get_catalog(index)
            id = self.point_source_catalog.get_id(index)
            ra_error = self.point_source_catalog.get_ra_error(index)
            dec_error = self.point_source_catalog.get_dec_error(index)

            # Create FWHM table
            fwhms = FWHMTable()

            # Loop over the frames
            for name in self.frames:

                # If point sources not searched in this frame
                if name not in self.point_tables: continue

                # Get the FWHM
                fwhm = self.point_tables[name].get_fwhm(index)

                # Add an entry to the FWHM table
                if fwhm is not None: fwhms.add_fwhm(self.frames[name].filter, fwhm)

            # Create the star object
            star = Star(index=index, position=position, sed=sed, catalog=catalog, id=id, ra_error=ra_error, dec_error=dec_error, fwhms=fwhms)

            # Add the star to the list
            self.stars.append(star)

    # -----------------------------------------------------------------

    def create_star_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the star table ...")

        # Loop over the stars
        for star in self.stars: self.star_table.add_star(star)

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

        # Parallel execution
        with ParallelTarget(detect_other, self.config.nprocesses) as target:

            # Loop over the frames
            for name in self.frames:

                # Ignore if requested
                if name in self.ignore: continue
                if name in self.ignore_other_sources: continue

                # Get the frame
                frame = self.frames[name]

                # If the wavelength of this image is greater than 25 micron, don't classify the sources that are found
                #if frame.wavelength is not None and frame.wavelength > wavelengths.ranges.ir.mir.max: self.trained_finder.config.classify = False
                #else: self.trained_finder.config.classify = True

                # Get masks
                special_mask = self.special_masks[name] if name in self.special_masks else None
                ignore_mask = self.ignore_masks[name] if name in self.ignore_masks else None
                bad_mask = None

                # Create the configuration
                config = self.config.other.copy()

                # Get other input
                #galaxies = self.galaxies[name]
                #stars = self.stars[name]

                # Get the segmentation maps
                galaxy_segments = self.segments[name].frames["extended"]
                star_segments = self.segments[name].frames["point"] if "point" in self.segments[name].frames else None

                # Get the PSF kernel
                kernel = self.psfs[name]

                # Call the target function
                result = target(frame, config, self.galaxies, self.stars, galaxy_segments, star_segments, kernel,
                                special_mask, ignore_mask, bad_mask, self.principal_masks[name], self.companion_masks[name])

                # Add the result
                results[name] = result

        # Process results
        for name in results:

            # Request
            results[name].request()

            # Get the result
            region_list, segments = results[name].output

            # Add the region
            self.other_regions[name] = region_list

            # Add the segmentation map of the other sources
            self.segments[name].add_frame(segments, "other")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the catalogs
        self.write_catalogs()

        # 2. Write the tables
        self.write_tables()

        # 3. Write region lists
        self.write_regions()

        # 4. Write segmentation maps
        self.write_segments()

    # -----------------------------------------------------------------

    def write_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing catalogs ...")

        # Write extended sources catalog
        self.write_extended_source_catalog()

        # Write point sources catalog
        self.write_point_source_catalog()

    # -----------------------------------------------------------------

    def write_extended_source_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the catalog of extended sources ...")

        # Write
        path = self.output_path_file("extended_sources.cat")
        self.extended_source_catalog.saveto(path)

    # -----------------------------------------------------------------

    def write_point_source_catalog(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the catalog of point sources ...")

        # Write
        path = self.output_path_file("point_sources.cat")
        self.point_source_catalog.saveto(path)

    # -----------------------------------------------------------------

    def write_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing tables ...")

        # Write galaxy table
        self.write_galaxy_table()

        # Write star table
        self.write_star_table()

    # -----------------------------------------------------------------

    def write_galaxy_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing galaxy table ...")

        # Write
        path = self.output_path_file("galaxies.dat")
        self.galaxy_table.saveto(path)

    # -----------------------------------------------------------------

    def write_star_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing star table ...")

        # Write
        path = self.output_path_file("stars.dat")
        self.star_table.saveto(path)

    # -----------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the regions ...")

        # Write galaxy regions
        self.write_galaxy_regions()

        # Write star regions
        self.write_star_regions()

        # Write saturation regions
        self.write_saturation_regions()

        # Write other regions
        self.write_other_regions()

    # -----------------------------------------------------------------

    def find_extended_region(self, name, index):

        """
        This function ...
        :param name: 
        :param index: 
        :return: 
        """

        for region in self.extended_regions[name]:
            #if int(region.meta["text"]) == index: return region
            #print(region.meta)
            if "index" not in region.meta: continue
            if region.meta["index"] == index: return region
        return None

    # -----------------------------------------------------------------

    def write_galaxy_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the galaxy regions ...")

        # Loop over the extended regions
        #for name in self.extended_regions:

        # Loop over the images
        for name in self.extended_regions:

            # Initialize region list
            #galaxy_regions = PixelRegionList()
            galaxy_regions = SkyRegionList()

            # Loop over the galaxies
            for galaxy in self.galaxies:

                # Get the index
                index = galaxy.index

                #print("index", index)

                # Get the corresponding region
                #region = self.extended_regions[name][index]
                region = self.find_extended_region(name, index)

                #print("region", region)

                # Add the region
                galaxy_regions.append(region)

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "galaxies.reg")
            else: path = self.output_path_file("galaxies_" + name + ".reg")

            # Save
            galaxy_regions.saveto(path)

    # -----------------------------------------------------------------

    def find_star_region(self, name, index):

        """
        This function ...
        :param name: 
        :param index: 
        :return: 
        """

        for region in self.point_regions[name]:
            #print(region.meta)
            #if int(region.meta["text"]) == index: return region
            if "index" not in region.meta: continue
            if region.meta["index"] == index: return region
        #print(name, index)
        #print([region.meta["index"] for region in self.point_regions[name] if "index" in region.meta])
        return None

    # -----------------------------------------------------------------

    def write_star_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the star regions ...")

        # Loop over the images
        for name in self.point_regions:

            # Initialize region list
            #star_regions = PixelRegionList()
            star_regions = SkyRegionList()

            # Loop over the stars
            for star in self.stars:

                # Get the index
                index = star.index

                # Get the corresponding region
                #region = self.point_regions[name][index]
                region = self.find_star_region(name, index)

                # Add the region
                if region is not None: star_regions.append(region)

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "stars.reg")
            else: path = self.output_path_file("stars_" + name + ".reg")

            # Save
            star_regions.saveto(path)

    # -----------------------------------------------------------------

    def write_saturation_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the saturation regions ...")

        # Loop over the images
        for name in self.saturation_regions:

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "saturation.reg")
            else: path = self.output_path_file("saturation_" + name + ".reg")

            # Save
            if self.saturation_regions[name] is not None: self.saturation_regions[name].saveto(path)

    # -----------------------------------------------------------------

    def write_other_regions(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the regions of other sources ...")

        # Loop over the images
        for name in self.other_regions:
            
            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "other_sources.reg")
            else: path = self.output_path_file("other_sources_" + name + ".reg")

            # Save
            self.other_regions[name].saveto(path)

    # -----------------------------------------------------------------

    def write_segments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the segmentation maps ...")

        # Loop over the images
        for name in self.segments:

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "segments.fits")
            else: path = self.output_path_file("segments_" + name + ".fits")
            
            # Save
            self.segments[name].saveto(path)

# -----------------------------------------------------------------

def serialize_input(*args):

    """
    This function ...
    :param args:
    :return:
    """

    pass

# -----------------------------------------------------------------

def detect_extended_sources_wrapper(frame_path, catalog, config, special_mask_path, ignore_mask_path, bad_mask_path):

    """
    This function ...
    :param frame_path:
    :param catalog:
    :param config:
    :param special_mask_path:
    :param ignore_mask_path:
    :param bad_mask_path:
    :return:
    """

    # Open the frame
    frame = Frame.from_file(frame_path)

    # Open the masks
    special_mask = Mask.from_file(special_mask_path) if special_mask_path is not None else None
    ignore_mask = Mask.from_file(ignore_mask_path) if ignore_mask_path is not None else None
    bad_mask = Mask.from_file(bad_mask_path) if bad_mask_path is not None else None

    # Implementation
    table, regions, segments, principal_mask, companion_mask = detect_extended_sources(frame, catalog, config, special_mask, ignore_mask, bad_mask)

    # Return
    return table, regions, segments, principal_mask, companion_mask

# -----------------------------------------------------------------

def detect_extended_sources(frame, catalog, config, special_mask, ignore_mask, bad_mask):

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
    finder = ExtendedSourceFinder(config)

    # Inform the user
    log.info("Starting detection of extended sources for '" + frame.name  + "' ...")

    # Run the finder
    finder.run(frame=frame, catalog=catalog, special_mask=special_mask, ignore_mask=ignore_mask, bad_mask=bad_mask)

    # Set the name of the principal galaxy
    # self.galaxy_name = self.galaxy_finder.principal.name

    # Get the galaxy region
    # galaxy_sky_region = self.galaxy_finder.finder.galaxy_sky_region
    # if galaxy_sky_region is not None:
    #    galaxy_region = galaxy_sky_region.to_pixel(image.wcs)

    if finder.regions is not None:

        galaxy_sky_regions = finder.regions.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.galaxy_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # ele: return self.galaxy_finder.region

        # Get region list
        regions = galaxy_sky_regions

    else: regions = None

    if finder.segments is not None:

        # if self.galaxy_finder.segments is None: return None
        # if self.downsampled:

        #    segments = self.galaxy_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled

        # else: return self.galaxy_finder.segments

        # Get the segments
        segments = finder.segments

    else: segments = None

    # Get the source table
    table = finder.table

    # Inform the user
    log.success("Finished finding the extended sources for '" + frame.name + "' ...")

    # Return the output
    #return galaxies, region_list, segments

    # Get mask of principal galaxy
    principal_mask = finder.principal_mask

    # Get mask of companion galaxies
    companion_mask = finder.companion_mask

    # Return the source table, regions and segments
    return table, regions, segments, principal_mask, companion_mask

# -----------------------------------------------------------------

def detect_point_sources_wrapper(frame_path, galaxies, catalog, config, special_mask, ignore_mask, bad_mask):

    """
    This function ...
    :param frame_path:
    :param galaxies:
    :param catalog:
    :param config:
    :param special_mask:
    :param ignore_mask:
    :param bad_mask:
    :return:
    """

    # Open the frame
    frame = Frame.from_file(frame_path)

    # Call the implementation
    table, regions, saturation_regions, segments, fwhm = detect_point_sources(frame, galaxies, catalog, config, special_mask, ignore_mask, bad_mask)

    return

# -----------------------------------------------------------------

def detect_point_sources(frame, galaxies, catalog, config, special_mask, ignore_mask, bad_mask, principal_mask):

    """
    This function ...
    :param frame:
    :param galaxies:
    :param catalog:
    :param config:
    :param special_mask:
    :param ignore_mask:
    :param bad_mask:
    :param principal_mask:
    :return:
    """

    # Create the star finder
    finder = PointSourceFinder(config)

    # Inform the user
    log.info("Starting detection of point sources ...")

    # Run the finder
    finder.run(frame=frame, galaxies=galaxies, catalog=catalog, special_mask=special_mask, ignore_mask=ignore_mask,
               bad_mask=bad_mask, principal_mask=principal_mask)

    if finder.regions is not None:

        star_sky_region = finder.regions.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.star_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.star_finder.star_region

        regions = star_sky_region

    else: regions = None

    if finder.saturation_regions is not None:

        saturation_sky_region = finder.saturation_regions.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.saturation_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.star_finder.saturation_region

        saturation_regions = saturation_sky_region

    else: saturation_regions = None

    if finder.segments is not None:

        # if self.star_finder.segments is None: return None
        # if self.downsampled:
        #    segments = self.star_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled
        # else: return self.star_finder.segments

        segments = finder.segments

    else: segments = None

    # Set the stars
    #stars = finder.stars

    # Get the table
    table = finder.table

    # kernel = self.star_finder.kernel # doesn't work when there was no star extraction on the image, self.star_finder does not have attribute image thus cannot give image.fwhm
    # Set the kernel (PSF)
    #if finder.config.use_frame_fwhm and frame.fwhm is not None:
    #    fwhm = frame.fwhm.to("arcsec").value / frame.average_pixelscale.to("arcsec").value
    #    sigma = fwhm * statistics.fwhm_to_sigma
    #    kernel = Gaussian2DKernel(sigma)
    #else: kernel = finder.kernel

    # Get the FWHM
    fwhm = finder.fwhm

    # Inform the user
    log.success("Finished finding the point sources for '" + frame.name + "' ...")

    # Return the output
    #return stars, star_region_list, saturation_region_list, star_segments, kernel, statistics

    # Return the source table, regions, saturation regions, and segments
    return table, regions, saturation_regions, segments, fwhm

# -----------------------------------------------------------------

def detect_other_wrapper(frame_path, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask, bad_mask):

    """
    This function ...
    :param frame_path:
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

    # Open the frame
    frame = Frame.from_file(frame_path)

    # Implementation
    other_sky_region, other_segments = detect_other(frame, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask, bad_mask)

    return other_sky_region, other_segments

# -----------------------------------------------------------------

def detect_other(frame, config, galaxies, stars, galaxy_segments, star_segments, kernel, special_mask, ignore_mask,
                 bad_mask, principal_mask, companion_mask):

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
    :param principal_mask:
    :param companion_mask:
    :return:
    """

    # Create the other source finder
    finder = OtherSourceFinder(config)

    # Inform the user
    log.info("Starting detection of other sources ...")

    # Run the finder just to find sources
    finder.run(frame=frame, galaxies=galaxies, stars=stars, special_mask=special_mask,
                ignore_mask=ignore_mask, bad_mask=bad_mask,
                galaxy_segments=galaxy_segments,
                star_segments=star_segments, kernel=kernel, principal_mask=principal_mask, companion_mask=companion_mask)

    if finder.region is not None:

        other_sky_region = finder.region.to_sky(frame.wcs)

        # if self.downsampled:
        #    sky_region = self.other_sky_region
        #    return sky_region.to_pixel(self.original_wcs) if sky_region is not None else None
        # else: return self.trained_finder.region

    else: other_sky_region = None

    if finder.segments is not None:

        # if self.trained_finder.segments is None: return None
        # if self.downsampled:
        #    segments = self.trained_finder.segments
        #    upsampled = segments.upsampled(self.config.downsample_factor, integers=True)
        #    upsampled.unpad(self.pad_x, self.pad_y)
        #    return upsampled
        # else: return self.trained_finder.segments

        other_segments = finder.segments

    else: other_segments = None

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

def find_star_in_list(stars, index):

    """
    This function ...
    :param stars:
    :param index:
    :return:
    """

    for star in stars:

        if star.index == index: return star

    return None

# -----------------------------------------------------------------

def is_stellar_sed(sed):

    """
    This function ...
    :param sed:
    :return:
    """

    fuv = BroadBandFilter("FUV")
    nuv = BroadBandFilter("NUV")

    # Check the FUV-NUV colour, if possible
    if sed.has_filter(fuv) and sed.has_filter(nuv):

        fuv_nuv_colour = sed.colour(fuv, nuv)
        if abs(fuv_nuv_colour) < 0.75: return False

    # | FUV - NUV | > 0.75
    # if they are detected at the 1σ level in their particular wavelength
    # band.

    # IRAC colours (see below). At these wavelengths, however, the
    # non-M 31 point sources are a mix of foreground stars and background
    # galaxies. Furthermore, some bright sources may be associated
    # with HII regions in M 31 and must not be masked. We designed
    # a scheme based on the technique by Muñoz-Mateos et al.
    # (2009b), which was successfully applied to the SINGS galaxies.
    # Foreground stars have almost no PAH emission, while the diffuse
    # ISM in galaxies shows a roughly constant F5.8/F8 ratio (Draine
    # & Li 2007). Background galaxies are redshifted spirals or ellipticals
    # and can consequently have a wide range in F5.8/F8. It is
    # thus possible to construct a rough filter relying on the difference
    # in MIR flux ratios. First, it was checked which point source extracted
    # from the IRAC 3.6 µm had a non-detection at 8 µm. This
    # criterion proved to be sufficient to select the foreground stars
    # in the field. A second, colour-based, criterion disentangled the
    # background galaxies from the HII regions:

    irac_i1 = BroadBandFilter("IRAC I1")
    irac_i3 = BroadBandFilter("IRAC I3")
    irac_i4 = BroadBandFilter("IRAC I4")

    # 0.29 < F5.8 / F8 < 0.85
    # F3.6 / F5.8 < 1.58

    # Check i3 i4 color, if possible
    if sed.has_filter(irac_i3) and sed.has_filter(irac_i4):

        i3_i4_colour = sed.colour(irac_i3, irac_i4)
        if i3_i4_colour < 0.29: return False
        if i3_i4_colour > 0.85: return False

    # Check the i1 i3 color, if possible
    if sed.has_filter(irac_i1) and sed.has_filter(irac_i3):

        i1_i3_colour = sed.colour(irac_i1, irac_i3)
        if i1_i3_colour > 1.58: return False

    # All checks passed (or none were possible), return True
    return True

# -----------------------------------------------------------------
