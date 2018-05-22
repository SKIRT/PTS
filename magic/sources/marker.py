#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.marker Contains the SourceMarker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ..tools import statistics
from ...core.tools import filesystem as fs
from ..core.image import Image
from ..core.dataset import DataSet
from ..region.list import SkyRegionList
from ..core.mask import Mask
from ..catalog.extended import ExtendedSourceCatalog
from ..catalog.point import PointSourceCatalog
from ...core.filter.filter import parse_filter
from ...core.units.parsing import parse_quantity
from ..region.circle import PixelCircleRegion
from ..tools import wavelengths

# -----------------------------------------------------------------

class SourceMarker(Configurable):

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
        super(SourceMarker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frames
        self.frames = dict()

        # FWHMs
        self.fwhms = dict()

        # Ignore images
        self.ignore = []
        self.ignore_stars = []
        self.ignore_other_sources = []

        # The masks
        self.special_masks = dict()
        self.ignore_masks = dict()

        # The regions covering areas that should be ignored throughout the entire extraction procedure
        self.special_region = None
        self.ignore_region = None

        self.ignore_stars = []

        # Output paths for image names
        self.output_paths = None

        # Catalogs
        self.extended_source_catalog = None
        self.point_source_catalog = None

        # The output regions
        self.extended_regions = None
        self.point_regions = None

    # -----------------------------------------------------------------

    def add_frame(self, name, frame, output_path=None):

        """
        This function ...
        :param name:
        :param frame:
        :param output_path:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

        # If output path is given
        if output_path is not None:

            if self.output_paths is None: self.output_paths = dict()
            self.output_paths[name] = output_path

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 2. Mark extended regions
        self.create_extended_regions()

        # 3. Mark point sources
        self.create_point_regions()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(SourceMarker, self).setup(**kwargs)

        # Get FWHMS
        self.fwhms = kwargs.pop("fwhms")

        # Load the images (from config or input kwargs)
        if "frames" in kwargs:
            self.frames = kwargs.pop("frames")
            #if "error_maps" in kwargs: self.error_maps = kwargs.pop("error_maps")
        elif "dataset" in kwargs:
            dataset = kwargs.pop("dataset")
            self.frames = dataset.get_frames()
            #self.error_maps = dataset.get_errormaps()
        else: self.load_frames()

        # Get the output path
        if "output_paths" in kwargs: self.output_paths = kwargs.pop("output_paths")

        # Ignore certain images
        self.ignore = kwargs.pop("ignore", [])

        # Load special region
        self.special_region = SkyRegionList.from_file(self.config.special_region) if self.config.special_region is not None else None

        # Load ignore region
        self.ignore_region = SkyRegionList.from_file(self.config.ignore_region) if self.config.ignore_region is not None else None

        # Ignore stars in certain images
        if "ignore_stars" in kwargs: self.ignore_stars = kwargs.pop("ignore_stars")

        # Catalogs
        self.load_catalogs(**kwargs)

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
            #if image.has_errors: self.add_error_map(name, image.errors)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the frames
            self.frames = dataset.get_frames()

            # Get the error maps
            #self.error_maps = dataset.get_errormaps()

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

                # Create the mask
                special_mask = Mask.from_region(self.special_region, self.frames[name].xsize, self.frames[name].ysize)

                self.special_masks[name] = special_mask

        # Ignore mask
        if self.ignore_region is not None:

            # Loop over the frames
            for name in self.frames:

                # Create the mask
                ignore_mask = Mask.from_region(self.ignore_region, self.frames[name].xsize, self.frames[name].ysize)

                self.ignore_masks[name] = ignore_mask

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
        THis function ...
        :param kwargs:
        :return: 
        """

        # Inform the user
        log.info("Loading catalog of extended sources ...")

        # From kwargs
        if "extended_source_catalog" in kwargs and kwargs["extended_source_catalog"] is not None: self.extended_source_catalog = kwargs.pop("extended_source_catalog")

        # From file
        elif self.config.extended_sources_catalog is not None: self.extended_source_catalog = ExtendedSourceCatalog.from_file(self.config.extended_sources_catalog)

        else: raise ValueError("Catalog of extended sources has to be specified")

    # -----------------------------------------------------------------

    def load_point_sources_catalog(self, kwargs):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Loading catalog of point sources ...")

        # From kwargs
        if "point_source_catalog" in kwargs and kwargs["point_source_catalog"] is not None: self.point_source_catalog = kwargs.pop("point_source_catalog")

        # From file
        elif self.config.point_sources_catalog is not None: self.point_source_catalog = PointSourceCatalog.from_file(self.config.point_sources_catalog)

        # Not specified
        else: raise ValueError("Catalog of point sources has to be specified")

    # -----------------------------------------------------------------

    def create_extended_regions(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating regions of extended regions ...")

        self.extended_regions = dict()

        default_radius = parse_quantity("20 arcsec")

        # Create sky region
        regions = self.extended_source_catalog.create_regions(default_radius=default_radius, add_point=True)

        #print(regions)

        # Check for which images the source finding step has already been performed
        for name in self.frames:

            if self.ignore is not None and name in self.ignore: continue

            # Debugging
            log.debug("Creating extended source regions in pixel coordinates for the " + name + " image ...")

            # Get wcs
            wcs = self.frames[name].wcs

            # Convert to pixel
            pixel_regions = regions.to_pixel(wcs)

            # Add
            self.extended_regions[name] = pixel_regions

            ####
            # Temporary: write here

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "galaxies_original.reg")
            else: path = self.output_path_file("galaxies_original_" + name + ".reg")

            # Save
            pixel_regions.saveto(path)

    # -----------------------------------------------------------------

    def create_point_regions(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Creating regions of point sources ...")

        # Initialize dictionary for the point source regions
        self.point_regions = dict()

        # Create sky region (circles in blue)
        regions = self.point_source_catalog.create_regions(add_point=True, color="blue")

        for name in self.frames:

            # Ignore if requested
            #if name in self.ignore: continue
            if name in self.ignore_stars: continue
            if self.ignore is not None and name in self.ignore: continue

            # Get the frame
            frame = self.frames[name]

            # Don't run the star finder if the wavelength of this image is greater than 25 micron
            if frame.wavelength is None or frame.wavelength > wavelengths.ranges.ir.mir.max:

                # No star subtraction for this image
                log.info("Marking point sources will not be performed for the '" + name + "' image")
                continue

            # Debugging
            log.debug("Creating point source regions in pixel coordinates for the " + name + " image ...")

            # Get wcs
            wcs = frame.wcs

            # Get the filter
            fltr = parse_filter(name)

            # Get the FWHm
            if fltr in self.fwhms: fwhm = self.fwhms[fltr]
            elif self.frames[name].fwhm is not None: fwhm = self.frames[name].fwhm
            else: fwhm = self.config.default_fwhm

            #default_radius = 20.0
            sigma_level = self.config.sigma_level

            # Calculate the radius in pixels
            radius = fwhm * statistics.fwhm_to_sigma * sigma_level
            radius = radius.to("arcsec").value / wcs.average_pixelscale.to("arcsec").value

            # Convert to pixel
            pixel_regions = regions.to_pixel(wcs)

            # Change radii
            for region in pixel_regions:
                if isinstance(region, PixelCircleRegion): region.radius = radius

            # Add to dictionary
            self.point_regions[name] = pixel_regions

            ### TEmporary: write here

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "stars_original.reg")
            else: path = self.output_path_file("stars_original_" + name + ".reg")

            # Save
            pixel_regions.saveto(path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write regions
        self.write_regions()

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

    # -----------------------------------------------------------------

    def write_galaxy_regions(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the galaxy regions ...")

        # Loop over the images
        for name in self.extended_regions:

            # Get regions
            regions = self.extended_regions[name]

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "galaxies.reg")
            else: path = self.output_path_file("galaxies_" + name + ".reg")

            # Save
            regions.saveto(path)

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

            # Get regions
            regions = self.point_regions[name]

            # Determine path
            if self.output_paths is not None and name in self.output_paths: path = fs.join(self.output_paths[name], "stars.reg")
            else: path = self.output_path_file("stars_" + name + ".reg")

            # Save
            regions.saveto(path)

# -----------------------------------------------------------------
