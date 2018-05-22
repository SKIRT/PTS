#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.initialization Contains the PreparationInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...magic.sources.finder import SourceFinder
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.misc.imageimporter import ImageImporter
from ...magic.core.frame import Frame
from ...magic.core.dataset import DataSet
from ...core.launch.pts import PTSRemoteLauncher
from ...core.filter.filter import parse_filter
from ...magic.convolution.kernels import get_fwhm, has_variable_fwhm
from ...magic.catalog.extended import ExtendedSourceCatalog
from ...magic.catalog.point import PointSourceCatalog
from ...magic.catalog.fetcher import CatalogFetcher
from ...dustpedia.core.properties import DustPediaProperties
from ...magic.tools import statistics
from ...magic.sources.marker import SourceMarker
from ...core.units.parsing import parse_unit as u
from ...core.remote.remote import Remote
from ..core.steps import cached_directory_name_for_single_command
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class PreparationInitializer(PreparationComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(PreparationInitializer, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frame paths
        self.paths = dict()
        self.error_paths = dict()

        # The paths to the initialized files
        self.initialized_paths = dict()

        # The source finder
        self.finder = None

        # The initial dataset
        self.set = None

        # Create the PTS remote launcher
        self.launcher = PTSRemoteLauncher()

        # The FWHMs found by the source finder
        self.fwhms = None

        # Catalogs
        self.extended_sources = None
        self.point_sources = None

        # Sources output directories
        self.sources_output_paths = dict()

        # Ignore source finding
        self.ignore_images = []

        # The caching remote
        self.remote = None

        # The remote caching directory
        self.remote_data_path = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the image paths
        self.get_paths()

        # 3. Create the initialized images
        self.initialize_images()

        # Create directories
        self.create_directories()

        # Get the catalogs
        if self.needs_catalogs: self.get_catalogs()

        # 5. Find sources
        if self.needs_sources: self.get_sources()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PreparationInitializer, self).setup(**kwargs)

        # Load the dataset
        if fs.is_file(self.initial_dataset_path): self.set = DataSet.from_file(self.initial_dataset_path)
        else: self.set = DataSet() # create new dataset

        # Set dataset path (so that we can do 'save')
        self.set.path = self.initial_dataset_path

        # Setup the remote
        self.remote = Remote(host_id=self.environment.cache_host_id)

        # Create the cache directory
        directory_name = cached_directory_name_for_single_command(self.environment, self.command_name())
        self.remote_data_path = fs.join(self.remote.home_directory, directory_name)
        if self.config.cache:
            if not self.remote.is_directory(self.remote_data_path): self.remote.create_directory(self.remote_data_path)

        # Set the path for the source finder to the preparation path
        self.config.sources.path = self.prep_path

        # Setup the remote PTS launcher
        if self.config.remote is not None: self.launcher.setup(self.config.remote)
        else: self.finder = SourceFinder(self.config.sources) # Create the source finder

    # -----------------------------------------------------------------

    def get_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for images and error frames ...")

        # Get the paths
        self.paths, self.error_paths = self.get_data_image_and_error_paths()

    # -----------------------------------------------------------------

    def initialize_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the images ...")

        # Loop over all image paths
        for prep_name in self.paths:

            # Get the image path
            image_path = self.paths[prep_name]

            # Determine error path
            error_path = self.error_paths[prep_name] if prep_name in self.error_paths else None

            # Determine the output path for this image
            output_path = self.get_prep_path(prep_name)

            # Set initialized path
            initialized_path = fs.join(output_path, "initialized.fits")
            self.initialized_paths[prep_name] = initialized_path

            # Check whether this image already has an initialized image
            if fs.is_file(initialized_path):

                # Already present
                log.success("Initialized '" + prep_name + "' is already present")

                # Check wether the sources directory is present
                # -> NO, we are going to loop over self.paths again anyway in the create_directories and get_sources functions ...

                # Check whether in dataset
                if prep_name not in self.set.names:

                    # Give a warning
                    log.warning("Initialized '" + prep_name + "' was not yet in the dataset: adding it now ...")
                    # Add to the dataset
                    # Add entry to the dataset
                    self.set.add_path(prep_name, initialized_path)
                    self.set.save() # Save

                # Cache
                self.cache_image(prep_name, image_path, error_path)

                # Now skip the rest
                continue

            # Debugging
            log.debug("Initializing image '" + image_path + "' ...")

            # Set the path to the region of bad pixels
            bad_region_path = fs.join(self.data_path, "bad", prep_name + ".reg")
            if not fs.is_file(bad_region_path): bad_region_path = None

            # Get the filter
            fltr = parse_filter(prep_name)

            # Set the FWHM if the instrument has a fixed PSF
            if has_variable_fwhm(fltr): fwhm = None
            else: fwhm = get_fwhm(fltr)

            # Debugging
            log.debug("Loading image " + image_path + " as " + prep_name + " ...")

            # Import the image
            importer = ImageImporter()
            importer.run(image_path, bad_region_path, fwhm=fwhm, find_error_frame=False) # don't look for error frames

            # Get the imported image
            image = importer.image

            # Set the image name
            image.name = prep_name

            # Remove all frames except for the primary frame
            image.remove_frames_except("primary")

            # If a poisson error map was found, add it to the image
            if prep_name in self.error_paths:

                # Debugging
                log.debug("Adding the poisson error frame to the " + prep_name + " image ...")

                # Add the error frame
                error_map = Frame.from_file(self.error_paths[prep_name])
                image.add_frame(error_map, "errors")

            # Save the image
            image.saveto(initialized_path)

            # Success
            log.success("Initialized the '" + prep_name + "' image")

            # NEW: ADD TO THE DATASET NOW AND SAVE IT IMMEDIATELY

            # Inform the user
            log.info("Adding the '" + prep_name + "' image to the dataset ...")

            # Add entry to the dataset
            self.set.add_path(prep_name, initialized_path)

            # Cache the original image
            self.cache_image(prep_name, image_path, error_path)

    # -----------------------------------------------------------------

    def cache_image(self, prep_name, image_path, error_path=None):

        """
        This function ...
        :param prep_name:
        :param image_path:
        :param error_path:
        :return:
        """

        # Inform the user
        log.info("Caching the original '" + prep_name + "' image from '" + image_path + "' to remote host '" + self.remote.host_id + "' ...")

        # Determine filename
        filename = fs.strip_extension(fs.name(image_path))

        # Determine directory name
        dirname = fs.name(fs.directory_of(image_path))

        # Debugging
        log.debug("Caching " + filename + " image ...")

        # Determine the remote directory for this image
        remote_directory_path = fs.join(self.remote_data_path, dirname)
        if not self.remote.is_directory(remote_directory_path): self.remote.create_directory(remote_directory_path)

        # Debugging
        log.debug("Uploading the file to '" + remote_directory_path + "' ...")

        # Upload
        self.remote.upload(image_path, remote_directory_path)

        # Debugging
        log.debug("Removing the local file (" + image_path + ") ...")

        # Remove the file
        fs.remove_file(image_path)

        # Cache the error map
        if error_path is not None:

            # Debugging
            log.debug("Caching the " + prep_name + " error map ...")

            # Upload
            self.remote.upload(error_path, remote_directory_path)

            # Debugging
            log.debug("Removing the local file (" + error_path + ") ...")

            # Remove the file
            fs.remove_file(error_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def catalog_coordinate_box(self):

        """
        This function ...
        :return: 
        """

        # Determine the bounding box
        if self.config.catalog_overlapping: coordinate_box = self.set.get_overlap_box()
        else: coordinate_box = self.set.get_bounding_box()

        # Return
        return coordinate_box

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating directories ...")

        # Loop
        for name in self.paths:

            # Get output path
            output_path = self.get_prep_path(name)

            # Determine the path to the "sources" directory within the output path for this image
            sources_output_path = fs.join(output_path, "sources")
            if not fs.is_directory(sources_output_path): fs.create_directory(sources_output_path)

            # If the source finding step has already been performed on this image, don't do it again
            elif not fs.is_empty(sources_output_path):

                # Debugging
                log.debug("Source finder output has been found for the " + name + " image, skipping source finding step ...")

                # Ignore this image for the source finder
                self.ignore_images.append(name)

            # Set path
            self.sources_output_paths[name] = sources_output_path

    # -----------------------------------------------------------------

    def get_catalogs(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the catalogs ...")

        # Extended
        self.get_extended_sources_catalog()

        # Point
        self.get_point_sources_catalog()

    # -----------------------------------------------------------------

    def get_extended_sources_catalog(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the catalog of extended sources ...")

        # Search for catalogs that are saved in the prepration directory
        # extended_source_catalog and point_source_catalog
        extended_sources_path = fs.join(self.prep_path, "extended_sources.dat")

        # Load or fetch the catalog
        if fs.is_file(extended_sources_path): self.extended_sources = ExtendedSourceCatalog.from_file(extended_sources_path)
        else: self.fetch_extended_sources_catalog(extended_sources_path)

    # -----------------------------------------------------------------

    def fetch_extended_sources_catalog(self, path):

        """
        This function ...
        :param path:
        :return: 
        """

        # Inform the user
        log.info("Fetching the catalog of extended sources ...")

        # Fetch
        fetcher = CatalogFetcher()
        self.extended_sources = fetcher.get_extended_source_catalog(self.catalog_coordinate_box)

        # Save the catalog
        self.extended_sources.saveto(path)

    # -----------------------------------------------------------------

    def get_point_sources_catalog(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the catalog of point sources ...")

        # Search for catalogs that are saved in the prepration directory
        # Set the path
        point_sources_path = fs.join(self.prep_path, "point_sources.dat")

        # Load or fetch the catalog
        if fs.is_file(point_sources_path): self.point_sources = PointSourceCatalog.from_file(point_sources_path)
        else: self.fetch_point_sources_catalog(point_sources_path)

    # -----------------------------------------------------------------

    def fetch_point_sources_catalog(self, path):

        """
        This function ...
        :param path
        :return: 
        """

        # Inform the user
        log.info("Fetching the catalog of point sources ...")

        # Fetch
        fetcher = CatalogFetcher()
        min_pixelscale = self.set.min_pixelscale
        self.point_sources = fetcher.get_point_source_catalog(self.catalog_coordinate_box, min_pixelscale, self.config.catalogs)

        # Save the catalog
        self.point_sources.saveto(path)

    # -----------------------------------------------------------------

    @property
    def nignored_images(self):

        """
        This function ...
        :return: 
        """

        return len(self.ignore_images)

    # -----------------------------------------------------------------

    @property
    def has_non_ignored_images(self):

        """
        This function ...
        :return: 
        """

        return len(self.set) - self.nignored_images > 0

    # -----------------------------------------------------------------

    @property
    def needs_catalogs(self):

        """
        This function ...
        :return: 
        """

        return self.has_non_ignored_images

    # -----------------------------------------------------------------

    @property
    def needs_sources(self):

        """
        This function ...
        :return: 
        """

        return self.has_non_ignored_images

    # -----------------------------------------------------------------

    def get_sources(self):

        """
        This fucntion ...
        :return: 
        """

        # Mark or find
        if self.config.manual: self.mark_sources()
        else: self.find_sources()

    # -----------------------------------------------------------------

    def mark_sources(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Marking sources in the images ...")

        # Get the DustPedia properties instance
        properties = DustPediaProperties()

        # Get FWHMs
        fwhms = properties.fwhms

        # Create the source marker
        marker = SourceMarker()

        # Set the path for the source finder to the preparation path
        marker.config.path = self.prep_path

        # Don't look for stars in the Halpha image
        # Look for name of Halpha image
        ignore_stars = []
        for name in self.paths:
            fltr = parse_filter(name)
            if fltr == "Ha": ignore_stars.append(name)
        #ignore_stars = [parse_filter("Halpha")]

        default_fwhm = 2.0 * u("arcsec")

        marker.config.default_fwhm = default_fwhm

        # Run
        marker.run(fwhms=fwhms, dataset=self.set, ignore=self.ignore_images, extended_source_catalog=self.extended_sources,
                   point_source_catalog=self.point_sources, ignore_stars=ignore_stars, output_paths=self.sources_output_paths)

        # Set the FWHMs
        self.fwhms = dict()
        for name in self.paths:
            fltr = parse_filter(name)
            if fltr in fwhms: fwhm = fwhms[fltr]
            else:
                frame = Frame.from_file(self.paths[name])
                if frame.fwhm is not None: fwhm = frame.fwhm
                else: fwhm = default_fwhm
            self.fwhms[name] = fwhm

        # Set the FWHMs
        self.set_fwhms()

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources in the images ...")

        # Don't look for stars in the Halpha image
        ignore_stars = ["Mosaic Halpha", "Halpha", "Ha"]

        # Don't look for other sources in the IRAC images
        ignore_other_sources = ["IRAC I1", "IRAC I2", "IRAC I3", "IRAC I4"]

        # Create an animation for the source finder
        #if self.config.visualise: animation = Animation()
        #else: animation = None

        # Find sources locally or remotely
        if self.config.remote is not None: self.find_sources_remote(self.ignore_images, ignore_stars, ignore_other_sources)
        else: self.find_sources_local(self.ignore_images, ignore_stars, ignore_other_sources)

        # Set FWHM of optical images
        self.set_fwhms()

    # -----------------------------------------------------------------

    def find_sources_local(self, ignore_images, ignore_stars, ignore_other_sources):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources locally ...")

        # Run the source finder
        self.finder.run(dataset=self.set, ignore=ignore_images, ignore_stars=ignore_stars,
                        ignore_other_sources=ignore_other_sources, extended_source_catalog=self.extended_sources,
                        point_source_catalog=self.point_sources, output_paths=self.sources_output_paths)

        # Get the statistics
        #self.statistics = self.finder.statistics

        # Get the fwhms
        self.fwhms = self.finder.fwhms

    # -----------------------------------------------------------------

    def find_sources_remote(self, ignore_images, ignore_stars, ignore_other_sources):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding sources remotely on host '" + self.config.remote + "'...")

        # Initialize the input dictionary
        input_dict = dict()

        # Set the input dataset
        input_dict["dataset"] = self.set
        input_dict["ignore_images"] = ignore_images
        input_dict["ignore_stars"] = ignore_stars
        input_dict["ignore_other_sources"] = ignore_other_sources
        input_dict["extended_source_catalog"] = self.extended_sources
        input_dict["point_source_catalog"] = self.point_sources
        input_dict["output_paths"] = self.sources_output_paths

        # Run the PTS find_sources command remotely and get the output
        #self.statistics = self.launcher.run_attached("find_sources", self.config.sources, input_dict, return_output_names=["statistics"], unpack=True)
        self.fwhms = self.launcher.run_attached("find_sources", self.config.sources, input_dict, return_output_names=["fwhms"], unpack=True)

    # -----------------------------------------------------------------

    def set_fwhms(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the FWHMs of the images ...")

        # Loop over the images
        for prep_name in self.set:

            if self.ignore_images is not None and prep_name in self.ignore_images: continue

            # Get the filter
            #fltr = parse_filter(prep_name)

            # Set the FWHM if the instrument has a fixed PSF
            #if has_variable_fwhm(fltr):

            # Open the image, set the FWHM, and save again
            image = self.set.get_image(prep_name)
            image.fwhm = self.fwhms[prep_name]
            image.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the dataset
        self.write_dataset()

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dataset ...")

        # Save the dataset
        self.set.save()

        # If already present
        #if fs.is_file(self.initial_dataset_path): fs.remove_file(self.initial_dataset_path)

        # Save the dataset
        #self.set.saveto(self.initial_dataset_path)

# -----------------------------------------------------------------
