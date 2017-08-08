#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod
import math
import numpy as np

# Import astronomical modules
from astropy.modeling.models import Gaussian2D, AiryDisk2D
from photutils.datasets import make_noise_image

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.core.frame import Frame
from pts.core.tools import filesystem as fs
from pts.magic.core.dataset import DataSet
from pts.core.units.parsing import parse_unit as u
from pts.modeling.tests.base import m81_data_path
from pts.magic.catalog.fetcher import CatalogFetcher
from pts.magic.basics.coordinate import SkyCoordinate
from pts.core.tools import stringify
from pts.magic.tools import wavelengths
from pts.magic.tools import fitting, statistics
from pts.magic.convolution.kernels import has_variable_fwhm, get_fwhm
from pts.magic.basics.vector import Pixel
from pts.magic.basics.coordinate import PixelCoordinate
from pts.magic.core.list import CoordinateSystemList

# -----------------------------------------------------------------

description = "Test the source detection and extraction"

# -----------------------------------------------------------------

# For M81
fwhms = {"2MASS H": 4.640929858306589 * u("arcsec"),
          "2MASS J": 4.580828087551186 * u("arcsec"),
          "2MASS Ks": 4.662813601376219 * u("arcsec"),
          "SDSS g": 2.015917936060279 * u("arcsec"),
          "SDSS i": 1.85631074608032 * u("arcsec"),
          "SDSS r": 2.026862297071852 * u("arcsec"),
          "SDSS u": 2.327165667182196 * u("arcsec"),
          "SDSS z": 1.841443699129355 * u("arcsec")}

# -----------------------------------------------------------------

# Determine the path to the headers directory
headers_path = fs.join(m81_data_path, "headers")

# -----------------------------------------------------------------

class SourcesTestBase(TestImplementation):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SourcesTestBase, self).__init__(*args, **kwargs)

        # The remote
        self.remote = None

        # Paths
        self.data_path = None
        self.data_frames_path = None
        self.data_masks_path = None
        self.find_path = None
        self.find_paths = dict()
        self.extract_path = None
        self.extract_paths = dict()

        # The coordinate systems
        self.coordinate_systems = CoordinateSystemList()

        # The frames
        self.frames = dict()
        
        # The dataset
        self.dataset = None

        # The catalog fetcher
        self.fetcher = CatalogFetcher()

        # Catalogs
        self.point_source_catalog = None

        # The real FWHMs
        self.real_fwhms = dict()

        # The source finder
        self.finder = None

        # The source extractors
        self.extractors = dict()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SourcesTestBase, self).setup(**kwargs)

        # Set paths
        self.data_path = fs.create_directory_in(self.path, "data")
        self.data_frames_path = fs.create_directory_in(self.data_path, "frames")
        self.data_masks_path = fs.create_directory_in(self.data_path, "masks")
        self.find_path = fs.create_directory_in(self.path, "find")
        self.extract_path = fs.create_directory_in(self.path, "extract")

    # -----------------------------------------------------------------

    def initialize_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the frames ...")

        # Loop over the filters
        for fltr in self.coordinate_systems.filters:

            # Debugging
            log.debug("Initializing the '" + str(fltr) + "' frame ...")

            # Get the wcs
            wcs = self.coordinate_systems[fltr]

            # Create new frame
            frame = Frame.zeros(wcs.shape)

            # Add the wcs
            frame.wcs = wcs

            # Set the filter
            frame.filter = fltr

            # Set the unit
            frame.unit = "Jy"

            # Add the frame
            self.frames[fltr] = frame

    # -----------------------------------------------------------------

    def set_fwhms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the FWHMs ...")

        # Loop over the filters
        for fltr in self.frames:

            # Get the fwhm
            if has_variable_fwhm(fltr): fwhm = fwhms[str(fltr)]
            else: fwhm = get_fwhm(fltr)

            # Debugging
            log.debug("The FWHM of the '" + str(fltr) + "' image is " + stringify.stringify(fwhm)[1])

            # Set
            self.real_fwhms[fltr] = fwhm
            
    # -----------------------------------------------------------------

    @property
    def star_filters(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the filters
        for fltr in self.frames:

            # Get wavelength
            wavelength = fltr.effective if fltr.effective is not None else fltr.center

            # Check
            if wavelength > wavelengths.ranges.ir.mir.max: continue
            filters.append(fltr)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    @property
    def extra_filters(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the filters
        for fltr in self.frames:

            # Get wavelength
            wavelength = fltr.effective if fltr.effective is not None else fltr.center

            # Check
            if wavelength < wavelengths.ranges.ir.mir.max: continue
            filters.append(fltr)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    def create_random_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating random point sources ...")

        # Generate random coordinates
        right_ascensions = np.random.uniform(self.coordinate_systems.min_ra_deg, self.coordinate_systems.max_ra_deg, size=self.config.nrandom_sources)
        declinations = np.random.uniform(self.coordinate_systems.min_dec_deg, self.coordinate_systems.max_dec_deg, size=self.config.nrandom_sources)

        # Loop over the coordinates
        for ra, dec in zip(right_ascensions, declinations):

            # Create a sky coordinate
            coordinate = SkyCoordinate(ra=ra, dec=dec, unit="deg")

            # Add to the point source catalog
            self.point_source_catalog.add_coordinate(coordinate)

    # -----------------------------------------------------------------

    def make_point_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making point sources ...")

        # Call the appropriate function
        if self.config.vary_fwhm: self.make_point_sources_variable_fwhm()
        else: self.make_point_sources_fixed_fwhm()

    # -----------------------------------------------------------------

    def make_point_sources_variable_fwhm(self, masks=None):

        """
        THis function ...
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Making point sources with a variable FWHM in each band ...")

        # Loop over the 'star' filters
        for fltr in self.star_filters:

            # Debugging
            log.debug("Making point sources for the '" + str(fltr) + "' image ...")

            # Get y and x
            if not self.config.only_local:
                y, x = np.indices(self.frames[fltr].shape)
            else: y = x = None

            # Get pixelscale
            pixelscale = self.frames[fltr].average_pixelscale

            # Debugging
            log.debug("The pixelscale of the image is " + stringify.stringify(pixelscale)[1])

            # Determine the FWHM in pixel coordinates
            fwhm_pix = self.real_fwhms[fltr].to("arcsec").value / pixelscale.to("arcsec").value

            counter = 0

            # Loop over the coordinates in the point sources catalog
            for coordinate in self.point_source_catalog.coordinates():

                counter += 1

                # Debugging
                log.debug("Adding point source " + str(counter) + " of " + str(len(self.point_source_catalog)) + " ...")

                # Check whether it falls in the frame
                if not self.frames[fltr].contains(coordinate): continue

                # Convert into pixel coordinate
                pixel_coordinate = coordinate.to_pixel(self.frames[fltr].wcs)

                # Get the corresponding pixel
                pixel = Pixel.for_coordinate(pixel_coordinate)

                # Check whether not masked
                if masks is not None and fltr in masks and masks[fltr][pixel.y, pixel.x]: continue

                # Generate random deviation
                x_deviation = np.random.normal(0.0, 1.)
                y_deviation = np.random.normal(0.0, 1.)

                # Debugging
                log.debug("Random pixel position deviation is (" + str(x_deviation) + ", " + str(y_deviation) + ")")

                # Alter pixel coordinate
                pixel_coordinate.x += x_deviation
                pixel_coordinate.y += y_deviation

                # Generate random deviation from FWHM
                fwhm_deviation = np.random.normal(0.0, 0.05 * fwhm_pix)

                # Debugging
                log.debug("Random FWHM deviation (on a FWHM of " + str(fwhm_pix) + ") is " + str(fwhm_deviation))

                # Add the deviation
                fwhm_pix += fwhm_deviation

                # Generate a random amplitude (from 100 till 100 000)
                amplitude_exponent = np.random.uniform(2., 5.)
                amplitude = 10**amplitude_exponent

                # Determine sigma
                sigma = statistics.fwhm_to_sigma * fwhm_pix

                # Only 'render' the Gaussian locally
                if self.config.only_local: min_x, max_x, min_y, max_y, rel_center = determine_patch_for_psf(pixel_coordinate, amplitude, sigma, max_x_pixel=self.frames[fltr].xsize-1, max_y_pixel=self.frames[fltr].ysize-1, keep_in_frame=True)

                # Render each Gaussian over the entire frame
                else:
                    min_x = max_x = min_y = max_y = None
                    rel_center = pixel_coordinate

                # Make the model
                model = create_model(self.config.psf_model, amplitude, rel_center, sigma)

                # Evaluate the model
                if self.config.only_local: y, x = np.indices((max_y - min_y, max_x - min_x))
                data = model(x, y)

                # Add the data
                if self.config.only_local: self.frames[fltr][min_y:max_y, min_x:max_x] += data
                else: self.frames[fltr] += data

    # -----------------------------------------------------------------

    def make_point_sources_fixed_fwhm(self, masks=None):

        """
        This function ...
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Making point sources with a fixed FWHM in each band ...")

        # Loop over the 'star' filters
        for fltr in self.star_filters:

            # Debugging
            log.debug("Making point sources for the '" + str(fltr) + "' image ...")

            # Get pixelscale
            pixelscale = self.frames[fltr].average_pixelscale

            # Debugging
            log.debug("The pixelscale of the image is " + stringify.stringify(pixelscale)[1])

            # Determine the FWHM in pixel coordinates
            fwhm_pix = self.real_fwhms[fltr].to("arcsec").value / pixelscale.to("arcsec").value

            # Determine sigma
            sigma = statistics.fwhm_to_sigma * fwhm_pix

            # Only 'render' the Gaussian locally
            origin = PixelCoordinate(0.0, 0.0)
            amplitude = 1.
            min_x, max_x, min_y, max_y, rel_center = determine_patch_for_psf(origin, amplitude, sigma, keep_in_frame=False)

            # Make the model
            model = create_model(self.config.psf_model, amplitude, rel_center, sigma)

            # Evalute the model
            y, x = np.indices((max_y - min_y, max_x - min_x))
            data = model(x, y)

            # Keep track of the number of sources
            counter = 0

            # Loop over the coordinates in the point sources catalog
            for coordinate in self.point_source_catalog.coordinates():

                counter += 1

                # Debugging
                log.debug("Adding source " + str(counter) + " of " + str(len(self.point_source_catalog)) + " ...")

                # Check whether it falls in the frame
                if not self.frames[fltr].contains(coordinate): continue

                # Convert into pixel coordinate
                pixel_coordinate = coordinate.to_pixel(self.frames[fltr].wcs)

                # Get the corresponding pixel
                pixel = Pixel.for_coordinate(pixel_coordinate)

                # Check whether not masked
                if masks is not None and fltr in masks and masks[fltr][pixel.y, pixel.x]: continue

                # Generate a random amplitude (from 100 till 100 000)
                amplitude_exponent = np.random.uniform(2., 5.)
                amplitude = 10 ** amplitude_exponent

                # Calculate absolute minima and maxima
                source_min_x = min_x + pixel.x
                source_max_x = max_x + pixel.x
                source_min_y = min_y + pixel.y
                source_max_y = max_y + pixel.y

                source_xsize = source_max_x - source_min_x
                source_ysize = source_max_y - source_min_y

                # Correct
                if source_min_x < 0:
                    cut_x_min = - source_min_x
                    source_min_x = 0
                else: cut_x_min = 0

                if source_max_x >= self.frames[fltr].xsize:
                    cut_x_max = source_xsize - (source_max_x - self.frames[fltr].xsize)
                    source_max_x = self.frames[fltr].xsize
                else: cut_x_max = source_xsize

                if source_min_y < 0:
                    cut_y_min = - source_min_y
                    source_min_y = 0
                else: cut_y_min = 0

                if source_max_y >= self.frames[fltr].ysize:
                    cut_y_max = source_ysize - (source_max_y - self.frames[fltr].ysize)
                    source_max_y = self.frames[fltr].ysize
                else: cut_y_max = source_ysize

                # Create the final data for this source
                source_data = data[cut_y_min:cut_y_max, cut_x_min:cut_x_max] * amplitude

                # Add the data
                self.frames[fltr][source_min_y:source_max_y, source_min_x:source_max_x] += source_data

    # -----------------------------------------------------------------

    def make_noise(self, masks=None):

        """
        This function ...
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Adding noise ...")

        # Loop over the frames
        for fltr in self.frames:

            # Get the frame
            frame = self.frames[fltr]

            # Make noise and add it
            data = make_noise_image(frame.shape, type='gaussian', mean=0., stddev=self.config.noise_stddev)
            frame += data

            # Mask
            if masks is not None and fltr in masks: frame[masks[fltr]] = 0.0
            
    # -----------------------------------------------------------------

    def create_dataset(self, masks=None):

        """
        This function ...
        :param masks:
        :return:
        """

        # Inform the user
        log.info("Creating the dataset ...")

        # Initialize
        self.dataset = DataSet()

        # Add the frames
        for fltr in self.frames:

            # Determine name for the image
            name = str(fltr)

            # Determine path for this frame
            path = fs.join(self.data_frames_path, name + ".fits")

            # Debugging
            log.debug("Saving the frame ...")

            # Save the frame
            self.frames[fltr].saveto(path)

            # Debugging
            log.debug("Adding the '" + name + "' image to the dataset ...")

            # Add the frame to the dataset
            self.dataset.add_path(name, path)

            # Determine the path for the mask
            mask_path = fs.join(self.data_masks_path, name + ".fits")

            # Mask
            if masks is not None and fltr in masks:

                # Debugging
                log.debug("Saving the mask ...")

                # Save
                masks[fltr].saveto(mask_path)

                # Debugging
                log.debug("Adding mask ...")

                # Add the mask
                self.dataset.add_mask_path(name, mask_path)

        # Determine database path
        path = fs.join(self.path, "database.dat")

        # Write the dataset
        self.dataset.saveto(path)

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the directories ...")

        # Loop over the frames
        for name in self.dataset.names:

            # Find
            self.find_paths[name] = fs.create_directory_in(self.find_path, name)

            # Extract
            self.extract_paths[name] = fs.create_directory_in(self.extract_path, name)

    # -----------------------------------------------------------------

    @abstractmethod
    def find(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the sources ...")

        # Loop over the images
        for fltr in self.frames:

            name = str(fltr)

            # Inform the user
            log.info("Extracting the sources for the '" + name + "' image ...")

            # Settings
            settings = dict()
            settings["input"] = self.find_paths[name]
            settings["output"] = self.extract_paths[name]

            # Input
            input_dict = dict()
            input_dict["frame"] = self.frames[fltr]

            # Construct the command
            command = Command("extract", "extract the sources", settings, input_dict)

            # Run the command
            extractor = self.run_command(command, remote=self.remote)

            # Add the extractor
            self.extractors[fltr] = extractor

# -----------------------------------------------------------------

def determine_patch_for_psf(center, amplitude, sigma, cutoff_fraction=1e-4, max_x_pixel=None, max_y_pixel=None, keep_in_frame=False):

    """
    This function ...
    :param center:
    :param amplitude:
    :param sigma:
    :param cutoff_fraction:
    :param max_x_pixel:
    :param max_y_pixel:
    :param keep_in_frame:
    :return:
    """

    # Check at which sigma level the Gaussian has decreased 3 orders of magnitude
    sigma_level = statistics.inverse_gaussian(0.0, sigma, amplitude, amplitude * cutoff_fraction) / sigma

    # Determine x and y range for evaluation
    min_x = int(math.floor(center.x - sigma_level * sigma))
    if keep_in_frame: min_x = max(min_x, 0)

    max_x = int(math.ceil(center.x + sigma_level * sigma))
    if keep_in_frame: max_x = min(max_x, max_x_pixel)

    min_y = int(math.floor(center.y - sigma_level * sigma))
    if keep_in_frame: min_y = max(min_y, 0)

    max_y = int(math.ceil(center.y + sigma_level * sigma))
    if keep_in_frame: max_y = min(max_y, max_y_pixel)

    # Create relative centers
    rel_x = center.x - min_x
    rel_y = center.y - min_y
    rel_center = PixelCoordinate(rel_x, rel_y)

    # Return
    return min_x, max_x, min_y, max_y, rel_center

# -----------------------------------------------------------------

def create_model(model_name, amplitude, center, sigma):

    """
    This function ...
    :param model_name:
    :param amplitude:
    :param center:
    :param sigma:
    :return:
    """

    # 2D GAussian model
    if model_name == "gaussian":

        # Create the model
        model = Gaussian2D(amplitude=amplitude, x_mean=center.x, y_mean=center.y, x_stddev=sigma, y_stddev=sigma)

    # Airy disk model
    elif model_name == "airydisk":

        # Determine the radius
        radius = fitting.gaussian_sigma_to_airy_radius(sigma)

        # Create the model
        model = AiryDisk2D(amplitude=amplitude, x_0=center.x, y_0=center.y, radius=radius)

    # Invalid
    else: raise ValueError("Not a valid model")

    # Return the model
    return model

# -----------------------------------------------------------------
