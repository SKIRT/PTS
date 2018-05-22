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
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.core.filter.filter import parse_filter
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.magic.basics.coordinate import SkyCoordinate
from pts.magic.basics.vector import Extent
from pts.magic.basics.pixelscale import Pixelscale
from pts.magic.tests.base import SourcesTestBase
from pts.do.commandline import Command
from pts.magic.catalog.point import PointSourceCatalog

# -----------------------------------------------------------------

description = "Test the point source detection and extraction"

# -----------------------------------------------------------------

filter_names = ["GALEX FUV", "GALEX NUV", "SDSS u", "SDSS g", "SDSS r"]

# -----------------------------------------------------------------

class SourcesLightTest(SourcesTestBase):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SourcesLightTest, self).__init__(*args, **kwargs)

        # The center coordinate
        self.center = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set the center coordinate
        self.set_center_coordinate()

        # Create the coordinate system
        self.create_coordinate_systems()

        # Initialize the frames
        self.initialize_frames()

        # Set FWHMs
        self.set_fwhms()

        # Create the catalog of point sources
        self.create_catalog()

        # 3. Generate the sources
        self.make_sources()

        # 4. Make noise
        self.make_noise()

        # Create the dataset
        self.create_dataset()

        # Create directories
        self.create_directories()

        # Find
        self.find()

        # Extract
        self.extract()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SourcesLightTest, self).setup(**kwargs)

        # Set paths
        self.data_path = fs.create_directory_in(self.path, "data")
        self.data_frames_path = fs.create_directory_in(self.data_path, "frames")
        self.data_masks_path = fs.create_directory_in(self.data_path, "masks")
        self.find_path = fs.create_directory_in(self.path, "find")
        self.extract_path = fs.create_directory_in(self.path, "extract")

    # -----------------------------------------------------------------

    def set_center_coordinate(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the center coordinate ...")

        # Generate random ra and dec
        ra = np.random.rand() * 360.0
        dec = np.random.rand() * 180.0 - 90.0

        # Set the center coordinate
        self.center = SkyCoordinate(ra=ra, dec=dec, unit="deg")

    # -----------------------------------------------------------------

    def create_coordinate_systems(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the coordinate systems ...")

        # Loop over the number of frames
        for index in range(self.config.nframes):

            # Get the next filter
            fltr = parse_filter(filter_names[index])

            # Set properties
            size = Extent(self.config.npixels, self.config.npixels)
            center_pixel = size * 0.5
            pixelscale = Pixelscale(1.0, unit="arcsec")

            # Create the coordinate system
            wcs = CoordinateSystem.from_properties(size, center_pixel, self.center, pixelscale)

            # Add the wcs
            self.coordinate_systems.append(wcs, fltr)

    # -----------------------------------------------------------------

    def create_catalog(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the catalog ...")

        # Initialize catalogs
        self.initialize_catalog()

        # Add random sources
        self.create_random_sources()

    # -----------------------------------------------------------------

    def initialize_catalog(self):

        """
        This function ...
        :return: 
        """

        # Point sources
        self.point_source_catalog = PointSourceCatalog()

    # -----------------------------------------------------------------

    def make_sources(self):

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

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding the sources ...")

        # Settings
        settings = dict()
        # settings["input"] =
        settings["output"] = self.find_path
        settings["nprocesses"] = self.config.nprocesses

        # Input
        input_dict = dict()
        input_dict["dataset"] = self.dataset
        #input_dict["extended_source_catalog"] = self.extended_source_catalog
        input_dict["point_source_catalog"] = self.point_source_catalog
        input_dict["output_paths"] = self.find_paths

        # Construct the command
        command = Command("find_sources", "find sources", settings, input_dict)

        # Run the command
        self.finder = self.run_command(command)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------