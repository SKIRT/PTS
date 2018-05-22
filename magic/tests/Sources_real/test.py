#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.magic.services.ned import get_image

# -----------------------------------------------------------------

description = "Test the source detection and extraction on a real image"

# -----------------------------------------------------------------

class RealSourcesTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RealSourcesTest, self).__init__(*args, **kwargs)

        # The image frame
        self.frame = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the image
        self.get_image()

        # 13. Find sources
        if self.config.manual: self.mark_sources()
        else: self.find_sources()

        # 14. Extract sources
        self.extract()

        # 15. Write
        self.write()

        # 16. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RealSourcesTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def get_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the image ...")

        # Get the image
        self.frame = get_image(self.config.galaxy_name, self.config.filter, self.config.year)

    # -----------------------------------------------------------------

    def fetch_catalogs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching catalogs ...")

        # Get extended source catalog
        self.extended_source_catalog = self.fetcher.get_extended_source_catalog(self.coordinate_systems.bounding_box)

        # Fetch
        self.point_source_catalog = self.fetcher.get_point_source_catalog(self.coordinate_systems.bounding_box, self.coordinate_systems.min_pixelscale, self.config.point_source_catalogs)

    # -----------------------------------------------------------------

    def mark_sources(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def find_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding the sources ...")

        # Settings
        settings = dict()
        #settings["input"] =
        settings["output"] = self.find_path
        settings["nprocesses"] = self.config.nprocesses

        # Input
        input_dict = dict()
        input_dict["dataset"] = self.dataset
        input_dict["extended_source_catalog"] = self.extended_source_catalog
        input_dict["point_source_catalog"] = self.point_source_catalog
        input_dict["output_paths"] = self.find_paths

        # Construct the command
        command = Command("find_sources", "find sources", settings, input_dict)

        # Run the command
        self.finder = self.run_command(command, remote=self.remote)

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting the sources ...")

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
