#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.train.collection Contains the Collector class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

# Import the relevant PTS classes and modules
from pts.magic.core import Image, Source
from ...core.basics.configurable import Configurable
from ...core.tools import filesystem, inspection

# -----------------------------------------------------------------

class Collector(Configurable):

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
        super(Collector, self).__init__(config, "magic")

        # The current source
        self.current_source = None

        # Determine the path to the magic/collection user directory
        self.collection_user_path = os.path.join(inspection.pts_user_dir, "magic", "collection")

        # Create the user collection directory
        filesystem.create_directory(self.collection_user_path)

        self.yes_path = None
        self.no_path = None

        self.current_index_yes = -1
        self.current_index_no = -1

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Collector instance
        collector = cls()

        # -- Adjust the configuration settings according to the command-line arguments --

        collector.mode = arguments.mode

        # Return the collector
        return collector

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Collect
        self.collect()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        self.current_source = None
        self.yes_path = None
        self.no_path = None
        self.current_index_yes = -1
        self.current_index_no = -1

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Collector, self).setup()

        # Determine the paths to the star collection and saturation collection directories
        mode_path = os.path.join(self.collection_user_path, self.config.mode)

        # Create the star collection and saturation collection directories
        filesystem.create_directory(mode_path)

        # Determine the paths to the 'yes' and 'no' saturation collection directories
        self.yes_path = os.path.join(mode_path, "yes")
        self.no_path = os.path.join(mode_path, "no")

        # Current index of saturation sources
        self.current_index_yes = -1
        self.current_index_no = -1

        for path in filesystem.files_in_path(self.yes_path, extension="fits"):

            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_yes: self.current_index_yes = index

        for path in filesystem.files_in_path(self.no_path, extension="fits"):

            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_no: self.current_index_no = index

    # -----------------------------------------------------------------

    def collect(self):

        """
        This function ...
        :return:
        """

        # Loop over all FITS files found in the current directory
        for file_path in filesystem.files_in_path(os.getcwd(), extension="fits", contains=self.config.mode):

            # Get information
            #name = os.path.basename(file_path).split(".fits")[0]
            #object_type, level, counter = name.split("_")

            # Open the image, select all frames
            image = Image(file_path, always_call_first_primary=False)
            image.frames.select_all()

            # Create a source
            source = Source.from_image(image)

            # Create a plot for the source
            source.plot(title="Is this a saturated star?", show=False, scale="log")

            # Add 'yes' and 'no' buttons
            self.current_source = source
            axyes = plt.axes([0.7, 0.05, 0.1, 0.075])
            axno = plt.axes([0.81, 0.05, 0.1, 0.075])
            yes_button = Button(axyes, 'Yes')
            yes_button.on_clicked(self.save_yes)
            no_button = Button(axno, 'No')
            no_button.on_clicked(self.save_no)

            # Show the plot
            plt.show()

    # -----------------------------------------------------------------

    def save_yes(self, event):

        """
        This function ...
        :param event:
        :return:
        """

        # Increment index
        self.current_index_yes += 1

        # Determine the path to the new FITS file
        path = os.path.join(self.yes_path, str(self.current_index_yes) + ".fits")

        # Inform the user and save the source object
        self.log.info("Saving the saturated star to " + path)
        self.current_source.save(path)

        # Close the currently active plotting window for this source
        plt.close()

    # -----------------------------------------------------------------

    def save_no(self, event):

        """
        This function ...
        :param event:
        :return:
        """

        # Increment index
        self.current_index_no += 1

        # Determine the path to the new FITS file
        path = os.path.join(self.no_path, str(self.current_index_no) + ".fits")

        # Inform the user and save the source object
        self.log.info("Saving the star without saturation to " + path)
        self.current_source.save(path)

        # Close the currently active plotting window for this source
        plt.close()

# -----------------------------------------------------------------
