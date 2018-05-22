#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
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
from ..core.image import Image
from ..core.source import Source
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

description = {"star": "star", "saturation": "saturated star / diffraction pattern"}

# -----------------------------------------------------------------

class Collector(Configurable):

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
        super(Collector, self).__init__(*args, **kwargs)

        # The current and previous source
        self.previous_source = None
        self.current_source = None

        # Determine the path to the magic/collection user directory
        self.collection_user_path = os.path.join(introspection.pts_user_dir, "magic", "collection")

        # Create the user collection directory
        fs.create_directory(self.collection_user_path)

        self.yes_path = None
        self.no_path = None

        self.last_path = None
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

        collector.config.mode = arguments.mode

        # Return the collector
        return collector

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Collector, self).setup(**kwargs)

        # Get description
        self.description = description[self.config.mode]

        # Determine the paths to the star collection and saturation collection directories
        mode_path = os.path.join(self.collection_user_path, self.config.mode)

        # Create the star collection and saturation collection directories
        fs.create_directory(mode_path)

        # Determine the paths to the 'yes' and 'no' saturation collection directories
        self.yes_path = os.path.join(mode_path, "yes")
        self.no_path = os.path.join(mode_path, "no")

        # Current index of saturation sources
        self.current_index_yes = -1
        self.current_index_no = -1

        for path in fs.files_in_path(self.yes_path, extension="fits"):

            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_yes: self.current_index_yes = index

        for path in fs.files_in_path(self.no_path, extension="fits"):

            name = os.path.basename(path)
            index = int(name.split(".fits")[0])
            if index > self.current_index_no: self.current_index_no = index

    # -----------------------------------------------------------------

    def collect(self):

        """
        This function ...
        :return:
        """

        # Get a list of the filepaths for every FITS file in the current working directory
        file_paths = fs.files_in_path(os.getcwd(), extension="fits", contains=self.config.mode)

        # Keep track of how many files have been processed
        self.number_of_files = len(file_paths)
        self.processed = 0

        # Loop over all FITS files found in the current directory
        for file_path in file_paths:

            # Get information
            name = os.path.basename(file_path).split(".fits")[0]
            info, index = name.split("_")

            # Open the image, select all frames
            image = Image.from_file(file_path, always_call_first_primary=False)
            image.frames.select_all()

            # Create a source
            source = Source.from_image(image)

            self.show(source)

    # -----------------------------------------------------------------

    def show(self, source):

        """
        This function ...
        :param source:
        :return:
        """

        # Create a plot for the source
        source.plot(title="Is this a " + self.description + "? (" + str(self.processed) + " out of " + str(self.number_of_files) + ")", show=False, scale="log")

        # Set current and previous source
        self.previous_source = self.current_source
        self.current_source = source

        # Axes
        axyes = plt.axes([0.6, 0.05, 0.1, 0.075])
        axno = plt.axes([0.7, 0.05, 0.1, 0.075])
        axunsure = plt.axes([0.8, 0.05, 0.1, 0.075])
        axback = plt.axes([0.1, 0.05, 0.1, 0.075])

        # Buttons
        yes_button = Button(axyes, 'Yes')
        yes_button.on_clicked(self.save_yes)
        no_button = Button(axno, 'No')
        no_button.on_clicked(self.save_no)
        unsure_button = Button(axunsure, 'Unsure')
        unsure_button.on_clicked(self.dont_save)
        back_button = Button(axback, 'Back')
        back_button.on_clicked(self.go_back)

        # Show the plot
        plt.show()

        # Increment the counter
        self.processed += 1

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
        self.log.info("Saving the " + self.description + " to " + path)
        self.current_source.saveto(path)

        self.last_path = path

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
        self.log.info("Saving the source to " + path)
        self.current_source.saveto(path)

        self.last_path = path

        # Close the currently active plotting window for this source
        plt.close()

    # -----------------------------------------------------------------

    def show_frame(self, event):

        """
        This function ...
        :param event:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def go_back(self, event):

        """
        This function ...
        :param event:
        :return:
        """

        if self.previous_source is not None:

            # Inform the user
            self.log.info("Going back to the previous source")

            plt.close()

            # Remove saved file
            if self.last_path is not None: fs.remove_file(self.last_path)

            self.processed -= 1

            self.current_source = None
            self.show(self.previous_source)

        else: self.log.warning("Cannot go back")

    # -----------------------------------------------------------------

    def dont_save(self, event):

        """
        This function ...
        :param event:
        :return:
        """

        # Inform the user
        self.log.info("Ignoring source")

        self.last_path = None

        # Close the currently active plotting window for this source
        plt.close()

# -----------------------------------------------------------------
