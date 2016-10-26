#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.datasetcreator Contains the DataSetCreator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .component import DataComponent
from ...magic.core.dataset import DataSet
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class DataSetCreator(DataComponent):

    """
    This class ...
    """
    
    def __init__(self, config=None):
    
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(DataSetCreator, self).__init__(config)

        # -- Attributes --

        # Create a new data set
        self.set = DataSet()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2. Add the images
        self.add_images()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataSetCreator, self).setup()

    # -----------------------------------------------------------------

    def add_images(self):

        """
        This function ...
        :return:
        """

        # Loop over the different image origins
        for path, origin in fs.directories_in_path(self.data_images_path, returns=["path", "name"]):

            # Ignore the Planck data (for now)
            if origin == "Planck": continue

            # Loop over the FITS files in the current directory
            for image_path, image_name in fs.files_in_path(path, extension="fits", not_contains="poisson", returns=["path", "name"]):

                # Open the image frame
                frame = Frame.from_file(image_path)

                # Determine the preparation name
                if frame.filter is not None: prep_name = str(frame.filter)
                else: prep_name = image_name

                # Add the image path
                self.set.add_path(prep_name, image_path)

                # Determine path to poisson error map
                poisson_path = fs.join(path, image_name + "_poisson.fits")

                # Set the path to the poisson error map
                if fs.is_file(poisson_path): self.set.add_error_path(prep_name, poisson_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the dataset
        self.write_dataset()

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        # Save the dataset
        self.set.saveto(self.initial_dataset_path)

# -----------------------------------------------------------------
