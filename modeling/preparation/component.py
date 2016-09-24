#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.component Contains the PreparationComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import tables
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class PreparationComponent(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PreparationComponent, self).__init__(config)

        # -- Attributes --

        # The names of the different images for the preparation components
        self.prep_names = dict()

        # The paths to the preparation subdirectories for each image
        self.prep_paths = dict()

        # The original paths of the different images
        self.original_paths = dict()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PreparationComponent, self).setup()

        # Check for presence of dataset
        if not fs.is_file(self.initial_dataset_path): raise RuntimeError("Dataset has not been created yet. Run create_dataset first.")

        # Set the image names
        for prep_name in self.initial_dataset.paths:

            # Path
            image_path = self.initial_dataset.paths[prep_name]

            # Name
            image_name = fs.strip_extension(fs.name(image_path))

            self.prep_names[image_name] = prep_name
            self.prep_paths[prep_name] = fs.join(self.prep_path, prep_name)
            self.original_paths[prep_name] = image_path

        # Create the preparation subdirectories for each image
        #fs.create_directories(self.prep_paths.values())

# -----------------------------------------------------------------
