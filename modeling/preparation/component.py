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

        # The path to the preparation info table
        self.prep_info_path = None

        # The names of the different images for the preparation components
        self.prep_names = dict()

        # The original names of the different images, based on the preparation names
        self.original_names = dict()

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

        # Set the path to the preparation info table
        self.prep_info_path = fs.join(self.prep_path, "prep_info.dat")

        # If the table has already been created, load the info
        if fs.is_file(self.prep_info_path):

            # Load the info table
            info = tables.from_file(self.prep_info_path, format="ascii.ecsv")

            # Set the image names
            for i in range(len(info)):

                original_name = info["Image name"][i]
                original_path = info["Image path"][i]
                prep_name = info["Preparation name"][i]

                self.prep_names[original_name] = prep_name
                self.original_names[prep_name] = original_name
                self.prep_paths[prep_name] = fs.join(self.prep_path, prep_name)
                self.original_paths[prep_name] = original_path

        # Create the preparation subdirectories for each image
        #fs.create_directories(self.prep_paths.values())

# -----------------------------------------------------------------
