#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.component Contains the DataComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class DataComponent(ModelingComponent):
    
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
        super(DataComponent, self).__init__(config)

        # -- Attributes --

        # The path to the galaxy info file
        self.galaxy_info_path = None

        # Different origins
        self.data_origins = ["GALEX", "SDSS", "Halpha", "2MASS", "Spitzer", "WISE", "Herschel", "Planck"]

        # The path to the data/images directory
        self.data_images_path = None

        # The paths to the data/images/ directories for the different origins
        self.data_images_paths = dict()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataComponent, self).setup()

        # Set the path to the galaxy info file
        self.galaxy_info_path = fs.join(self.data_path, "info.dat")

        # Set ...
        self.data_images_path = fs.create_directory_in(self.data_path, "images")

        # Set ...
        for origin in self.data_origins: self.data_images_paths[origin] = fs.create_directory_in(self.data_images_path, origin)

# -----------------------------------------------------------------
